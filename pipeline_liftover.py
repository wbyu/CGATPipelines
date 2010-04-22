################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_liftover.py 2900 2010-04-13 14:38:00Z andreas $
#
#   Copyright (C) 2010 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""

:Author: Andreas Heger
:Release: $Id: pipeline_liftover.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Pipeline for mapping a set of intervals onto a reference or target genome.

.. todo::
   * make the merging step optional. Currently overlapping intervals are merged.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import csv, gzip
from ruffus import *
import sqlite3

import IOTools
import Experiment as E
import Pipeline as P
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools, Database, GFF, GTF
import Database

import PipelineGeneset as PGeneset
import PipelineAnnotator as PAnnotator

if not os.path.exists("conf.py"):
    raise IOError( "could not find configuration file conf.py" )

execfile("conf.py")

PARAMS = P.getParameters()

PARAMS.update( {
    "annotation" : "regions.gff",
    "genes": "genes.gtf",
    "transcripts": "transcripts.gtf.gz" } )

###################################################################
@transform(  [ x for x in glob.glob('*.gtf.gz') if not x.endswith(".mapped.gtf.gz") ]
             , suffix(".gtf.gz")
             , '.psl.gz' )
def convertGtf2Psl( infile, outfile ):
    """convert a gtf to a psl file."""
    
    track = outfile[:-len(".psl.gz")]
    genomefile = "%s_genome" % track
    
    statement = """gunzip 
    < %(infile)s 
    | awk '$3 == "exon"' 
    | python %(scriptsdir)s/gff2blat.py 
           --is-gtf 
           --log=%(outfile)s.log 
    | gzip > %(outfile)s
    """ 
    P.run()

###################################################################
@transform(  '*.bed.gz', 
             suffix(".bed.gz"), 
             '.psl.gz' )
def convertBed2Psl( infile, outfile ):
    """convert a bed to a psl file."""
    
    track = outfile[:-len(".bed.gz")]
    genomefile = PARAMS["%s_genome" % track]
    
    statement = """gunzip < %(infile)s 
    | python %(scriptsdir)s/bed2psl.py 
         --genome=%(genomefile)s
         --log=%(outfile)s.log 
    | gzip > %(outfile)s
    """ 
    P.run()

###################################################################
###################################################################
###################################################################
@transform( (convertGtf2Psl, convertBed2Psl)
            , suffix(".psl.gz")
            , '.transcripts' )
def mergeTranscripts( infile, outfile ):
    """merge transcripts before mapping.

    Overlapping transcripts are combined in order to
    speed up the mapping process.
    """

    track = outfile[:-len(".transcripts")]
    genomefile = PARAMS["%s_genome" % track]
    
    statement = """
        gunzip < %(infile)s 
        | awk '/psLayout/ { x = 4; next; } x > 0 { --x; next} { print; }' 
        | sort -k 14,14 -k 16,16n
	| %(cmd-farm)s 
		--split-at-column=14 
		--output-header 
		--renumber="%%06i" 
		--renumber-column=":id" 
		--log=%(outfile)s.log 
		--subdirs \
        "python %(scriptsdir)s/blat2assembly.py 
               --staggered=all 
               --method=region
               --method=transcript
                --threshold-merge-distance=0 
                --threshold-merge-overlap=3
		--genome=%(genomefile)s
		--mali-output-format=fasta 
		--log=%(outfile)s.log 
		--output-filename-pattern=%%DIR%%%(outfile)s.%%s" 
	> %(outfile)s"""

    P.run()

@transform( mergeTranscripts
            , suffix(".transcripts" )
            , '.merged.mapped.psl' )
def mapMergedTranscripts( infile, outfile ):
    """map transcripts from PSL.

    Mapping from PSL is equivalent to first converting to genePred format 
    and using the option -gp.
    """

    track = outfile[:-len(".merged.mapped.psl")]
    chainfile = PARAMS["%s_chain" % track]

    statement = """
        liftOver -minMatch=0.2 -minBlocks=0.01 -pslT 
                 %(infile)s.transcripts.psl 
                 <(gunzip < %(chainfile)s) 
                 %(outfile)s 
                 %(outfile)s.unmapped 
        >& %(outfile)s.log
        """
    P.run()

@transform( (convertGtf2Psl, convertBed2Psl)
            , suffix(".psl.gz" )
            , '.mapped.psl' )
def mapTranscripts( infile, outfile ):
    """map transcripts from PSL.

    Mapping from PSL is equivalent to first converting to genePred format 
    and using the option -gp.
    """

    track = outfile[:-len(".mapped.psl")]
    chainfile = PARAMS["%s_chain" % track]

    if not os.path.exists( chainfile ):
        raise OSError("chain %s not found" % chainfile)

    statement = """
        liftOver -minMatch=0.2 -minBlocks=0.01 -pslT 
                 <(gunzip < %(infile)s )
                 <(gunzip < %(chainfile)s) 
                 %(outfile)s 
                 %(outfile)s.unmapped 
        >& %(outfile)s.log
        """
    P.run()


@transform(  (mapMergedTranscripts, mapTranscripts)
             , suffix(".psl")
             , '.gtf.gz' )
def convertMappedPslToGtf( infile, outfile ):
    '''convert to gtf for export.'''
    statement = """
    python %(scriptsdir)s/blat2gff.py --as-gtf 
    < %(infile)s 
    | gzip
    > %(outfile)s
    """
    P.run()

@transform( convertMappedPslToGtf
            , suffix(".gtf.gz")
            , '.summary' )
def summary( infile, outfile ):
    '''compute mapping stats.'''

    def _getfiles( filename ):
    
        track = outfile[:-len(".mapped.summary")]
        if track.endswith( ".merged" ):
            xtrack = track[:-len(".merged")]
            finput = "%s.psl.gz" % xtrack
            fmerged = "%s.transcripts.transcripts.psl" % xtrack 
            fmapped = "%s.mapped.psl" % track
        else:
            finput = "%s.psl.gz" % track
            fmerged = finput
            fmapped = "%s.mapped.psl" % track
        return track, finput, fmerged, fmapped
    
    outf = open(outfile, "w")
    outf.write( "track\tinput\tmerged\tpmerged\tmapped\tpmapped\tpoutput\n" )
    
    def countPSL( filename ):
        if filename.endswith(".gz"):
            i = gzip.open( filename )
        else:
            i = open(filename)
        ll = [ x[:10] for x in i.readlines() if not x.startswith("#") ]
        if ll[0].startswith("psLayout"): return len(ll) - 5
        else: return len( ll)
    
    track, finput, fmerged, fmapped = _getfiles( outfile )
    ninput = countPSL( finput )
    # subtract header
    nmerged = countPSL( fmerged ) - 5
    nmapped = countPSL( fmapped )
    
    outf.write( "%s\t%i\t%i\t%s\t%i\t%s\t%s\n" % 
                ( track,
                  ninput,
                  nmerged,
                  IOTools.prettyPercent( nmerged, ninput),
                  nmapped,
                  IOTools.prettyPercent( nmapped, nmerged),
                  IOTools.prettyPercent( nmapped, ninput) ) )
                  

@follows( convertMappedPslToGtf, summary )
def full(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

    