################################################################################
#   Gene prediction pipeline 
#
#   $Id: predictions2disruptions.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, optparse

USAGE="""python %s < predictions > genes

Version: $Id: predictions2disruptions.py 2781 2009-09-10 11:33:14Z andreas $

calculate disruptions in predictions

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-g, --genome-pattern=           pattern for filenames with the genomic DNA (FASTA).
""" % sys.argv[0]

import Experiment
import Genomics
import IndexedFasta
import PredictionParser

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: predictions2disruptions.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome pattern."  )
    
    parser.add_option( "--start-codon-boundary", dest="start_codon_boundary", type="int",
                      help="maximum extension for start codon (make divisible by 3)."  )
    
    parser.add_option( "--stop-codon-boundary", dest="stop_codon_boundary", type="int",
                      help="maximum extension for stop codon (make divisible by 3)."  )

    
    parser.set_defaults(
        genome_file = "genome.fasta",
        stop_codons = ("TAG", "TAA", "TGA")
        )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    p = PredictionParser.PredictionParserEntry()

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    
    for line in sys.stdin:
        
        if line[0] == "#": continue

        p.Read(line)

        genomic_sequence = fasta.getSequence( p.mSbjctToken, p.mSbjctStrand,
                                              p.mSbjctGenomeFrom, p.mSbjctGenomeTo )
        
        if options.loglevel >= 2:
            options.stdlog.write ("# parsing alignment %s\n" % p.mAlignmentString)
        try:
            nintrons, nframeshifts, ngaps, nsplits, nstopcodons, disruptions =\
                      Genomics.CountGeneFeatures( 0,
                                                  p.mMapPeptide2Genome,
                                                  genomic_sequence,
                                                  border_stop_codon = 0,
                                                  stop_codons = options.stop_codons )
        except ValueError, msg:
            options.stderr.write( "# parsing error: %s in line %s\n" % (line[:-1], msg))
            sys.exit(1)

        for type, \
                cds_pos_from, cds_pos_to, \
                genome_pos_from, genome_pos_to in disruptions:
            options.stdout.write( "\t".join(map(str, (p.mPredictionId,
                                                      type,
                                                      cds_pos_from, cds_pos_to,
                                                      genome_pos_from + p.mSbjctGenomeFrom,
                                                      genome_pos_to + p.mSbjctGenomeFrom) ) )+ "\n")

        options.stdout.flush()
        
    Experiment.Stop()