################################################################################
#   Gene prediction pipeline 
#
#   $Id: nr2table.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2006 Tyler ???? and Andreas Heger 
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
import os, sys, string, re, optparse, math, time, tempfile, subprocess

USAGE="""python %s [OPTIONS] < nr.fasta

convert a nr fasta file into a more readable list of fields.
""" % sys.argv[0]

import Experiment

import Bio, Bio.Fasta

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: nr2table.py 2782 2009-09-10 11:40:29Z andreas $", usage = USAGE)

    parser.set_defaults()

    (options, args) = Experiment.Start( parser )
    
    parser = Bio.Fasta.RecordParser()
    iterator = Bio.Fasta.Iterator( sys.stdin, parser)

    sequences = []
    
    ninput, noutput, nentries = 0, 0, 0

    options.stdout.write( "gid\tsrc\tacc\tprotid\tannotation\tspecies\n" )
    while 1:
        
        cur_record = iterator.next()

        if cur_record == None: break

        ninput += 1

        records = cur_record.title.split( chr(1) )
        for record in records:

            a, anno = re.search("(\S+)\s+(.+)", record).groups()

            vals = a.split("|")
            
            try:
                gi, gid, src, acc, protid = vals
            except ValueError:
                raise "parsing error for record: %s" % record
            
            try:
                annotation, species = re.search( "(.+)\s+\[(.*)\]", anno).groups()
            except AttributeError:
                annotation = anno.strip()
                species = ""

            annotation = re.sub( "\t", " ", annotation)

            options.stdout.write( "\t".join( (gid, src, acc, protid,
                                              annotation, species) ) + "\n" )

            nentries += 1
            
        noutput += 1
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nentries=%i\n" % (ninput, noutput, nentries ))
        
    Experiment.Stop()
    
    