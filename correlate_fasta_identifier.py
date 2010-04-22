################################################################################
#   Gene prediction pipeline 
#
#   $Id: correlate_fasta_identifier.py 14 2005-08-09 15:24:07Z andreas $
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
import os, sys, string, re, getopt, tempfile

USAGE="""python %s [OPTIONS] filename < stdin > stdout

Version: $Id: correlate_fasta_identifier.py 14 2005-08-09 15:24:07Z andreas $

Given a two files fasta sequences substitute
identifiers in stream with those given in file with filename.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
"""

param_loglevel = 0

param_long_options=["verbose=", "help"]
param_short_options="v:h"

##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)
        
    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)

    if len(args) != 1:
        print "please supply filename with replacement identifiers."
        print USAGE
        sys.exit(1)


    identifiers = []
    infile = open(args[0], "r")
    for line in infile:
        if line[0] == ">":
            identifiers.append( line[1:string.find(" ", line)] )
    infile.close()

    x = 0
    for line in sys.stdin:
        if line[0] == ">":
            if x >= len(identifiers): raise "different number of sequences."            
            line = ">" + identifiers[x] + "\n"
            x += 1
        print line[:-1]
    
            