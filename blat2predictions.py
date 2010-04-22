################################################################################
#   Gene prediction pipeline 
#
#   $Id: blat2predictions.py 14 2005-08-09 15:24:07Z andreas $
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
import sys, string, re, getopt

USAGE="""python %s [OPTIONS] < blat > predictions

Convert BLAT output to predictions.

Version: $Id: blat2predictions.py 14 2005-08-09 15:24:07Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-t, --trans                     input is translated DNA
""" % sys.argv[0]


param_long_options=["verbose=", "help", "trans"]
param_short_options="v:ht"

param_trans = None

import Experiment
import PredictionParser

if __name__ == "__main__":

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
        elif o in ( "-t", "--trans"):
            param_trans = 1

    print Experiment.GetHeader()
    print Experiment.GetParams()

    if param_trans:
        parser = PredictionParser.PredictionParserBlatTrans()
    else:
        parser = PredictionParser.PredictionParserBlatCDNA()


    nmatches = 1
    for line in sys.stdin:
        if line[0] == "#": continue
        if not re.match("^[0-9]", line): continue

        try:
            entries = parser.Parse((line,))
        except PredictionParser.AlignmentError, e:
            print "# %s" % str(e)
            print "#", line[:-1]
            sys.exit(1)
        
        for entry in entries:
            entry.mPredictionId = nmatches
            nmatches += 1
            
        print str(entries)

    print Experiment.GetFooter()

    