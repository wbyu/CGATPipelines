####
####
##
## Project PairsDBTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: calculate_histogram_2D.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####


import sys, re, string, os, getopt, time, optparse

import Experiment
import Histogram2D

USAGE = """python calculate_histogram.py < stdin > stdout
read in data and build histogram of column
-1, --column1    column to take [default = 0]
-2, --column2    column to take [default = 1]
# at start of line is a comment
"""

param_column1 = 0
param_column2 = 1
param_bin_size1 = 1
param_bin_size2 = 1
param_titles = True

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: calculate_histogram_2D.py 2782 2009-09-10 11:40:29Z andreas $", usage = USAGE)

    parser.add_option("-t", "--titles", dest="titles", action="store_true",
                      help="input data has title in first row [default=%default]." )

    parser.add_option( "--no-titles", dest="titles", action="store_false",
                      help="input data has no title in first row [default=%default]." )

    parser.add_option("-1", "--column1", dest="column1", type="int",
                      help="first column to use [default=%default]." )
 
    parser.add_option("-2", "--column2", dest="column2", type="int",
                      help="second column to use [default=%default]." )
 
    parser.add_option("--bin-size1", dest="bin_size1", type="float",
                      help="bin size for first column [default=%default]." )

    parser.add_option("--bin-size2", dest="bin_size2", type="float",
                      help="bin size for second column [default=%default]." )

    parser.set_defaults( 
        column1 = 1,
        column2 = 2,
        bin_size1 = 1.0,
        bin_size2 = 1.0,
        titles = True )


    (options, args) = Experiment.Start( parser )                                  
    options.column1 -= 1
    options.column2 -= 1

    histograms = []

    vals = []
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    if options.titles:
        data = string.split(lines[0][:-1], "\t")        
        print "\t".join((data[options.column1], data[options.column2], "counts"))
        del lines[0]

    ninput, noutput, nskipped = 0, 0, 0

    for l in lines:
        ninput += 1
        data = string.split(l[:-1], "\t")

        try:
            val = map( string.atof, (data[options.column1],data[options.column2]))
        except IndexError:
            nskipped += 1
            continue
        except ValueError:
            nskipped += 1
            continue

        vals.append( val )
        noutput += 1

    lines = None
    
    h = Histogram2D.Calculate( vals, bin_function=lambda x:(int(x[0] / options.bin_size1), int(x[1] / options.bin_size2)) )

    Histogram2D.Print(h, bin_function=lambda x: (x[0] * options.bin_size1, x[1] * options.bin_size2, x[2]))

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput,noutput,nskipped) )

    Experiment.Stop()
