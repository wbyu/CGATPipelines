################################################################################
#   Gene prediction pipeline 
#
#   $Id: csv2xls.py 2782 2009-09-10 11:40:29Z andreas $
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
import os, sys, string, re, getopt, time, optparse, math, tempfile

import Experiment
import csv, CSV

import pyExcelerator 

def GetHeaderStyle():
    fnt = pyExcelerator.Font()
    fnt.name = 'Arial'
    fnt.colour_index = 10
    fnt.outline = False

    borders = pyExcelerator.Borders()
    borders.bottom = 10

    style = pyExcelerator.XFStyle()
    style.font = fnt
    style.borders = borders
    return style


def GetDataStyle():
    fnt = pyExcelerator.Font()
    fnt.name = 'Arial'
    fnt.colour_index = 0
    fnt.outline = False

    style = pyExcelerator.XFStyle()
    style.font = fnt
    return style


if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: csv2xls.py 2782 2009-09-10 11:40:29Z andreas $")
    
    parser.add_option( "-o", "--outfile=", dest="output_filename", type="string",
                       help="write to output filename." )

    parser.set_defaults(
        output_filename=None,
        )

    (options, args) = Experiment.Start( parser, add_csv_options  = True)

    if not options.output_filename:
        raise "please specify an output filename."

    w = pyExcelerator.Workbook()

    ## create styles
    header_style = GetHeaderStyle()
    data_style = GetDataStyle()

    for filename in args:

        lines = filter( lambda x: x[0] != "#", open(filename, "r").readlines())

        if len(lines) == 0: continue

        if options.loglevel >= 2:
            print "# read %i rows" % len(lines)
            sys.stdout.flush()

        headers = lines[0][:-1].split("\t")

        ws = w.add_sheet(os.path.basename(filename))
        
        cur_row = 0
        
        for x in range(len(headers)):
            ws.write( cur_row, x, headers[x], header_style )
            
        cur_row += 1
        
        reader = csv.DictReader( lines, dialect=options.csv_dialect )

        for row in reader:
            row = CSV.ConvertDictionary( row )

            for x in range(len(headers)):
                if headers[x] in row:
                    ws.write( cur_row, x, row[headers[x]], data_style )
                
            cur_row += 1

    
    w.save(options.output_filename)
        
    Experiment.Stop()

