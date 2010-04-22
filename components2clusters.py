import sys, os, re, getopt

################################################################################
#   Gene prediction pipeline 
#
#   $Id: components2clusters.py 14 2005-08-09 15:24:07Z andreas $
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
import os, sys, string, re, getopt, tempfile, time, popen2

USAGE="""python %s [OPTIONS] < components > out

Version: $Id: components2clusters.py 14 2005-08-09 15:24:07Z andreas $


Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-e, --equivalences=             input: equivalences between transcript, gene, and cluster
-c, --components=               input: components
-q, --queries=                  output: map_cluster2queries
-t, --transcripts=              output: map_transcripts2cluster
-f, --filter=                   only write transcripts in filter
""" % sys.argv[0]

param_long_options=["verbose=", "help", 
                    "equivalences=", "components=", "queries=", "transcripts=", "filter="]
param_short_options="v:he:c:q:t:f:"

param_filename_equivalences = None
param_filename_components = None
param_loglevel = 1
param_filename_map_cluster2queries = None
param_filename_map_transcript2cluster = None
param_filename_filter = None

import Experiment

if __name__ == "__main__":

    print Experiment.GetHeader()
    print Experiment.GetParams()    

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
        elif o in ("-e", "--equivalences"):
            param_filename_equivalences = a
        elif o in ("-c", "--components"):
            param_filename_components = a
        elif o in ("-q", "--queries"):
            param_filename_map_cluster2queries = a
        elif o in ("-t", "--transcripts"):
            param_filename_map_transcript2cluster = a
        elif o in ("-f", "--filter"):
            param_filename_filter = a

    filter_transcripts = {}
    if param_filename_filter:
        for line in open(param_filename_filter, "r"):
            if line[0] == "#": continue
            filter_transcripts[line[:-1].split("\t")[0]] = 1
        
    map_gid2tid = {}
    map_tid2cid = {}
    for line in open(param_filename_equivalences, "r"):
        if line[0] == "#": continue
        tid, gid, cid = line[:-1].split("\t")
        map_tid2cid[tid] = cid
        if gid not in map_gid2tid: map_gid2tid[gid] = []
        map_gid2tid[gid].append( tid )

    if param_loglevel >= 1:
        print "# read %i transcripts" % len(map_tid2cid)
    
    map_nr2gid = {}
    
    for line in open(param_filename_components,"r"):
        if line[0] == "#": continue
        if line[:2] == "//": continue
        if string.find( line, "# total number of") != -1: break

        xid, nr = line[:-1].split("\t")

        ## save all gene ids per cluster
        if nr not in map_nr2gid: map_nr2gid[nr] = []
        if xid in map_gid2tid:
            map_nr2gid[nr].append( xid )

    if param_loglevel >= 1:
        print "# read %i clusters" % len(map_nr2gid)

    outfile_clusters = open(param_filename_map_cluster2queries, "w")
    outfile_transcripts = open(param_filename_map_transcript2cluster, "w")    

    ntranscripts = 0
    nclusters = 0
    for nr, gids in map_nr2gid.items():
        queries = {}
        master = None
        for gid in gids:
            for tid in map_gid2tid[gid]:
                if filter_transcripts and tid not in filter_transcripts: continue
                queries[map_tid2cid[tid]] = 1
                if not master: master = map_tid2cid[tid]
                ntranscripts += 1
                outfile_transcripts.write( string.join( map(str,(
                    tid, master )), "\t") + "\n")
        outfile_clusters.write( string.join( map(str, ( \
            master, nr, len(gids), len(queries), string.join( queries.keys(), ";"))), "\t") + "\n")
        nclusters += 1
        
    if param_loglevel >= 1:
        print "# written %i transcripts" % ntranscripts        
        print "# written %i clusters" % nclusters
        
    outfile_clusters.close()
    outfile_transcripts.close()

    print Experiment.GetFooter()