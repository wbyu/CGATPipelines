"""
PipelineGtfsubset.py - Tasks for GTF subsetting 
===============================================

Reference
---------

"""

import CGAT.Experiment as E
import os
import re
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGATPipelines.Pipeline as P
import pandas as pd


class SubsetGTF():

    def __init__(self, *args, **kwargs):
        pass

    def makeLineDict(self, line):
        D = line.asDict()
        D['chrom'] = line.contig
        D['source'] = line.source
        D['feature'] = line.feature
        D['start'] = line.start
        D['end'] = line.end
        D['score'] = line.score
        D['strand'] = line.strand
        D['frame'] = line.frame

        return D

    def filterCDS(self, infile, outfile, filteroption):
        '''
        '''
        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        with IOTools.openFile(outfile, "w") as outf:
             for line in gtf:
                 D = self.makeLineDict(line)
                 if D[filteroption] == "CDS":
                     outf.write("%s\n" % str(line))

        outf.close()

    def filterExon(self, infile, outfile, filteroption):
        '''
        
        '''
        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        with IOTools.openFile(outfile, "w") as outf:
            for line in gtf:
                D = self.makeLineDict(line)
                if D[filteroption] == "exon":
                    outf.write("%s\n" % str(line))

    def filterExonCoding(self, infile, outfile, filteroptions):
        '''
        
        '''

        feature, gene_biotype = filteroptions

        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        with IOTools.openFile(outfile, "w") as outf:
            for line in gtf:
                D = self.makeLineDict(line)
                if D[feature] == "exon" and D[gene_biotype] == "protein_coding":
                    outf.write("%s\n" % str(line))

    def filterLincExonCoding(self, infile, outfile, filteroptions):
        '''
        
        '''

        feature, gene_biotype = filteroptions

        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        with IOTools.openFile(outfile, "w") as outf:
            for line in gtf:
                D = self.makeLineDict(line)
                if D[feature] == "exon" and D[gene_biotype] == "lincRNA":
                    outf.write("%s\n" % str(line))

    def filterNonCodingExonCoding(self, infile, outfile, filteroptions):
        '''
        
        '''
        feature, gene_biotype = filteroptions

        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        with IOTools.openFile(outfile, "w") as outf:
            for line in gtf:
                D = self.makeLineDict(line)
                if D[feature] == "exon" and not D[gene_biotype] == "protein_coding":
                    outf.write("%s\n" % str(line))
