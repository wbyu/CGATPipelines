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
    def gtf2dataframe(self, infile):

        # using the cgat GTF iterator, iterate over the GTF file
        gtf = GTF.iterator(IOTools.openFile(infile, "r"))

        row = []
        for line in gtf:
            # access each of the gtf attributes and then store them as
            # dictionary key/value pairs so a pd DataFrame can be made
            D = line.asDict()
            D['chrom'] = line.contig
            D['source'] = line.source
            D['feature'] = line.feature
            D['start'] = line.start
            D['end'] = line.end
            D['score'] = line.score
            D['strand'] = line.strand
            D['frame'] = line.frame
            D['attributes'] = line.attributes
            row.append(D)

        return pd.DataFrame(row)

    def filterDataFrame(self, dataframe):

        return dataframe[dataframe['feature'] == 'CDS']

    def reorderDataFrame(self, dataframe):

        return dataframe[['chrom','source','feature','start','end',
                   'score','strand','frame', 'attributes']]

    def dataframe2Gtf(self, infile, outfile):

        dataframe = self.gtf2dataframe(infile)
        filtering = self.filterDataFrame(dataframe=dataframe)
        reorder = self.reorderDataFrame(filtering)

        reorder.to_csv(outfile, header=False, index=False, compression='gzip')
