"""
===================
Gtf_subset pipeline
===================

This pipeline generates a number of GTF files that can be used with downstream
CGAT pipelines. The user will download a GTF from ENSEMBL and then the GTF is
parsed and filtered. 

====
Code
====
"""
import sys
import re
import os
import glob
import sqlite3
from ruffus import follows, transform, merge, mkdir, files, jobs_limit,\
    suffix, regex, add_inputs

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGATPipelines.PipelineGtfsubset as PipelineGtfsubset

###################################################
# Pipeline configuration
###################################################
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


def connect():
    '''connect to database.'''

    dbh = sqlite3.connect(PARAMS["database_name"])
    return dbh


# -----------------------------------------------------------------
# ENSEMBL gene set


@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"], PARAMS['interface_geneset_all_gtf'])
def buildUCSCGeneSet(infile, outfile):
    '''output sanitized ENSEMBL geneset.

    This method outputs an ENSEMBL gene set after some sanitizing steps:

    1. Chromosome names are changed to the UCSC convention.
    2. Chromosomes that match the regular expression specified in
       the configuration file are removed.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       geneset in :term:`gtf` format.

    '''
    statement = ['''zcat %(infile)s
    | grep 'transcript_id'
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=ucsc
    --log=%(outfile)s.log
    ''']

    if PARAMS["ensembl_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --contig-pattern="%(ensembl_remove_contigs)s" ''')

    statement.append(
        '''
        | cgat gtf2gtf
        --method=set-gene_biotype-to-source
        --log=%(outfile)s.log
        | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()

@follows(buildUCSCGeneSet)
@files(PARAMS["ensembl_filename_gtf"],
       PARAMS['interface_geneset_ucsc_genome_gtf'])
def buildGenomeGeneSet(infile, outfile):
    '''output sanitized ENSEMBL geneset and removal of all contigs
    that are not present in the downloaded UCSC genome.

    This method outputs an ENSEMBL gene set after some sanitizing steps:

    1. Chromosome names are changed to the UCSC convention.
    2. Transcripts that are not part of the chosen genome assembly
       are removed.
    3. Chromosomes that match the regular expression specified in
       the configuration file are removed.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       geneset in :term:`gtf` format.

    '''

    statement = ['''zcat %(infile)s
    | grep 'transcript_id'
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    ''']

    if PARAMS["ensembl_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --contig-pattern="%(ensembl_remove_contigs)s" ''')

    statement.append(
        '''
        | cgat gtf2gtf
        --method=set-gene_biotype-to-source
        --log=%(outfile)s.log
        | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()


@files(PARAMS['interface_geneset_ucsc_genome_gtf'], PARAMS['interface_geneset_cds_gtf'])
def buildCdsGeneSet(infile, outfile):
    m = PipelineGtfsubset.SubsetGTF()

    m.dataframe2Gtf(infile, outfile)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
