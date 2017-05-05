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
import os
import sqlite3
from ruffus import follows, transform, merge, mkdir, files, jobs_limit,\
    suffix, regex, add_inputs
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineGtfsubset as PipelineGtfsubset
import CGATPipelines.PipelineUCSC as PipelineUCSC

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


def connectToUCSC():
    return PipelineGtfsubset.connectToUCSC(
        host=PARAMS["ucsc_host"],
        user=PARAMS["ucsc_user"],
        database=PARAMS["ucsc_database"])
# -----------------------------------------------------------------
# ENSEMBL gene set

#### check this code because there may be a problem using set-genebiotype_to source

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


@files(PARAMS['interface_geneset_ucsc_genome_gtf'],
       PARAMS['interface_geneset_cds_gtf'])
def buildCdsTranscript(infile, outfile):
    '''
    Output the CDS features from an ENSEMBL gene set

    takes all of the features from a :term:`gtf` file
    that are feature types of ``CDS``.

    Note - we have not filtered on gene_biotype because some of the CDS
    are classified as polymorphic_pseudogene.
        

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["CDS"]

    m.filterGTF(outfile, filteroption, filteritem, operators=None)


@files(PARAMS['interface_geneset_ucsc_genome_gtf'],
       PARAMS['interface_geneset_exons_gtf'])
def buildExonTranscript(infile, outfile):
    '''
    Output of the exon features from abn ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``exon``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["exon"]

    m.filterGTF(outfile, filteroption, filteritem, operators=None)

@files(PARAMS['interface_geneset_ucsc_genome_gtf'],
       PARAMS['interface_geneset_coding_exons_gtf'])
def buildCodingExonTranscript(infile, outfile):
    '''
    Output of the coding exon features from abn ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``exon``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = [PARAMS['ensembl_cgat_feature'],
                    PARAMS['ensembl_cgat_gene_biotype']]
    filteritem = ["exon", "protein_coding"]

    m.filterGTF(outfile, filteroption, filteritem, operators="and")


@files(PARAMS['interface_geneset_ucsc_genome_gtf'],
       PARAMS['interface_geneset_lincrna_exons_gtf'])
def buildLincRNAExonTranscript(infile, outfile):
    '''
    Output of the lincRNA features from an ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``lincRNA``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroptions = [PARAMS['ensembl_cgat_feature'],
                     PARAMS['ensembl_cgat_gene_biotype']]

    filteritem = ["exon", "lincRNA"]

    m.filterGTF(outfile, filteroptions, filteritem, operators="and")


@files(PARAMS['interface_geneset_ucsc_genome_gtf'],
       PARAMS['interface_geneset_noncoding_exons_gtf'])
def buildNonCodingExonTranscript(infile, outfile):
    '''
    Output of the non-coding exon features from an ENSEMBL gene set

    Remove all of the features from a :term:`gtf` file
    that are features of ``exon`` and are protein-coding

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroptions = [PARAMS['ensembl_cgat_feature'],
                     PARAMS['ensembl_cgat_gene_biotype']]
    filteritem = ["exon", "protein_coding"]

    m.filterGTF(outfile, filteroptions, filteritem, operators="and not")



# ---------------------------------------------------------------
# UCSC derived annotations
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_rna_gff"]), ))
def importRNAAnnotationFromUCSC(infile, outfile):
    """This task downloads UCSC repetetive RNA types.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_rnatypes"]),
        outfile=outfile,
        remove_contigs_regex=PARAMS["ensembl_remove_contigs"])


@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_repeats_gff"]), ))
def importRepeatsFromUCSC(infile, outfile):
    """This task downloads UCSC repeats types as identified
    in the configuration file.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_repeattypes"]),
        outfile=outfile)

# ---------------------------------------------------------------
# miRBase annotations

@files(PARAMS['mirbase_filename_mir_gff'],
       PARAMS['interface_geneset_primary_mir_gff'])
def buildmiRPrimaryTranscript(infile, outfile):

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA_primary_transcript"]

    m.filterGFF3(outfile, filteroption, filteritem)


@files(PARAMS['mirbase_filename_mir_gff'],
       PARAMS['interface_geneset_mir_gff'])
def buildmiRNonPrimaryTranscript(infile, outfile):

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA"]

    m.filterGFF3(outfile, filteroption, filteritem)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
