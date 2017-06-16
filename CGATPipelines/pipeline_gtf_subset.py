"""
===================
Gtf_subset pipeline
===================

Overview
========

This pipeline generates a number of annotations that can be used with
downstream CGAT pipelines. The user will download a GTF from ENSEMBL
and then the GTF is parsed and filtered.In addition to downloading an
ensembl GTF the user will need to download an assembly report for their
specific genome and add it to the directory the pipeline is ran.

Common to all of the annotations generated in this pipeline is that they
are genomic - i.e. they are genomic intervals or relate to genomic intervals.
Thus, annotations are tied to a particular version of the genome. This is
parameterised within the pipeline.ini configuration file. The pipeline
follows two principle releases: the UCSC_ genome assembly and an ENSEMBL_
geneset version.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Principle targets
-----------------


full
   This will run the entire pipeline


Configuration
-------------

The :file:`pipeline.ini` needs to be edited so that it points to the
appropriate locations of the auxiliary files.


On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Input
-----

This pipeline requires a Ensembl GTF, mirBase GFF3 file and an assembly report.

Ensembl GTF:
    This can be downloaded from http://www.ensembl.org/info/data/ftp/index.html
    Note: CGAT pipelines use the UCSC GTF convention (chr naming of contigs)
    and therefore the GTF is sanitized to the UCSC convention. As part of this
    process an ncbi assembly report needs to be specified (see below).

Assembly report:
    This is downloaded from the ncbi assembly page for your specific genome.
    Using hg19 as an example:
        Navigate to www......
        From the database tab select assembly and add your genome into the
        search bar i.e. hg19.
        Then click the link "Download the full sequence report"
        Add it to the folder where the pipeline will be ran, the file is
        for hg38 is called "GRCh38.p10_assembly_report.txt".

miRbase GFF3:
   This can be downloaded from miRbase http://www.mirbase.org/ftp.shtml.
   A path to the :term:`GFF3` file needs to be specified in the pipelin.ini
   configuration file. Make sure that the genome build version of the GFF3
   annotation file matches the ENSEMBL genome.

Running
-------

The pipeline can be run as any other CGAT pipeline, but as its purpose
is to provide a set of annotation that can be used by other pipelines
therefore there is an etiquette to be followed:

Using the pipeline results
--------------------------

The gtf_subset pipeline provides an interface for presenting its
results to other pipelines. The interface is defined in the file
:file:`pipeline.ini`. For example::

   [interface]
   # fasta file with cdna sequences
   cdna_fasta=ensembl.dir/cdna.fasta

The ini file of pipeline annotations can be loaded into the parameter
dictionary of your own pipeline::

    PARAMS.update(P.peekParameters(
         PARAMS["annotations_dir"],
         "pipeline_annotations.py",
         prefix="annotations_"),
         update_interface=True)

Parameters from the gtf_subset pipeline are now accessible via the
``annotations_`` prefix. As a result, the file
:file:`ensembl.dir/cdna.fasta` can be accessed as::

    PARAMS['annotations_cdna_fasta']


Working with non-ENSEMBL species
--------------------------------

:doc:`pipeline_gtf_subset` is very much wedded to annotations in ENSEMBL-
and UCSC_. Using a non-ENSEMBL species or non-UCSC species is possible by
building ENSEMBL- or UCSC-like input files. Even so, annotations that are
downloaded from the ENSEMBL or UCSC database will not be built. You will
thus need to ask if it is worth the effort.

As other pipelines will depend on the annotations in this pipeline it is
necessary to set up a :doc:`pipeline_gtf_subset` stub. To do so, simply
build the config files by running::

   python <SRC>pipeline_annotations.py config

and create the files that are being used in the downstream pipeline
explicitely (for example, for protein coding genes)::

   mkdir ensembl.dir
   cp <MYDATADIR>/my_gtf_geneset.gtf.gz ensembl.dir/geneset_coding.gtf.gz


Pipeline output
===============

The results of the computation are all stored in an sqlite relational
database file csvdb or as compressed files in genomic formats in the pipeline
directories. Output files are grouped by sections listed below.

The sections correspond to primary targets in the pipeline, i.e., to
build all annotations in the section ``assembly`` type::

   python <SRC>pipeline_annotations.py make assembly


section: ensembl
----------------

geneset_all.gtf.gz
   The full gene set after reconciling with assembly. Chromosomes names are
   renamed to be consistent with the assembly and some chromosomes
   are optionally removed. This file is the starting point for
   all annotations derived from the ENSEMBL geneset.

geneset_cds.gtf.gz
   A :term:`gtf` formatted file with only the CDS parts of transcripts.
   This set will naturally include only coding transcripts. UTR regions
   have been removed.

geneset_exons.gtf.gz
   A :term:`gtf` formatted file with only the exon parts of transcripts.
   This set includes both coding and non-coding transcripts. Coding
   transcripts span both the UTR and the CDS.

geneset_coding_exons.gtf.gz
   :term:`gtf` file with exon parts of protein coding transcripts.
   All other features are removed. These are all features annotated
   as "protein_coding" in the ENSEMBL gtf file.

geneset_noncoding_exons.gtf.gz
   :term:`gtf` file with exon parts of non-coding transcripts
   all other features are removed. These are all transcripts not
   annotated as "protein_coding" in the ENSEMBL gtf file.

geneset_lincrna_exons.gtf.gz
   :term:`gtf` file with exon parts of lincRNA transcripts. These
   are transcripts annotated as "lincRNA" in the ENSEMBL gtf file.

geneset_flat.gtf.gz
   A :term:`gtf` formatted file of flattened gene
   models. All overlapping transcripts have been merged. This set
   includes both coding and non-coding transcripts.

geneset_introns.gtf.gz
   A :term:`gtf` formatted file containing all intron features. All
   protein coding genes are retained and their exonic sequences are
   removed to retain introns from nested genes that may overlap.

section: mirbase
----------------

miRNA_non_primary_transcripts.gff3.gz
   A :term:`gff3` formatted file containing all of the non primary miRNA
   transcripts from mirbase

miRNA_primary_transcripts.gff3.gz
   A :term:`GFF3` formatted file containing all of the primery miRNA
   transcripts from miRbase.

section: ucsc
-------------

repeats.gff.gz
   :term:`gff` formatted file with structural/complex repeats

rna.gff.gz
   :term:`gff` formatted file with ribosomal rna annotations

section: geneset
----------------
Annotations derived from the ENSEMBL gene set. Annotations in
this section have been computed from the ENSEMBL gene set.
Results are in the directory :file:`geneset.dir`.

ref_flat.txt
   This creates a flat reference file from geneset_flat.gtf.gz
   for use in picard tools RNAseqmetrics.

section: bed
------------

This directory contains bed files that are generated from other annotations
in this pipeline.

genomic_context.bed.gz
   bed-formatted file with genomic context

Database design
---------------

Tables in the database usually represent genomic features such as
transcripts, genes or chromosomes. These are identified by the
following columns:

+--------------------+-----------------------------------------+
|*Column*            |*Content*                                |
+--------------------+-----------------------------------------+
|transcript_id       |ENSEMBL transcript identifier            |
+--------------------+-----------------------------------------+
|gene_id             |ENSEMBL gene id                          |
+--------------------+-----------------------------------------+
|contig              |Chromosome name                          |
+--------------------+-----------------------------------------+

Example
=======

**Supply example data**


====
Code
====
"""
import sys
import os
import sqlite3
from ruffus import follows, transform, merge, mkdir, files, jobs_limit,\
    suffix, regex, add_inputs, originate
import CGAT.IOTools as IOTools
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

###################################################################
# ENSEMBL gene set
###################################################################


@follows(mkdir('ensembl.dir'))
@transform(PARAMS["ensembl_filename_gtf"],
           regex("(\S+)"),
           r"%s" % PARAMS['interface_geneset_all_gtf'])
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
    --assembly-report="%(ncbi_assembly_report)s"
    --log=%(outfile)s.log
    ''']

    if PARAMS["ncbi_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --contig-pattern="%(ncbi_remove_contigs)s" ''')

    statement.append(
        '''
        | cgat gtf2gtf
        --method=set-gene_biotype-to-source
        --log=%(outfile)s.log
        | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
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


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_exons_gtf'])
def buildExonTranscript(infile, outfile):
    '''
    Output of the exon features from an ENSEMBL gene set

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


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
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


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
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


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
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


@transform((buildUCSCGeneSet,
            buildCdsTranscript,
            buildExonTranscript,
            buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           suffix(".gtf.gz"), "_gtf.load")
def loadTranscripts(infile, outfile):
    '''load transcripts from a GTF file into the database.

    The table will be indexed on ``gene_id`` and ``transcript_id``

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id "
        "--allow-empty-file ")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2tsv
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


@P.add_doc(PipelineGtfsubset.buildFlatGeneSet)
@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_flat_gtf'])
def buildFlatGeneSet(infile, outfile):
    PipelineGtfsubset.buildFlatGeneSet(infile, outfile)


@follows(mkdir("geneset.dir"))
@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS["interface_ref_flat"])
def buildRefFlat(infile, outfile):
    '''build flat geneset for Picard RnaSeqMetrics.
    '''

    tmpflat = P.getTempFilename(".")

    statement = '''
    gtfToGenePred -genePredExt -geneNameAsName2 %(infile)s %(tmpflat)s;
    paste <(cut -f 12 %(tmpflat)s) <(cut -f 1-10 %(tmpflat)s)
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmpflat)


@P.add_doc(PipelineGtfsubset.loadGeneInformation)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(mkdir('ensembl.dir'))
@transform(PARAMS["ensembl_filename_gtf"],
           suffix(PARAMS["ensembl_filename_gtf"]),
           "ensembl.dir/gene_info.load")
def loadGeneInformation(infile, outfile):
    '''load the transcript set.'''
    PipelineGtfsubset.loadGeneInformation(infile, outfile)


@originate("protein_coding_gene_ids.tsv")
def identifyProteinCodingGenes(outfile):
    '''Output a list of proteing coding gene identifiers

    Identify protein coding genes from the annotation database table
    and output the gene identifiers

    Parameters
    ----------
    oufile : str
       Output file of :term:`gtf` format
    annotations_interface_table_gene_info : str
       :term:`PARAMS`. Database table name for gene information

    '''

    dbh = connect()

    select = dbh.execute("""SELECT DISTINCT gene_id
    FROM gene_info
    WHERE gene_biotype = 'protein_coding'""" % locals())

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("gene_id\n")
        outf.write("\n".join((x[0] for x in select)) + "\n")

# check this function output to make sure it is outputting introns- I changed some inputs


@follows(mkdir("geneset.dir"))
@transform(buildFlatGeneSet,
           regex(".*"),
           add_inputs(identifyProteinCodingGenes,
                      buildExonTranscript),
           PARAMS['interface_geneset_intron_gtf'])
def buildIntronGeneModels(infiles, outfile):
    '''build protein-coding intron-transcipts

    Retain the protein coding genes from the input gene set and
    convert the exonic sequences to intronic sequences. 10 bp is
    truncated on either end of an intron and need to have a minimum
    length of 100. Introns from nested genes might overlap, but all
    exons are removed.

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`gtf` format
    infiles[1] : str
       Input filename in :term:`tsv` format

    outfile: str
       Output filename in :term:`gtf` format

    annotations_interface_geneset_exons_gtf: str, PARAMS
       Filename for :term:`gtf` format file containing gene set exons

    '''

    infile, genes_tsv, filename_exons = infiles

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter
    --map-tsv-file=%(genes_tsv)s
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort
    --sort-order=gene
    | cgat gtf2gtf
    --method=exons2introns
    --intron-min-length=100
    --intron-border=10
    --log=%(outfile)s.log
    | cgat gff2gff
    --method=crop
    --crop-gff-file=%(filename_exons)s
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=set-transcript-to-gene
    --log=%(outfile)s.log
    | awk -v OFS="\\t" -v FS="\\t" '{$3="exon"; print}'
    | gzip
    > %(outfile)s
    '''
    P.run()

# Next need to add identifyProteinCodingGenes, buildIntronGeneModels
# aim is to generate the intron gtf here for use in bamstats

################################################################
# UCSC derived annotations
################################################################


@follows(mkdir('ucsc.dir'))
@originate(PARAMS["interface_rna_gff"])
def importRNAAnnotationFromUCSC(outfile):
    """This task downloads UCSC repetetive RNA types.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_rnatypes"]),
        outfile=outfile,
        remove_contigs_regex=PARAMS["ncbi_remove_contigs"])


@follows(mkdir('ucsc.dir'))
@originate(PARAMS["interface_repeats_gff"])
def importRepeatsFromUCSC(outfile):
    """This task downloads UCSC repeats types as identified
    in the configuration file.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_repeattypes"]),
        outfile=outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((importRepeatsFromUCSC,
            importRNAAnnotationFromUCSC),
           suffix(".gff.gz"), "_gff.load")
def loadRepeats(infile, outfile):
    """load genomic locations of repeats into database.

    This method loads the genomic coordinates (contig, start, end)
    and the repeat name into the database.

    Arguments
    ---------
    infile : string
        Input filename in :term:`gff` with repeat annotations.
    outfile : string
        Output filename with logging information. The table name is
        derived from outfile.

    """
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=class "
        "--header-names=contig,start,stop,class")

    statement = """zcat %(infile)s
    | cgat gff2bed --set-name=class
    | grep -v "#"
    | cut -f1,2,3,4
    | %(load_statement)s
    > %(outfile)s"""
    P.run()


# ---------------------------------------------------------------
# miRBase annotations


@transform(PARAMS['mirbase_filename_mir_gff'],
           suffix(PARAMS['mirbase_filename_mir_gff']),
           PARAMS['interface_geneset_primary_mir_gff'])
def buildmiRPrimaryTranscript(infile, outfile):

    '''
    This function will subset a miRbase annotation gff3 file.The GFF3
    file can be downloaded from miRbase. Make sure the annotation matches
    the genome build that you are using.

    This function will subset the GFF3 file by selecting annotations that are
    labled "miRNA_primary_transcript"
    '''

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA_primary_transcript"]

    m.filterGFF3(outfile, filteroption, filteritem)


@transform(PARAMS['mirbase_filename_mir_gff'],
           suffix(PARAMS['mirbase_filename_mir_gff']),
           PARAMS['interface_geneset_mir_gff'])
def buildmiRNonPrimaryTranscript(infile, outfile):

    '''
    This function will subset a miRbase annotation gff3 file.The GFF3
    file can be downloaded from miRbase. Make sure the annotation matches
    the genome build that you are using.

    This function will subset the GFF3 file by selecting annotations that are
    labled "miRNA". This will subset all of the non primary transcripts.
    '''

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA"]

    m.filterGFF3(outfile, filteroption, filteritem)


# Need to write this once andreas has sorted out the GTF parsing option in
# pysam
@transform((buildmiRPrimaryTranscript,
            buildmiRNonPrimaryTranscript),
           suffix(".gff3.gz"), "_gff3.load")
def loadmiRNATranscripts(infile, outfile):
    '''load transcripts from a GTF file into the database.

    The table will be indexed on ``gene_id`` and ``transcript_id``

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id "
        "--allow-empty-file ")

    statement = '''
    zcat %(infile)s
    | cgat gtf2tsv
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


##################################################################
# Generation of BED files from GTF files
##################################################################

@P.add_doc(PipelineGtfsubset.buildGenomicContext)
@follows(mkdir('bed.dir'))
#  also add tasks from miRNA
@merge((importRepeatsFromUCSC,
        importRNAAnnotationFromUCSC,
        buildUCSCGeneSet),
       PARAMS["interface_genomic_context_bed"])
def buildGenomicContext(infiles, outfile):
    PipelineGtfsubset.buildGenomicContext(infiles, outfile)


@follows(buildUCSCGeneSet,
         buildCdsTranscript,
         buildExonTranscript,
         buildCodingExonTranscript,
         buildNonCodingExonTranscript,
         buildLincRNAExonTranscript,
         loadTranscripts,
         loadRepeats,
         buildGenomicContext)
def full():
    '''build all targets - A dummy task to run the pipeline to
    completion.'''
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
