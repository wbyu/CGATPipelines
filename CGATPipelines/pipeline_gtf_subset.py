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

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA_primary_transcript"]

    m.filterGFF3(outfile, filteroption, filteritem)


@transform(PARAMS['mirbase_filename_mir_gff'],
           suffix(PARAMS['mirbase_filename_mir_gff']),
           PARAMS['interface_geneset_mir_gff'])
def buildmiRNonPrimaryTranscript(infile, outfile):

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
         loadRepeats)
def full():
    '''build all targets - A dummy task to run the pipeline to
    completion.'''
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
