
"""
======================
Exome Cancer pipeline
======================


.. todo::

   *Final filtering if SNPs/INDELs is currently done in the
   reporting. This should be handled by the pipeline. The SNP output
   would also then be passed to the mutational signature task
   *Document
   *fully make phone home/key option work - GATK public key?  Summarise
   *Indel calling (size of indels called) Example



The exome cancer pipeline imports unmapped reads from matched sample fastqs or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called and filtered


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNPs) on control samples using muTect to generate
      a "panel of normal" variants
   5a. Variant calling (SNPs) with tumour samples using muTect including
      filtering
   5b. Variant calling (indels) using Strelka
   6a. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   6b. Variant annotation with data from eBIO
   6c. Load Network of Cancer Genes (NCG) for Variant annotation in reporting


.. note::

   An optional downsampling analysis can also be performed to assess how
   coverage a control sample affects the called variants

   1. Currently the pipeline is not able to deal with replicates, i.e
      replicates will be treated seperately.



Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <patientID>-<tissue>-<replicate>.<suffix>

``patientID`` and ``tissue`` make up an :term:`experiment`, while ``replicate``
denotes the :term:`replicate` within an :term:`experiment`.
The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted
   by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

If you would like to annotate genes of interest with a particular
value in the results table, create a file call [label]_annotations.tsv
in your working directory listing all the genes. For example, to
annotate all genes identified in a previous shRNA screen, add a file
called shRNA_annoations.tsv listing the genes and the results table
will contain a column called "shRNA" with values "shRNA" and "null".

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+------------+-------------------------------------------+
|*Program*           |*Version*   |*Purpose*                                  |
+--------------------+------------+-------------------------------------------+
|Stampy              |>=0.9.0     |read mapping                               |
+--------------------+------------+-------------------------------------------+
|BWA                 |            |read mapping                               |
+--------------------+------------+-------------------------------------------+
|SAMtools            |            |filtering, SNV / indel calling             |
+--------------------+------------+-------------------------------------------+
|BEDTools            |            |filtering                                  |
+--------------------+------------+-------------------------------------------+
|sra-tools           |            |extracting reads from .sra files           |
+--------------------+------------+-------------------------------------------+
|picard              |>=1.38      |bam/sam files. The .jar files need to be in|
|                    |            |your CLASSPATH environment variable.       |
+--------------------+------------+-------------------------------------------+
|vcf-tools           |            |VCF filtering                              |
+--------------------+------------+-------------------------------------------+
|GATK                | 2.5-2      |local realignment, BQSR, variant calling   |
+--------------------+------------+-------------------------------------------+
|SNPeff              | 3.3        |                                           |
+--------------------+------------+-------------------------------------------+

Pipeline output
===============

The major output is a csvdb containing quality control information
and variant information by patientID and an html report with
similar information.

Example
=======


Code
====

"""

# load modules
from ruffus import *
# from rpy2.robjects import r as R

import numpy
import pandas as pd
import CGAT.Experiment as E
import sys
import os
import sqlite3
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.Pipeline as P
import re
import CGATPipelines.PipelineExome as PipelineExome
import CGATPipelines.PipelineBamStats as PipelineBamStats

USECLUSTER = True

#########################################################################
#########################################################################


def connect():
    '''connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''
    dbh = sqlite3.connect(PARAMS["database_name"])

    return dbh


def getGATKOptions():
    return "-l job_memory=1.4G"


#########################################################################
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False},
    only_import=__name__ != "__main__")

PARAMS = P.PARAMS

PipelineMapping.PARAMS = PARAMS
PipelineMappingQC.PARAMS = PARAMS
PipelineExome.PARAMS = PARAMS
#########################################################################

#######################################################################
# Check for design file to Match Control and Tumor input BAMs ########
#######################################################################

# This section checks for the design table and generates:
# 1. A dictionary, inputD, linking each tumour input file and its matched control,
# as specified in the design table
# 2. A pandas dataframe, df, containing the information from the
#    design table.
# 3. BAM_tumour: a list of tumour (input) bam file names following the naming guidelines
# 4. BAM_control: a list of patient matched control bam files.

# if design table is missing the input bams the df will be empty. This gets
# round the import tests

if os.path.exists("design.tsv"):
    # read the design table
    df = pd.read_csv("design.tsv", sep="\t")

    TUMOUR = list(df['BAM_tumour'].values)
    CONTROL = list(df['BAM_control'].values)

    SAMPLE_DICT = {}
    
    for key, value in zip(TUMOUR, CONTROL):
        SAMPLE_DICT[key] = value
else:
    E.warn("design.tsv is not located within the folder")
    TUMOUR = []
    CONTROL = []

#########################################################################
# Load manual annotations
#########################################################################

@transform("*_annotations.tsv",
           suffix(".tsv"),
           ".load")
def loadManualAnnotations(infile, outfile):

    tmp = P.getTempFilename(".")

    annotation = P.snip(infile, "_annotations.tsv")

    with IOTools.openFile(tmp, "w") as outf:
        outf.write("%s\tgene_id\n" % annotation)
        with IOTools.openFile(infile, "r") as inf:
            for line in inf:
                outf.write("%s\t%s" % (annotation, line))

    P.load(tmp, outfile, options="--add-index=gene_id")
    os.unlink(tmp)

#########################################################################
# Alignment to a reference genome
#########################################################################


@follows(mkdir("bam"))
@transform(("*.fastq.1.gz", "*.fastq.gz", "*.sra"),
           regex(r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
           r"bam/\1.bam")
def mapReads(infile, outfile):
    '''Map reads to the genome using BWA, sort and index BAM file,
    generate alignment statistics and deduplicate using Picard'''

    job_threads = PARAMS["bwa_threads"]
    job_memory = PARAMS["bwa_memory"]

    if PARAMS["bwa_algorithm"] == "aln":
        m = PipelineMapping.BWA(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False)

    elif PARAMS["bwa_algorithm"] == "mem":
        m = PipelineMapping.BWAMEM(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False)
    else:
        raise ValueError("bwa algorithm '%s' not known" % algorithm)

    statement = m.build((infile,), outfile)
    print(statement)
    P.run()

#########################################################################
# Post-alignment QC
#########################################################################

@transform(mapReads,
           regex("bam/(.*).bam$"),
           r"bam/\1.dedup.bam")
def dedup(infile, outfile):
    '''Get duplicate stats from picard MarkDuplicates '''
    PipelineBamStats.buildPicardDuplicateStats(infile, outfile)

@follows(dedup)
@merge("bam/*.bam", "picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineBamStats.loadPicardDuplicateStats(infiles, outfile)

#########################################################################

@transform(dedup,
           regex("bam/(.*).bam$"),
           add_inputs(os.path.join(PARAMS["bwa_index_dir"],
                                   PARAMS["genome"])),
           r"bam/\1.picard_stats")
def buildPicardAlignStats(infiles, outfile):
    ''' build Picard alignment stats '''
    infile, reffile = infiles

    PipelineBamStats.buildPicardAlignmentStats(infile,
                                               outfile,
                                               reffile)


@follows(buildPicardAlignStats)
@merge("bam/*.picard_stats", "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineBamStats.loadPicardAlignmentStats(infiles, outfile)

#########################################################################


@transform(dedup, regex(r"bam/(\S+).dedup.bam"), r"bam/\1.dedup.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a
       bed file using Picard'''

    # TS check whether this is always required or specific to current baits
    # file

    # baits file requires modification to make picard accept it
    # this is performed before CalculateHsMetrics
    to_cluster = USECLUSTER
    baits = PARAMS["roi_baits"]
    modified_baits = infile + "_temp_baits_final.bed"
    regions = PARAMS["roi_regions"]
    statement = '''samtools view -H %(infile)s > %(infile)s_temp_header.txt;
                awk 'NR>2' %(baits)s |
                awk -F '\\t' 'BEGIN { OFS="\\t" } {print $1,$2,$3,"+",$4;}'
                > %(infile)s_temp_baits.bed;
                cat  %(infile)s_temp_header.txt %(infile)s_temp_baits.bed
                > %(modified_baits)s; checkpoint ;
                rm -rf %(infile)s_temp_baits.bed %(infile)s_temp_header.txt
                '''
    P.run()

    PipelineBamStats.buildPicardCoverageStats(
        infile, outfile, modified_baits, modified_baits)

    #IOTools.zapFile(modified_baits)


@follows(buildCoverageStats)
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    PipelineBamStats.loadPicardCoverageStats(infiles, outfile)

#########################################################################
#########################################################################
#########################################################################
# GATK realign bams
#########################################################################

@follows(loadCoverageStats, mkdir("gatk"))
@transform(dedup,
           regex(r"bam/(\S+).dedup.bam"),
            r"gatk/\1.readgroups.bam")
def GATKReadGroups(infile, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = re.sub(r'-\w+-\w+\.bam', '', os.path.basename(infile))
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = PARAMS["gatk_threads"]
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    genome = PARAMS["bwa_index_dir"] + PARAMS["genome"] 
    
    PipelineExome.GATKReadGroups(infile, outfile, genome,
                                 library, platform,
                                 platform_unit, track)
    #IOTools.zapFile(infile)
    
    
###############################################################################


@transform(GATKReadGroups,
           regex(r"gatk/(\S+).readgroups.bam"),
           r"gatk/\1.bqsr.bam")
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    #intrack = P.snip(os.path.basename(infile), ".bam")
    #outtrack = P.snip(os.path.basename(outfile), ".bam")
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["bwa_index_dir"] + PARAMS["genome"] 
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    
    PipelineExome.GATKBaseRecal(infile, outfile, genome, intervals,
                                    padding, dbsnp, solid_options)
   # IOTools.zapFile(infile)


#########################################################################
#########################################################################
#########################################################################
# Variant Calling
#########################################################################


@follows(mkdir("normal_panel_variants", GATKBaseRecal))
@transform(GATKBaseRecal,
           regex(r"gatk/(\S+)-%s-(\S+).bqsr.bam" %
                 PARAMS["sample_control"]),
           r"normal_panel_variants/\1_normal_mutect.vcf")
def callControlVariants(infile, outfile):
    '''run mutect to call snps in control sample'''

    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    cosmic = PARAMS["mutect_cosmic"]
    dbsnp = PARAMS["mutect_dbsnp"]
    roi_intervals = PARAMS["roi_intervals"]
    threads = PARAMS['mutect_threads']
    memory = PARAMS['mutect_memory']

    genome = PARAMS["bwa_index_dir"] + PARAMS["genome"]

    PipelineExome.MuTect2Caller(infile, outfile, mutect_log, genome, roi_intervals, 
                    cosmic, dbsnp, call_stats_out, memory, threads, artifact=True)


@transform(callControlVariants,
           regex(r"normal_panel_variants/(\S+).vcf"),
           r"normal_panel_variants/\1_slim.vcf.gz")
def indexControlVariants(infile, outfile):
    '''index control vcf for intersection by vcftools'''
    '''tabix is a tool that allows to perform fast interval queries based on tab delimited interval file'''

    outfile = P.snip(outfile, ".gz")

    statement = '''cut -f1-8 %(infile)s > %(outfile)s;
                   bgzip -f %(outfile)s;
                   tabix -f %(outfile)s.gz'''
    P.run()


# paramaterise vcf intersection (number of req. observations - currently 1)
@merge(indexControlVariants,
       "normal_panel_variants/combined.vcf")
def mergeControlVariants(infiles, outfile):
    ''' intersect control vcfs to generate a panel of normals for mutect'''
    '''outputs positions present in at least one file'''
    infiles = " ".join(infiles)

    statement = '''module load vcftools/0.1.14;
                   module load perl;
                   export PERL5LIB=/package/vcftools/0.1.14/lib/site_perl/5.18.1;
                   vcf-isec -n +1 %(infiles)s
                   > %(outfile)s'''
    P.run()

@follows(callControlVariants)
@collate(GATKBaseRecal,
           regex(r"gatk/(CM[0-9]{4})(\S+).bqsr.bam"),
           r"gatk/\1.pt")
def patientID(infiles, outfile):
    '''makes and empty file for patient ID'''
    '''patient sample names should start with CM'''
    '''might need to change it for different patient names'''

    to_cluster = False
    statement = '''touch %(outfile)s'''
    
    P.run()

@follows(mkdir("variants"), patientID)
@transform(patientID,
           regex(r"gatk/(.*).pt"),
           r"variants/\1.mutect.vcf")
def runMutect(infile, outfile):
    '''calls somatic SNPs and indels using MuTect2'''
    
    base = P.snip(os.path.basename(infile), ".pt")
    infile_tumour = base + "-Tumour-R1.bqsr.bam"
    infile_tumour_key = base + "-Tumour-R1"
    infile_control = SAMPLE_DICT[infile_tumour_key] + ".bqsr.bam"

    basename = P.snip(outfile, ".mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, strand_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_lod"], PARAMS["mutect_strand_lod"])
    
    threads = PARAMS['mutect_threads']
    memory = PARAMS['mutect_memory']  

    genome = PARAMS["bwa_index_dir"] + PARAMS["genome"]

    PipelineExome.MuTect2Caller(
        infile_tumour, infile_control, outfile, mutect_log, genome,
        cosmic, dbsnp, call_stats_out,
        memory, threads,
        quality, max_alt_qual,
        max_alt, max_fraction, tumor_LOD, strand_LOD,
        infile_matched=True)


@transform(runMutect,
           regex(r"variants/(\S+).mutect.snp.vcf"),
           r"variants/\1_call_stats.load")
def loadMutectExtendedOutput(infile, outfile):
    '''Load mutect extended output into database'''

    infile = infile.replace(".mutect.snp.vcf", "_call_stats.out")

    indices = "contig,position"
    P.load(infile, outfile, options="--add-index=%(indices)s" % locals())


@transform(GATKBaseRecal,
           regex(r"bam/(\S+)-%s-(\S+).realigned.split.bqsr.bam" %
                 PARAMS["sample_control"]),
           r"variants/\1/results/all.somatic.indels.vcf")
def indelCaller(infile, outfile):
    '''Call somatic indels using Strelka'''

    infile_tumour = infile.replace(
        PARAMS["sample_control"], PARAMS["sample_tumour"])
    outdir = "/".join(outfile.split("/")[0:2])
    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineExome.strelkaINDELCaller(infile, infile_tumour, outfile,
                                     genome, PARAMS['strelka_config'], outdir,
                                     PARAMS['strelka_memory'],
                                     PARAMS['strelka_threads'])

##########################################################################
##########################################################################
##########################################################################
# repeat mutect in reverse and on subsampled control bam as quality control
##########################################################################
# this analysis should be part of an optional check of mutect parameters
# mutect paramters should be identical to the runMutect function above
# splitMergedRealigned replaced by GATKBaseRecal

@follows(mergeControlVariants)
@transform(GATKBaseRecal,
           regex(r"bam/(\S+)-%s-(\S+).realigned.split.bqsr.bam" %
                 PARAMS["sample_control"]),
           add_inputs(mergeControlVariants),
           r"variants/\1.mutect.reverse.snp.vcf")
def runMutectReverse(infiles, outfile):
    '''Use control as tumor and vis versa to estimate false positive rate'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        PARAMS["sample_control"], PARAMS["sample_tumour"])

    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    basename = P.snip(outfile, ".mutect.reverse.snp.vcf")
    call_stats_out = basename + "_call_stats.reverse.out"
    coverage_wig_out = basename + "_coverage.reverse.wig"
    mutect_log = basename + ".reverse.log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineExome.mutectSNPCaller(infile, outfile, mutect_log, genome,
                                  cosmic, dbsnp, call_stats_out,
                                  PARAMS['mutect_memory'],
                                  PARAMS['mutect_threads'],
                                  quality, max_alt_qual,
                                  max_alt, max_fraction, tumor_LOD,
                                  normal_panel, infile_tumour)


# generalise the functions below
# 1. identify sample with highest coverage in control
# - should this check coverage in tumour also?
# 2. subset control bam
# 3. run mutect calling function with subset against unsubsetted tumour
# 4. summary table

adeno_bam = "bam/NU16C-Control-1.realigned.bqsr.bam"


@subdivide(adeno_bam,
           regex("(\S+).bqsr.bam"),
           [r"\1.0.1.bqsr.bam",
            r"\1.0.2.bqsr.bam",
            r"\1.0.3.bqsr.bam",
            r"\1.0.4.bqsr.bam",
            r"\1.0.5.bqsr.bam",
            r"\1.0.6.bqsr.bam",
            r"\1.0.7.bqsr.bam",
            r"\1.0.8.bqsr.bam",
            r"\1.0.9.bqsr.bam",
            r"\1.1.0.bqsr.bam"])
def subsetControlBam(infile, outfiles):
    statements = []
    n = 0
    for fraction in numpy.arange(0.1, 1.1, 0.1):
        outfile = outfiles[n]
        n += 1
        statement = '''samtools view -s %(fraction)s -b %(infile)s
                     > %(outfile)s'''
        P.run()


@transform(subsetControlBam,
           suffix(".bam"),
           ".bam.bai")
def indexSubsets(infile, outfile):
    statement = '''samtools index %(infile)s'''
    P.run()


@follows(indexSubsets)
@transform(subsetControlBam,
           regex(r"bam/(\S+)-%s-1.realigned.(\S+).bqsr.bam" %
                 PARAMS["sample_control"]),
           add_inputs(mergeControlVariants),
           r"variants/\1-downsampled-\2.mutect.snp.vcf")
def runMutectOnDownsampled(infiles, outfile):
    '''call somatic SNPs using MuTect on downsampled bams'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        PARAMS["sample_control"], PARAMS["sample_tumour"])
    basename = P.snip(outfile, "_normal_mutect.vcf")

    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineExome.mutectSNPCaller(infile_tumour, outfile, mutect_log, genome,
                                  cosmic, dbsnp, call_stats_out,
                                  PARAMS['mutect_memory'], PARAMS[
                                      'mutect_threads'],
                                  quality, max_alt_qual,
                                  max_alt, max_fraction, tumor_LOD,
                                  normal_panel, infile)

##############################################################################
##############################################################################
##############################################################################
# Variant Annotation and Recalibration
##############################################################################


@collate(GATKBaseRecal,
         regex(r"bam/(\S+)-(\S+)-(\S+).realigned.split.bqsr.bam"),
         r"bam/\1.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for each patient,
       for use in variant calling'''
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            infile_tumour = infile.replace(
                PARAMS["sample_control"], PARAMS["sample_tumour"])
            outf.write(infile + '\n')
            outf.write(infile_tumour + '\n')


@transform(runMutect,
           regex(r"variants/(\S+).mutect.snp.vcf"),
           r"variants/\1.mutect.snp.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate SNP variants using SNPeff'''
    to_cluster = USECLUSTER
    job_memory = "4G"
    job_threads = 2

    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''java -Xmx4G -jar /ifs/apps/bio/snpEff-3.3-dev/snpEff.jar
                   -c %(config)s -v %(snpeff_genome)s -o gatk
                   %(infile)s > %(outfile)s'''
    P.run()


@transform(indelCaller,
           regex("variants/(\S+)/results/all.somatic.indels.vcf"),
           r"variants/\1.indels.snpeff.vcf")
def annotateVariantsINDELsSNPeff(infile, outfile):
    '''Annotate INDEL variants using SNPeff'''
    to_cluster = USECLUSTER
    job_memory = "4G"
    job_threads = 2

    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''java -Xmx4G -jar /ifs/apps/bio/snpEff-3.3-dev/snpEff.jar
                   -c %(config)s -v %(snpeff_genome)s -o gatk
                   %(infile)s > %(outfile)s'''
    P.run()


#########################################################################
# Annotate SNP and INDEL variants
#########################################################################

# Need to check whether variant annotatot is using both bams
# from a single patient?
# should just be the tumour bam or else scores will be wrong!

@follows(annotateVariantsSNPeff, listOfBAMs)
@transform(runMutect,
           regex(r"variants/(\S+).mutect.snp.vcf"),
           add_inputs(r"bam/\1.list",
                      r"variants/\1.mutect.snp.snpeff.vcf"),
           r"variants/\1.mutect.snp.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK
                   -T VariantAnnotator
                   -R %(bwa_index_dir)s/%(genome)s.fa
                   -I %(bamlist)s
                   -A SnpEff --snpEffFile %(effFile)s
                   -o %(outfile)s
                   --variant %(infile)s
                   -L %(infile)s
                   --dbsnp %(dbsnp)s
                   -A HaplotypeScore
                   -A MappingQualityRankSumTest
                   -A ReadPosRankSumTest
                   -A AlleleBalanceBySample'''
    P.run()


@follows(annotateVariantsINDELsSNPeff, listOfBAMs)
@transform(indelCaller,
           regex("variants/(\S+)/results/all.somatic.indels.vcf"),
           add_inputs(r"bam/\1.list", r"variants/\1.indels.snpeff.vcf"),
           r"variants/\1.indels.annotated.vcf")
def variantAnnotatorIndels(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    statement = '''GenomeAnalysisTK
                   -T VariantAnnotator
                   -R %(bwa_index_dir)s/%(genome)s.fa
                   -I %(bamlist)s
                   -A SnpEff --snpEffFile %(effFile)s
                   -o %(outfile)s
                   --variant %(infile)s
                   -L %(infile)s
                   -A Coverage
                   -A FisherStrand
                   -A HaplotypeScore
                   -A MappingQualityRankSumTest
                   -A ReadPosRankSumTest
                   -A AlleleBalanceBySample
                   -A RMSMappingQuality'''
    P.run()


######################################################################

# this does not work - insufficient number of indels in mills+
# therefore this task is not a dependency of task full
@transform(variantAnnotatorIndels,
           suffix(".annotated.vcf"),
           ".annotated.recalibrated.vcf")
def variantRecalibrator(infile, outfile):
    '''Create variant recalibration file for indels'''
    to_cluster = USECLUSTER
    job_memory = PARAMS["gatk_memory"]
    job_threads = 6
    track = P.snip(os.path.basename(outfile), ".annotated.recalibrated.vcf")
    mills = PARAMS["gatk_mills"]

    statement = '''GenomeAnalysisTK
                   -T VariantRecalibrator
                   -R %(bwa_index_dir)s/%(genome)s.fa
                   -input %(infile)s
                   -resource:mills,known=true,training=true,truth=true,prior=12.0
                   %(mills)s
                   -an DP -an MQRankSum -an ReadPosRankSum
                   -mode INDEL
                   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
                   --maxGaussians 4
                   -recalFile %(outfile)s
                   -tranchesFile variants/%(track)s.tranches
                   -rscriptFile variants/%(track)s.plots.R'''
    P.run()

##############################################################################
# Filter SNPs and INDELs
##############################################################################


@transform(variantAnnotatorIndels,
           suffix(".annotated.vcf"),
           ".annotated.filtered.vcf")
def filterIndels(infile, outfile):
    ''' use SnpSift to filter INDELS using VCF fields'''

    statement = '''cat %(infile)s |
                   java -Xmx2g -jar /ifs/apps/bio/snpEff-3.1/SnpSift.jar filter
                   "(QSI_NT>%(filter_indel_nt)s &
                     IHP<%(filter_indel_ihp)s &
                     RC<%(filter_indel_rc)s &
                     IC<%(filter_indel_rc)s) "
                   > %(outfile)s '''
    P.run()


@transform(variantAnnotator,
           regex("variants/(\S+).mutect.snp.annotated.vcf"),
           r"variants/\1.mutect.snp.annotated.filtered.vcf")
def filterMutect(infile, outfile):
    ''' filter mutect snps using allele frequencies '''

    logfile = outfile.replace(".vcf", ".log")

    min_t_alt = PARAMS["filter_minimum_tumor_allele"]
    min_t_alt_freq = PARAMS["filter_minimum_tumor_allele_frequency"]
    min_n_depth = PARAMS["filter_minimum_normal_depth"]
    max_n_alt_freq = PARAMS["filter_maximum_normal_allele_frequency"]
    min_ratio = PARAMS["filter_minimum_ratio"]

    PipelineExome.filterMutect(
        infile, outfile, logfile,
        PARAMS["sample_control"], PARAMS["sample_tumour"],
        min_t_alt, min_n_depth, max_n_alt_freq,
        min_t_alt_freq, min_ratio)

##############################################################################
# Intersect filtered SNPs and INDELs
##############################################################################


@mkdir("intersection.dir")
@collate((filterIndels, filterMutect),
         regex(r"variants/(\S+)\.(\S+).annotated.filtered.vcf"),
         r"intersection.dir/overlap_\2_heatmap.png")
def intersectHeatmap(infiles, outfile):
    ''' intersect DE test_ids across the different quantifiers'''

    PipelineExome.intersectionHeatmap(infiles, outfile)

#########################################################################
#########################################################################
# convert vcf to tsv files and load into database


@transform(filterMutect,
           regex("variants/(\S+).annotated.filtered.vcf"),
           r"variants/\1.annotated.filtered.tsv")
def snpvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    statement = '''GenomeAnalysisTK
                   -T VariantsToTable -R %(bwa_index_dir)s/%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -GF GT -GF AD -GF SS -GF FA -GF AB -GF DP
                   -o %(outfile)s'''
    P.run()


@transform(filterIndels,
           regex("variants/(\S+).annotated.filtered.vcf"),
           r"variants/\1.annotated.filtered.tsv")
def indelvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    statement = '''GenomeAnalysisTK
                   -T VariantsToTable -R %(bwa_index_dir)s/%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -F TQSI -F TSQI_NT -F DP -F IC -F IHP -F NT
                   -F QSI -F QSI_NT -F RC -F RU -F SGT
                   -GF DP -GF DP2 -GF DP50 -GF SUBDP50 -GF TAR -GF TIR -GF TOR
                   -o %(outfile)s'''
    P.run()


@transform([snpvcfToTable,
            indelvcfToTable],
           regex(r"variants/(\S+).annotated.filtered.tsv"),
           r"variants/\1_annotated.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''

    if infile.endswith("indels.annotated.filtered.tsv"):
        indices = "CHROM,POS,SNPEFF_GENE_NAME"
    elif infile.endswith("mutect.snp.annotated.filtered.tsv"):
        indices = "CHROM,POS,SNPEFF_GENE_NAME"

    P.load(infile, outfile, options="--add-index=%(indices)s" % locals())

#########################################################################
# Genes of interest
# check this will run in the correct position if option selected

# @active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
# @transform((annotateVariantsSNPsift),
#           regex(r"variants/(\S+).haplotypeCaller.snpsift.vcf"),
#           r"variants/\1.genes.vcf")
# def findGenes(infile, outfile):
#    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
#    if variant is within a gene of interest as defined in the ini
#    file'''
#
#    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
#    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
#    statement = '''GenomeAnalysisTK -T VariantFiltration
#    -R %%(bwa_index_dir)s/%%(genome)s.fa
#    --variant %(infile)s
#    --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
#    --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
#    P.run()

#########################################################################
#########################################################################
#########################################################################
# vcf statistics -   this only summarises the nucleotide changes
# this currently does not provide useful output!


@transform((variantAnnotator,
            variantAnnotatorIndels),
           regex(r"variants/(\S+).vcf"),
           r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    to_cluster = USECLUSTER
    statement = '''vcf-stats %(infile)s
                   > %(outfile)s 2>>%(outfile)s.log;'''
    P.run()


@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    csv2db_options = PARAMS["csv2db_options"]
    E.info("Loading vcf stats...")
    statement = '''cgat vcfstats2db
                   %(filenames)s >> %(outfile)s; '''
    statement += '''cat vcfstats.txt |
                    cgat csv2db %(csv2db_options)s
                    --allow-empty-file --add-index=track --table=vcf_stats
                    >> %(outfile)s; '''
    P.run()

#########################################################################


@transform(runMutect,
           suffix(".mutect.snp.vcf"),
           "_mutect_filtering_summary.tsv")
def summariseFiltering(infile, outfile):
    infile = infile.replace(".mutect.snp.vcf", "_call_stats.out")

    PipelineExome.parseMutectCallStats(infile, outfile, submit=True)


@transform(summariseFiltering,
           regex(r"variants/(\S+)_mutect_filtering_summary.tsv"),
           r"variants/\1_mutect_filtering_summary.load")
def loadMutectFilteringSummary(infile, outfile):
    '''Load mutect extended output into database'''

    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   cgat csv2db
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@originate("eBio_studies.tsv")
def defineEBioStudies(outfile):
    ''' For the cancer types specified in pipeline.ini, identify the
    relevent studies in eBio '''

    cancer_types = PARAMS["annotation_ebio_cancer_types"]

    PipelineExome.defineEBioStudies(cancer_types, outfile, submit=False)


@transform(defineEBioStudies,
           suffix("eBio_studies.tsv"),
           add_inputs(filterIndels, filterMutect),
           "eBio_studies_gene_frequencies.tsv")
def extractEBioinfo(infiles, outfile):
    '''find the number of mutations identitified in previous studies (ebio_ids)
    for the mutated genes in the annotated vcfs'''

    eBio_ids = infiles[0]
    vcfs = infiles[1:]

    PipelineExome.extractEBioinfo(eBio_ids, vcfs, outfile, submit=False)


@transform(extractEBioinfo,
           suffix(".tsv"),
           ".load")
def loadEBioInfo(infile, outfile):
    '''load the frequencies from the eBIO portal'''

    P.load(infile, outfile, options="--add-index=gene")

#########################################################################
#########################################################################
#########################################################################
# load Network of Cancer Genes table

# parameterise file location:


@originate("cancergenes.load")
def loadNCG(outfile):
    '''Load NCG into database'''

    infile = PARAMS["cancergenes_table"]
    # infile = "/ifs/projects/proj053/backup/NCG/cancergenes2016.tsv"

    P.load(infile, outfile, options="--add-index=symbol")

#########################################################################
#########################################################################
#########################################################################
# analyse mutational siganture of filtered variants


@merge(filterMutect,
       ["variants/mutational_signature.tsv",
        "variants/mutational_signature_table.tsv"])
def mutationalSignature(infiles, outfiles):

    PipelineExome.compileMutationalSignature(
        infiles, outfiles)


@transform(mutationalSignature,
           suffix(".tsv"),
           ".load")
def loadMutationalSignature(infiles, outfile):
    outfile2 = re.sub(".load", "_table.load", outfile)
    P.load(infiles[0], outfile)
    P.load(infiles[1], outfile2)


#########################################################################
#########################################################################
#########################################################################


@follows(defineEBioStudies)
def test():
    pass


@follows(runMutectOnDownsampled,
         runMutectReverse)
def TestMutect():
    '''This target runs function which can be used to assess the chosen
    mutect parameters'''


# @follows(loadROI,
#         loadROI2Gene)
# def loadMetadata():
#    pass


@follows(mapReads)
def mapping():
    pass


@follows(dedup,
         loadPicardDuplicateStats,
         buildPicardAlignStats,
         loadPicardAlignStats,
         buildCoverageStats,
         loadCoverageStats)
def postMappingQC():
    pass


@follows(GATKReadGroups,
         GATKBaseRecal)
def gatk():
    pass


@follows(runMutect,
         indelCaller)
def callVariants():
    pass


@follows(loadVariantAnnotation)
def tabulation():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(postMappingQC,
         loadManualAnnotations,
         loadMutectFilteringSummary,
         loadMutectExtendedOutput,
         loadVariantAnnotation,
         loadCoverageStats,
         loadPicardAlignStats,
         loadNCG,
         loadMutationalSignature,
         loadEBioInfo,
         intersectHeatmap)
def full():
    pass

#########################################################################
#########################################################################
#########################################################################


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''
    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''
    E.info("updating documentation")
    P.run_report(clean=False)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
