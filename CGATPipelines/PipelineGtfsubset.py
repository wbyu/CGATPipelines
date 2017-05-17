"""
PipelineGtfsubset.py - Tasks for GTF subsetting
===============================================

Reference
---------

"""

import CGAT.Experiment as E
import os
import MySQLdb
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGATPipelines.Pipeline as P
import CGAT.GFF3 as GFF3


class SubsetGTF():

    def __init__(self, infile, *args, **kwargs):

        self.gtf = GTF.iterator(IOTools.openFile(infile, "r"))

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

    def filterGTF(self, outfile, filteroption, filteritem, operators):
        '''

        '''

        with IOTools.openFile(outfile, "w") as outf:
            for line in self.gtf:
                D = self.makeLineDict(line)
                if len(filteritem) == 1:
                    if D[filteroption] == filteritem[0]:
                        outf.write("%s\n" % str(line))

                elif len(filteritem) == 2:
                    filter1, filter2 = filteritem
                    filteroption1, filteroption2 = filteroption
                    if operators == "and":
                        if D[filteroption1] == filter1 and \
                           D[filteroption2] == filter2:
                            outf.write("%s\n" % str(line))
                    elif operators == "and not":
                        if D[filteroption1] == filter1 and not\
                           D[filteroption2] == filter2:
                            outf.write("%s\n" % str(line))

                else:
                    pass

        outf.close()


class SubsetGFF3():

    def __init__(self, infile, *args, **kwargs):
        self.gff = GFF3.flat_file_iterator(IOTools.openFile(infile, "r"))

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

    def filterGFF3(self, outfile, filteroption, filteritem):

        with IOTools.openFile(outfile, "w") as outf:
            for line in self.gff:
                D = self.makeLineDict(line)
                if len(filteritem) == 1:
                    if D[filteroption] == filteritem[0]:
                        outf.write("%s\n" % str(line))

        outf.close()


def connectToUCSC(host="genome-mysql.cse.ucsc.edu",
                  user="genome",
                  database=None):
    """connect to UCSC database.

    Arguments
    ---------
    host : string
        Host to connect to
    user : string
        Username to connect with
    Database : string
        database to use

    Returns
    -------
    Database handle

    """
    dbhandle = MySQLdb.Connect(host=host,
                               user=user)

    cc = dbhandle.cursor()
    cc.execute("USE %s " % database)

    return dbhandle


def getRepeatDataFromUCSC(dbhandle,
                          repclasses,
                          outfile,
                          remove_contigs_regex=None):
    '''download data from UCSC database and write to `outfile` in
    :term:`gff` format.

    This method downloads repeats from the repeatmasker track at
    the UCSC.

    Arguments
    ---------
    dbhandle : object
       Database handle to UCSC mysql database
    repclasses : list
       List of repeat classes to select. If empty, all repeat classes
       will be collected.
    outfile : string
       Filename of output file in :term:`gff` format.
    remove_contigs_regex : string
       If given, remove repeats on contigs matching the regular
       expression given.

    '''
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    # now collect repeats
    tmpfile = P.getTempFile(".")

    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd,
        '.', strand, '.',
        CONCAT('class \\"', repClass, '\\"; family \\"',
        repFamily, '\\"; repName \\"', repName, '\\";')
        FROM %(table)s"""

        if repclasses:
            repclasses_str = ",".join(
                ["'" + x.strip() + "'" for x in repclasses])
            sql += ''' WHERE repClass in (%(repclasses_str)s) ''' % locals()

        sql = sql % locals()

        E.debug("executing sql statement: %s" % sql)
        cc.execute(sql)
        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")

    tmpfile.close()

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = ['''cat %(tmpfilename)s
    | %(pipeline_scriptsdir)s/gff_sort pos
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log ''']

    if remove_contigs_regex:
        statement.append(
            ''' --contig-pattern="%(remove_contigs_regex)s" ''')

    statement.append('''| gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()

    os.unlink(tmpfilename)
