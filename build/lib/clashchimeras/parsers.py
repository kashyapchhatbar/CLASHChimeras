import csv
import gzip
import logging
import mmap
import os
import sys
import textwrap
from collections import Counter, defaultdict
from itertools import groupby
from operator import itemgetter

import pandas as pd
import pyfaidx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import clashchimeras

logger = logging.getLogger('root')


class GFF:
    """GFF file parser for mirbase gff3 file

    This class uses memory-mapped file object to read a mirbase gff3 file. It
    contains methods to read, process a gff3 file and return genomic coordinates

    Attributes:
      fileName: A mirbase gff3 file path
    """

    def __init__(self, fileName=None):
        self.features = {}
        self.fileName = fileName

    def read(self, featureType='miRNA_primary_transcript'):
        """Reads gff3 file provided during class initialization

        Stores the byte positions of every feature in a dict object named
        self.features

        Keyword Args:
          featureType: Feature type of a gff3 record, the third element of every
            record in the file. Please change this if you want to store mature
            form of microRNA, by default it uses primary transcript
            (default 'miRNA_primary_transcript')
        """
        logger.info('Reading %s' % self.fileName)
        self.fileHandle = open(self.fileName, 'r+b')
        bytePosition = self.fileHandle.tell()
        for line in self.fileHandle:
            row = line.decode('utf-8').rstrip().split("\t")
            if not row[0].startswith("#") and row[2] == featureType:
                attributes = row[-1].split(";")
                for attribute in attributes:
                    if attribute.startswith('Name'):
                        mirbase_name = attribute.split("=")[-1]
                        self.features[mirbase_name] = bytePosition
            bytePosition = self.fileHandle.tell()
        self.fileHandle.close()
        logger.debug('Reading %s finished' % self.fileName)

    def process(self, name):
        """A method to return a Record object providing genomic information

        Args:
          name: A valid miRNA_primary_transcript name

        Returns:
          An object Record containing scaffold, start, end, strand, mirbase_id and
            mirbase_name as its variables for access
        """
        self.fileHandle = open(self.fileName, 'r+b')
        self.mm = mmap.mmap(self.fileHandle.fileno(), 0)
        self.mm.seek(self.features[name])
        row = self.mm.readline().decode('utf-8').rstrip().split("\t")
        attributes = row[-1].split(";")
        for attribute in attributes:
            if attribute.startswith("ID"):
                _id = attribute.split("=")[-1]
            elif attribute.startswith("Name"):
                _name = attribute.split("=")[-1]
        record = Record(scaffold=row[0], start=int(row[3]), end=int(row[4]),
                        strand=row[6], mirbase_id=_id, mirbase_name=_name)
        self.fileHandle.close()
        return record

    def coordinates(self, name, start=None, end=None):
        """A method to return a bed record containing genomic coordinates for the
          aligned segment

        Keyword Args:
          start: The alignment start position of the cDNA molecule or the relative
            start of the particular molecule
          end: The alignment end position in the cDNA molecule or the relative end
            of the particular molecule

        Args:
          name: A valid miRNA_primary_transcript name

        Returns:
          A tuple of strings containing elements for a bed record
        """
        record = self.process(name)
        if not start and not end:
            start = 1
            end = record.end - record.start + 1
        positions = {}
        match_positions = []

        if record.strand == '+':
            _start = 1
            for relative, actual in enumerate(range(record.start - 1, record.end),
                                              start=_start):
                positions[relative] = actual
            for pos in range(start, end + 1):
                match_positions.append(positions[pos])
            return [(record.scaffold, min(match_positions), max(match_positions) + 1,
                     record.mirbase_name, 0, record.strand)]

        elif record.strand == '-':
            _start = 1
            for relative, actual in enumerate(reversed(range(record.start - 1,
                                                             record.end)), start=_start):
                positions[relative] = actual
            for pos in range(start, end + 1):
                match_positions.append(positions[pos])
            return [(record.scaffold, min(match_positions), max(match_positions) + 1,
                     record.mirbase_name, 0, record.strand)]


class GTF:
    """GTF file parser for gencode gtf file

    This class uses memory-mapped file object to read a gencode gtf file. It
    contains methods to read, process a gtf file and return genomic coordinates

    Attributes:
      fileName: A gencode gtf file path
    """

    def __init__(self, fileName=None):
        self.features = defaultdict(list)
        self.biotypeFeatures = defaultdict(list)
        self.geneFeatures = defaultdict(list)
        self.fileName = fileName
        self.geneIds = {}

    def readBiotype(self, featureType='exon', biotype=None):
        logger.info('Reading %s' % self.fileName)
        self.fileHandle = open(self.fileName, 'r+b')
        for line in self.fileHandle:
            row = line.decode('utf-8').rstrip().split("\t")

            if not row[0].startswith("#") and row[2] == featureType:
                attributes = row[-1].split("; ")

                havana_transcript = '-'
                havana_gene = '-'
                exon_number = '0'

                for attribute in attributes:

                    if attribute.startswith("transcript_id"):
                        transcript_id = attribute.split(" ")[-1][1:-1]

                    elif attribute.startswith("transcript_type"):
                        transcript_type = attribute.split(" ")[-1][1:-1]

                    elif attribute.startswith("exon_number"):
                        exon_number = int(attribute.split(" ")[-1])

                    elif attribute.startswith("havana_gene"):
                        havana_gene = attribute.split(" ")[-1][1:-1]

                    elif attribute.startswith("havana_transcript"):
                        havana_transcript = attribute.split(" ")[-1][1:-2]

                    elif attribute.startswith("gene_id"):
                        gene_id = attribute.split(" ")[-1][1:-1]

                    elif attribute.startswith("gene_name"):
                        gene_name = attribute.split(" ")[-1][1:-1]

                    elif attribute.startswith("transcript_name"):
                        transcript_name = attribute.split(" ")[-1][1:-1]

                if biotype == 'tRNA':
                    if transcript_type == "tRNAscan":
                        self.biotypeFeatures[transcript_id].append((exon_number, row[0],
                                                                    int(row[3]), int(row[4]),
                                                                    row[6], gene_id,
                                                                    havana_gene,
                                                                    havana_transcript,
                                                                    transcript_name,
                                                                    gene_name))
                else:
                    if transcript_type == biotype:
                        self.biotypeFeatures[transcript_id].append((exon_number, row[0],
                                                                    int(row[3]), int(row[4]),
                                                                    row[6], gene_id,
                                                                    havana_gene,
                                                                    havana_transcript,
                                                                    transcript_name,
                                                                    gene_name))

        self.fileHandle.close()

    def read(self, featureType='exon'):
        """Reads gtf file provided during class initialization

        Stores the byte positions of every feature in a defaultdict(list) object
        named self.features

        Keyword Args:
          featureType: Feature type of a gtf record, the third element of every
            record in the file. Please change this if you want to get specific
            records (e.g. 'UTR') (default 'exon')
        """

        logger.info('Reading %s' % self.fileName)
        self.fileHandle = open(self.fileName, 'r+b')
        bytePosition = self.fileHandle.tell()
        for line in self.fileHandle:
            row = line.decode('utf-8').rstrip().split("\t")
            if not row[0].startswith("#") and row[2] == featureType:
                attributes = row[-1].split("; ")
                for attribute in attributes:
                    if attribute.startswith("transcript_id"):
                        transcript_id = attribute.split(" ")[-1][1:-1]
                        self.features[transcript_id].append(bytePosition)
                        self.geneIds[transcript_id] = gene_id
                    if attribute.startswith("gene_id"):
                        gene_id = attribute.split(" ")[-1][1:-1]
            bytePosition = self.fileHandle.tell()
        self.fileHandle.close()

        self.fileHandle = open(self.fileName, 'r+b')
        bytePosition = self.fileHandle.tell()
        for line in self.fileHandle:
            row = line.decode('utf-8').rstrip().split("\t")
            if not row[0].startswith("#") and row[2] == featureType:
                attributes = row[-1].split("; ")
                for attribute in attributes:
                    if attribute.startswith("gene_id"):
                        gene_id = attribute.split(" ")[-1][1:-1]
                        self.geneFeatures[gene_id].append(bytePosition)
            bytePosition = self.fileHandle.tell()
        self.fileHandle.close()

        logger.debug('Reading %s finished' % self.fileName)

    def process(self, name):
        """A method to return a Record object providing genomic information

        Args:
          name: A valid gencode transcript_id

        Returns:
          An object Record containing scaffold, start, end, strand, mirbase_id and
            mirbase_name as its variables for access
        """
        self.fileHandle = open(self.fileName, 'r+b')
        self.mm = mmap.mmap(self.fileHandle.fileno(), 0)
        positions = self.features[name]
        for position in positions:
            self.mm.seek(position)
            row = self.mm.readline().decode('utf-8').rstrip().split("\t")
            attributes = row[-1].split("; ")
            _eid = '-'
            _enb = '0'
            for attribute in attributes:
                if attribute.startswith("transcript_type"):
                    _tt = attribute.split(" ")[-1][1:-1]
                elif attribute.startswith("transcript_id"):
                    _tid = attribute.split(" ")[-1][1:-1]
                elif attribute.startswith("exon_id"):
                    _eid = attribute.split(" ")[-1][1:-1]
                elif attribute.startswith("exon_number"):
                    _enb = int(attribute.split(" ")[-1])
                elif attribute.startswith("gene_name"):
                    _gn = attribute.split(" ")[-1][1:-1]
            record = Record(scaffold=row[0], start=int(row[3]), end=int(row[4]),
                            strand=row[6], transcript_type=_tt, transcript_id=_tid, exon_id=_eid,
                            exon_number=_enb, gene_name=_gn)
            yield record
        self.fileHandle.close()

    def geneExonicRegions(self, df):
        """Given a DataFrame with the exon coordinates from Gencode for a single
        gene, return the total number of coding regions in that gene.
        """
        scaffold = df.iloc[0].scaffold
        strand = df.iloc[0].strand
        gene_type = df.iloc[0].gene_type
        gene_id = df.iloc[0].gene_id
        gene_name = df.iloc[0].gene_name
        start = df.start.min()
        end = df.end.max()
        bp = [False] * (end - start + 1)
        for i in range(df.shape[0]):
            s = df.iloc[i]['start'] - start
            e = df.iloc[i]['end'] - start + 1
            bp[s:e] = [True] * (e - s)
        regions = list(range(start, end + 1))
        groups = []

        for i, j in groupby(bp):
            groups.append((i, len(list(j))))
        e_start = 0

        for i in groups:
            e_end = e_start + i[1]
            if i[0]:
                record = Record(scaffold=scaffold, start=regions[e_start],
                                end=regions[e_end - 1], gene_type=gene_type, gene_id=gene_id,
                                gene_name=gene_name, strand=strand)
                yield record
            e_start += i[1]

    def geneProcess(self, name):
        """A method to return a Record object providing genomic information

        Args:
          name: A valid gencode gene_id

        Returns:
          An object Record containing scaffold, start, end, strand, mirbase_id and
            mirbase_name as its variables for access
        """
        self.fileHandle = open(self.fileName, 'r+b')
        self.mm = mmap.mmap(self.fileHandle.fileno(), 0)
        positions = self.geneFeatures[name]
        exons = []
        for position in positions:
            self.mm.seek(position)
            row = self.mm.readline().decode('utf-8').rstrip().split("\t")
            attributes = row[-1].split("; ")
            for attribute in attributes:
                if attribute.startswith("gene_type"):
                    _gt = attribute.split(" ")[-1][1:-1]
                elif attribute.startswith("gene_id"):
                    _gid = attribute.split(" ")[-1][1:-1]
                elif attribute.startswith("gene_name"):
                    _gn = attribute.split(" ")[-1][1:-1]
            exons.append((row[0], int(row[3]), int(row[4]), row[6], _gt, _gid, _gn))
        self.fileHandle.close()
        exons_df = pd.DataFrame(exons, columns=['scaffold', 'start', 'end',
                                                'strand', 'gene_type', 'gene_id', 'gene_name'])

        for record in self.geneExonicRegions(exons_df):
            yield record

    def coordinates(self, name, start=None, end=None):
        """A generator to return a bed record containing genomic coordinates for the
          aligned segment

        Keyword Args:
          start: The alignment start position of the cDNA molecule or the relative
            start of the particular molecule
          end: The alignment end position in the cDNA molecule or the relative end
            of the particular molecule

        Args:
          name: A valid miRNA_primary_transcript name

        Returns:
          A list of tuple(s) of strings containing elements for a bed record. There
          may be more than one because of alternate splicing.
        """
        if "|" in name:
            self.name = name.split("|")[0]
        else:
            self.name = name
        positions = {}
        match_positions = []
        records = []
        segments = []
        result_segments = []
        for record in self.process(self.name):
            records.append(record)
        records.sort(key=lambda x: int(x.exon_number))

        if records[0].strand == '+':
            _start = 1
            for record in records:
                for relative, actual in enumerate(range(record.start, record.end + 1),
                                                  start=_start):
                    positions[relative] = actual
                _start = relative + 1
            for pos in range(start, end):
                match_positions.append(positions[pos])
            for key, group in groupby(enumerate(match_positions),
                                      lambda x: x[0] - x[-1]):
                segment = list(map(itemgetter(1), group))
                segments.append([segment[0], segment[-1]])
            for segment in segments:
                for record in records:
                    if segment[0] >= record.start and segment[1] <= record.end:
                        result_segments.append((record.scaffold, segment[0], segment[1],
                                                record.transcript_id + '|' + record.gene_name, 0, record.strand))

        elif records[0].strand == '-':
            _start = 1
            for record in records:
                for relative, actual in enumerate(reversed(range(record.start,
                                                                 record.end + 1)), start=_start):
                    positions[relative] = actual
                _start = relative + 1
            for pos in range(start, end):
                match_positions.append(positions[pos])
            for key, group in groupby(enumerate(reversed(match_positions)),
                                      lambda x: x[0] - x[-1]):
                segment = list(map(itemgetter(1), group))
                segments.append([segment[0], segment[-1]])
            for segment in segments:
                for record in records:
                    if segment[0] >= record.start and segment[1] <= record.end:
                        result_segments.append((record.scaffold, segment[0], segment[1],
                                                record.transcript_id + '|' + record.gene_name, 0, record.strand))

        if len(result_segments) == 0:
            logger.debug('%s, %s, %s' % (name, start, end))
            logger.debug('%s' % str(segments))
            for r in records:
                logger.debug('%s %s %s %s' % (r.scaffold, r.strand,
                                              r.start, r.end))

        return result_segments


class SAM:
    """SAM file parser for parsing bowtie2 generated files

    This class uses memory-mapped file object to read a sam file

    Attributes:
      fileName: A sam file path
    """

    def __init__(self, fileName=None):
        self.fileName = fileName
        self.records = {}

    def read(self, flag=0):
        """Reads sam file provided during class initialization

        Stores the byte position of every record based on the keyword arg flag
        provided, to a dict object named self.records

        Keyword Args:
          flag: The SAM alignment flag for a record. For default, it uses the
          primary alignment for every record and ignores secondary alignments
          (default 0)
        """
        logger.info('Reading %s' % self.fileName)
        self.fileHandle = open(self.fileName, 'r+b')
        bytePosition = self.fileHandle.tell()
        for line in self.fileHandle:
            read = line.decode('utf-8').split("\t")
            if not read[0].startswith("@") and read[1] == str(flag):
                self.records[read[0]] = bytePosition
            bytePosition = self.fileHandle.tell()
        self.fileHandle.close()
        logger.debug('Reading %s finished' % self.fileName)

    def access(self, queryName):
        """Provides random access of a record from the sam file

        Args:
          queryName: The query name of the read from the sam file

        Returns:
          A list generated after splitting the record line from sam file
        """
        self.fileHandle = open(self.fileName, 'r+b')
        self.mm = mmap.mmap(self.fileHandle.fileno(), 0)
        self.mm.seek(self.records[queryName])
        row = self.mm.readline().decode('utf-8').rstrip().split("\t")
        self.fileHandle.close()

        return self.pretty(row)

    def filterPotentialChimeras(self, min_length=30, flag=0, target=None):
        """Generated a filtered fasta file from a sam file

        This filtered fasta file contains reads that can be potentially chimeras.
        The criteria for filtering is based on the minimum length

        Keyword Args:
          min_length: To be selected as a potential chimera, this is the minimum
            read length (default 30)
          flag: The SAM alignment flag describing the type of alignment (default 0)
          target: The prefix for output file
        """
        logger.debug('Filtering {} for potential chimeras'.format(target))
        target = '{}.filter.fasta'.format(target.rpartition(".")[0])
        if os.path.exists(target):
            logger.info('Skipping filtering for {}'.format(target))
        else:
            with open(target, 'w') as oH:
                with open(self.fileName) as iH:
                    for row in csv.reader(iH, delimiter="\t"):
                        if not row[0].startswith('@') and row[1] == str(flag):
                            if len(row[9]) >= 30:
                                print(textwrap.fill('>%s' % row[0], width=80), file=oH)
                                print(textwrap.fill('%s' % row[9], width=80), file=oH)
            logger.debug('Filtering finished')
        return target

    def pretty(self, row):
        refId = row[2]
        start = int(row[3])
        for i in row[10:]:
            if i.startswith('MD'):
                mismatchInfo = i
        sequence = row[9]
        cigar = row[5]
        cigarString = clashchimeras.methods.convertCigar(row[5])
        matchLength = cigarString.count("M") + cigarString.count("D")
        end = start + matchLength - 1
        record = Record(refId=refId, start=start, mismatchInfo=mismatchInfo,
                        sequence=sequence, cigarString=cigarString, matchLength=matchLength,
                        cigar=cigar, end=end)
        return record


class Output:
    """Contains methods for writing output files

    This class is used to generate every kind of output generated by this
    package which includes plain text, ansi colored text and bed file

    Attributes:
      target: A prefix for output file which will be automatically followed by
        extension (default 'wip')
      overlap: Minimum overlap to be set between two molecules when determining
        chimera (default 4)
      gap: Maximum gap (number of unknown nucleotides) to be allowed between
        two molecules within a chimera (default 9)
    """

    def __init__(self,
                 target=None,
                 smallRNABed=False,
                 targetRNABed=False,
                 overlap=4,
                 gap=9):
        self.target = target
        self.overlap = overlap
        self.gap = gap

        if smallRNABed:
            self.smallRNABedHandle = open('{}.smallRNA.bed'.format(self.target), 'w')
            print('# BED locations of smallRNA part of the identified chimera',
                  file=self.smallRNABedHandle)
            self.smallRNABedCSV = csv.writer(self.smallRNABedHandle, delimiter="\t")
            self.smallRNABedCSV.writerow(
                ['# The name field represents the following:'])
            self.smallRNABedCSV.writerow(
                ['#  E.g. 201980-1-48|hsa-mir-100==PAPSS1'])
            self.smallRNABedCSV.writerow(
                ['#    201980-1-48 is the fasta identifier'])
            self.smallRNABedCSV.writerow(
                ["#      201980 is the unique identifier"])
            self.smallRNABedCSV.writerow(
                ["#      1 is the number of times that sequence was observed in raw "
                 "fastq "])
            self.smallRNABedCSV.writerow(
                ["#      48 is the length of the sequence"])
            self.smallRNABedCSV.writerow(
                ['#    hsa-mir-100 represents the smallRNA transcript'])
            self.smallRNABedCSV.writerow(
                ['#    PAPSS1 represents the gene symbol for targetRNA transcript '
                 'transcript '])
        if targetRNABed:
            self.targetRNABedHandle = open('{}.targetRNA.bed'.format(self.target),
                                           'w')
            self.targetRNABedCSV = csv.writer(self.targetRNABedHandle, delimiter="\t")
            self.targetRNABedCSV.writerow(
                ['# The name field represents the following:'])
            self.targetRNABedCSV.writerow(
                ['#  E.g. 136019-1-48|ENST00000375759.6|SPEN==hsa-mir-103a-2'])
            self.targetRNABedCSV.writerow(
                ['#    136019-1-48 is the fasta identifier'])
            self.targetRNABedCSV.writerow(
                ["#      136019 is the unique identifier"])
            self.targetRNABedCSV.writerow(
                ["#      1 is the number of times that sequence was observed in raw "
                 "fastq "])
            self.targetRNABedCSV.writerow(
                ["#      48 is the length of the sequence"])
            self.targetRNABedCSV.writerow(
                ["#    ENST00000375759.6 is the targetRNA transcript identifier"])
            self.targetRNABedCSV.writerow(
                ['#    SPEN is the gene symbol for for targetRNA transcript '
                 'ENST00000375759.6'])
            self.targetRNABedCSV.writerow(
                ['#    hsa-mir-103a-2 represents the smallRNA transcript '])

        self.hybWriter = open('%s.chimeras.tsv' % self.target, 'w')
        self.hybComments()

    def hybComments(self):
        print("# fasta Identifier: The identifier in <sample>.unique.fasta. ",
              "#\tE.g. 123456-3-68 ",
              "#\t123456 is the unique identifier",
              "#\t3 is the number of times that sequence was observed in raw "
              "fastq ",
              "#\t68 is the length of the sequence", sep="\n", file=self.hybWriter)

        print("# smallRNA: The cDNA ID of the type of RNA labelled as smallRNA in "
              "the analysis",
              "#\tE.g. hsa-let-7b (miRBase identifier)",
              "#\tE.g. ENST00000619178.1|SNORD3D| (Gencode snoRNA identifier)",
              sep="\n", file=self.hybWriter)

        print("# smallRNA_start: cDNA alignment start position of the smallRNA "
              "part of the chimera", file=self.hybWriter)

        print("# smallRNA_MDtag: Showing the MD tag from the smallRNA SAM "
              "alignment for the chimera",
              "#\tSAM file format specification",
              "#\thttp://samtools.github.io/hts-specs/SAMv1.pdf",
              "#\tMD Z String for mismatching positions.Regex:[0-9]+((["
              "A-Z]|\^[A-Z]+)[0-9]+)*9", sep="\n", file=self.hybWriter)

        print('# smallRNA_cigar: Cigar string from the smallRNA SAM alignment for '
              'the chimera',
              "#\tSAM file format specification",
              "#\thttp://samtools.github.io/hts-specs/SAMv1.pdf",
              '#\tSee CIGAR in the file', sep="\n", file=self.hybWriter)

        print('# arbitrary_chimera: The chimera representation indicating what '
              'part of the sequence represents smallRNA and targetRNA',
              '#\t{ is representing a match with smallRNA',
              '#\t} is representing a match with targetRNA',
              '#\t# is representing unaligned sequences (identified as --gap -ga)',
              '#\t- is representing a deletion (D in cigar string)',
              '#\t+ is representing a deletion (I in cigar string)',
              '#\tE.g {{{{{{{{-{{{{{{{{{{{{{##}}}}}}}}}}+}}}}}}}}}}}}}}}}}}}}}}'
              '#\tE.g The first 22 nucleotides are aligning to smallRNA cDNA',
              '#\tE.g The last 33 nucleotides are aligning to targetRNA cDNA',
              sep="\n", file=self.hybWriter)

        print('# read_sequence: The actual sequence that is appeared in raw '
              'reads', file=self.hybWriter)

        print("# targetRNA: The cDNA ID of the type of RNA labelled as targetRNA "
              "in "
              "the analysis",
              "#\tE.g. hsa-let-7b (miRBase identifier)",
              "#\tE.g. ENST00000619178.1|SNORD3D| (Gencode snoRNA identifier)",
              sep="\n", file=self.hybWriter)

        print("# targetRNA_start: cDNA alignment start position of the targetRNA "
              "part of the chimera", file=self.hybWriter)

        print("# targetRNA_MDtag: Showing the MD tag from the targetRNA SAM "
              "alignment for the chimera",
              "#\tSAM file format specification",
              "#\thttp://samtools.github.io/hts-specs/SAMv1.pdf",
              "#\tMD Z String for mismatching positions.Regex:[0-9]+((["
              "A-Z]|\^[A-Z]+)[0-9]+)*9", sep="\n", file=self.hybWriter)

        print('# targetRNA_cigar: Cigar string from the targetRNA SAM alignment '
              'for '
              'the chimera',
              "#\tSAM file format specification",
              "#\thttp://samtools.github.io/hts-specs/SAMv1.pdf",
              '#\tSee CIGAR in the file', sep="\n", file=self.hybWriter)

        print("# fasta_Identifier", "smallRNA", "smallRNA_start", "smallRNA_MDtag",
              "smallRNA_cigar", "arbitrary_chimera", "read_sequence", "targetRNA",
              "targetRNA_start", "targetRNA_MDtag", "targetRNA_cigar", sep="\t",
              file=self.hybWriter)

    def writeTargetRNABed(self, query, targetRNASegments, smallRNA):
        if "ENS" in smallRNA and "|" in smallRNA:
            _smallRNA = smallRNA.split("|")[5]
        else:
            _smallRNA = smallRNA
        for segment in targetRNASegments:
            _segment = list(segment)
            _segment[3] = query + "|" + _segment[3] + "==" + _smallRNA
            self.targetRNABedCSV.writerow(_segment)

    def writeSmallRNABed(self, query, smallRNASegments, targetRNA):
        if "ENS" in targetRNA and "|" in targetRNA:
            _targetRNA = targetRNA.split("|")[5]
        else:
            _targetRNA = targetRNA
        for segment in smallRNASegments:
            _segment = list(segment)
            _segment[3] = query + "|" + _segment[3] + "==" + _targetRNA
            self.smallRNABedCSV.writerow(_segment)

    def write(self, queryName, smallRNA, targetRNA):
        chimeraString = clashchimeras.methods.chimeraOrNot(smallRNA.cigarString,
                                                           targetRNA.cigarString, overlap=self.overlap, gap=self.gap)
        smallRNARegion = clashchimeras.methods.findRegion(smallRNA)
        targetRNARegion = clashchimeras.methods.findRegion(targetRNA)
        print(queryName, smallRNARegion, smallRNA.start, smallRNA.mismatchInfo,
              smallRNA.cigar, chimeraString, smallRNA.sequence,
              targetRNARegion, targetRNA.start,
              targetRNA.mismatchInfo, targetRNA.cigar, sep="\t", file=self.hybWriter)

    def __del__(self):
        self.hybWriter.close()


class Fasta:

    def __init__(self, genome=None, gtf=None):
        self.genome = genome
        self.gtf = gtf
        self.faidx = pyfaidx.Fasta(self.genome)

    def getBiotype(self, output=None, biotype=None):
        self.sequences = []
        g = GTF(fileName=self.gtf)
        if biotype == 'tRNA':
            g.readBiotype(biotype=biotype, featureType='tRNAscan')
        else:
            g.readBiotype(biotype=biotype)
        for transcript_id, exons in g.biotypeFeatures.items():
            temp_seq = ''
            exons.sort(key=itemgetter(0))
            for exon in exons:
                if exon[4] == '-':
                    temp_seq += (-self.faidx[exon[1]][exon[2] - 1:exon[3]]).seq
                elif exon[4] == '+':
                    temp_seq += self.faidx[exon[1]][exon[2] - 1:exon[3]].seq

            _id = '{}|{}|{}|{}|{}|{}|{}'.format(transcript_id,
                                                exons[0][5],
                                                exons[0][6],
                                                exons[0][7],
                                                exons[0][8],
                                                exons[0][9],
                                                len(temp_seq))

            temp_rec = SeqRecord(seq=Seq(temp_seq), id=_id,
                                 description='')

            self.sequences.append(temp_rec)

        if not output:
            logger.error('Please provide output file..')
            sys.exit()
        else:
            logger.info('Writing {}'.format(output))
            SeqIO.write(self.sequences, output, 'fasta')


class Fastq:

    def __init__(self, fileName=None, compressed=False):
        self.fileName = fileName
        self.compressed = compressed
        self.n = 4
        self.sequences = Counter()
        self.uniqueOutput = fileName.rpartition(".")[0] + '.unique.fasta'

    def recordIterator(self):
        record = []
        record_length = 0
        for line in self.fileHandle:
            if record_length == self.n:
                yield record
                record_length = 0
                record = []
            record.append(line.decode().rstrip())
            record_length += 1
        yield record

    def createUnique(self):
        if self.compressed:
            self.fileHandle = gzip.open(self.fileName, 'rb')
        else:
            self.fileHandle = open(self.fileName, 'rb')
        logger.info('Reading {}'.format(self.fileName))
        for record in self.recordIterator():
            self.sequences[record[1]] += 1
        logger.info('Writing {}'.format(self.uniqueOutput))
        with open(self.uniqueOutput, 'w') as wH:
            for index, (sequence, counts) in enumerate(sorted(self.sequences.items(),
                                                              key=itemgetter(1), reverse=True), start=1):
                print('>{}-{}-{}'.format(index, counts, len(sequence)), file=wH)
                print(textwrap.fill(sequence, width=80), file=wH)
        logger.debug('Finished writing {}'.format(self.uniqueOutput))
        self.fileHandle.close()


class Record:
    """A custom object (preferred over dict) for easy access using variables

    It's a dependency for GTF and GFF classes
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)
