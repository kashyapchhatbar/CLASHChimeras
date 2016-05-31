import glob
import gzip
import hashlib
import logging
import os
import re
import subprocess
import sys
import textwrap
import tempfile
from collections import defaultdict

import ftputil
import pandas as pd
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from colors import green
from bs4 import BeautifulSoup
from tabulate import tabulate

from clashchimeras.parsers import Fasta

logger = logging.getLogger('root')


class Releases:

    def __init__(self, gencodeOrganism='H.sapiens', mirbaseOrganism='hsa',
                 path=os.path.join(os.path.expanduser("~"), '.CLASHchimeras')):
        organisms = {'H.sapiens': 'http://www.gencodegenes.org/releases/',
                     'M.musculus': 'http://www.gencodegenes.org/mouse_releases/'}
        shortOrganisms = {'H.sapiens': 'human',
                          'M.musculus': 'mouse'}
        self.gencodeOrganism = gencodeOrganism
        self.mirbaseOrganism = mirbaseOrganism
        self.path = path
        self.url = organisms[self.gencodeOrganism]
        self.gencodeShortOrganism = shortOrganisms[self.gencodeOrganism]
        self.mirbaseUrl = 'ftp://mirbase.org/pub/mirbase/CURRENT/README'
        self.df = self.openUrl()
        self.mirbaseDf = self.openMirbaseReadme()

        self.selectGencodeRelease()

        self.selectMirbaseRelease()

        g = Gencode(release=self.gencodeRelease,
                    path=self.gencodePath,
                    ftpDir=self.gencodeFtpDir)

        g.download()

        m = Mirbase(release=self.mirbaseRelease,
                    path=self.mirbasePath,
                    ftpDir=self.mirbaseFtpDir,
                    organism=self.mirbaseOrganism)

        m.download()

        m.process()

    def openUrl(self):
        r = requests.get(self.url)
        data = r.text
        soup = BeautifulSoup(data)
        table = soup.find("table")
        headings = [th.get_text().strip() for th in table.find("tr").find_all("th")]
        dataset = defaultdict(dict)
        for index, row in enumerate(table.find_all("tr")[1:]):
            for i, j in zip(headings,
                            (td.get_text().strip() for td in row.find_all("td"))):
                dataset[index][i] = j

        return pd.DataFrame(dataset).transpose()

    def openMirbaseReadme(self):
        with ftputil.FTPHost('mirbase.org', 'anonymous', 'anonymous') as fH:
            fobj = fH.open('pub/mirbase/CURRENT/README')
            store = False
            dataset = defaultdict(dict)
            index = 0
            for line in fobj.readlines():
                if store:
                    row = line.strip().split()
                    if len(row) == 3 and row[1][0].isdigit():
                        dataset[index]['Version'] = row[0]
                        dataset[index]['Date'] = row[1]
                        dataset[index]['Entries'] = row[2]
                        index += 1
                if 'HISTORY' in line:
                    store = True

            return pd.DataFrame(dataset).transpose()

    def downloadMirbaseRelease(self, key):
        logger.warn('Are you sure??')
        decision = input('Y/n: ')
        if decision.lower() == 'y' or decision == '':

            self.mirbaseRelease = self.mirbaseDf.loc[int(key)]['Version']

            os.makedirs(os.path.join(self.path, 'Mirbase',
                                     str(self.mirbaseRelease)), exist_ok=True)

            self.mirbasePath = os.path.join(self.path, 'Mirbase',
                                            str(self.mirbaseRelease))

            self.mirbaseFtpDir = os.path.join('pub/mirbase', self.mirbaseRelease)

            return True

    def downloadGencodeRelease(self, key):
        logger.warn('Are you sure??')
        decision = input('Y/n: ')
        if decision.lower() == 'y' or decision == '':

            self.gencodeRelease = self.df.loc[int(key)]['GENCODE release']

            os.makedirs(os.path.join(self.path, 'Gencode', self.gencodeOrganism,
                                     str(self.gencodeRelease)), exist_ok=True)

            self.gencodePath = os.path.join(self.path, 'Gencode',
                                            self.gencodeOrganism,
                                            str(self.gencodeRelease))

            self.gencodeFtpDir = os.path.join('pub/gencode/Gencode_' +
                                              self.gencodeShortOrganism,
                                              'release_' + self.gencodeRelease)

            return True

    def selectGencodeRelease(self):
        cols = self.df.columns

        logger.info("Select the release that you want \n{}".format(
                    tabulate(self.df[[cols[2], cols[0],
                                      cols[1], cols[3]]],
                             headers=['Index', cols[2], cols[0], cols[1], cols[3]],
                             tablefmt="simple")))

        logger.warn('Please bare in mind that the automatic download relies on '
                    'regex search which is known to work for Gencode release 17 '
                    'and higher ')

        releaseKeys = self.df.index

        while True:
            logger.warn('Please select which Gencode release to use (select index ' +
                        '[{}]): '.format(min(self.df.index)))
            key = input('Index: ')
            if key.isdigit() and int(key) in releaseKeys:
                logger.warn('This will download the Gencode release {} which '
                            'corresponds to Ensembl release {} and Genome assembly '
                            'version {} freezed on {}'.format(
                                self.df.loc[int(key)]['GENCODE release'],
                                self.df.loc[int(key)]['Ensembl release'],
                                self.df.loc[int(key)]['Genome assembly version'],
                                self.df.loc[int(key)]['Freeze date *']))
                if self.downloadGencodeRelease(key):
                    break

            elif key == '':
                logger.warn('This will download the latest Gencode release {} which '
                            'corresponds to Ensembl release {} and Genome assembly '
                            'version {} freezed on {}'.format(
                                self.df.loc[min(self.df.index)]['GENCODE release'],
                                self.df.loc[min(self.df.index)]['Ensembl release'],
                                self.df.loc[min(self.df.index)]['Genome assembly '
                                                                'version'],
                                self.df.loc[min(self.df.index)]['Freeze date *']))
                if self.downloadGencodeRelease(min(self.df.index)):
                    break
            else:
                logger.warn('Not a valid release index:')

    def selectMirbaseRelease(self):
        cols = self.mirbaseDf.columns

        logger.info("Select the Mirbase release you want \n{}".format(
                    tabulate(self.mirbaseDf, headers=['Index', cols[0], cols[1],
                                                      cols[2]], tablefmt="simple")))

        releaseKeys = self.mirbaseDf.index

        while True:
            logger.warn('Please select which Mirbase release to use (select index ' +
                        '[{}]): '.format(max(self.mirbaseDf.index)))
            key = input('Index: ')
            if key.isdigit() and int(key) in releaseKeys:
                logger.warn('This will download miRBase release {} which contains {} '
                            'entries and released on {}'.format(
                                self.mirbaseDf.loc[int(key)]['Version'],
                                self.mirbaseDf.loc[int(key)]['Entries'],
                                self.mirbaseDf.loc[int(key)]['Date']))
                if self.downloadMirbaseRelease(key):
                    break

            elif key == '':
                logger.warn('This will download miRBase release {} which contains {} '
                            'entries and released on {}'.format(
                                self.mirbaseDf.loc[max(self.mirbaseDf.index)]['Version'],
                                self.mirbaseDf.loc[max(self.mirbaseDf.index)]['Entries'],
                                self.mirbaseDf.loc[max(self.mirbaseDf.index)]['Date']))
                if self.downloadMirbaseRelease(max(self.mirbaseDf.index)):
                    break

            else:
                logger.warn('Not a valid release index:')


class Arguments:

    def __init__(self, args, type=None):
        self.args = args
        if type == 'align':
            if args.genomeIndex:
                self.args.genomeIndex = os.path.expanduser(args.genomeIndex)
            if args.transcriptomeIndex:
                self.args.transcriptomeIndex = os.path.expanduser(
                    args.transcriptomeIndex)
            if args.input:
                self.args.input = os.path.expanduser(args.input)
            if args.smallRNAIndex:
                self.args.smallRNAIndex = os.path.expanduser(args.smallRNAIndex)
            if self.args.targetRNAIndex:
                self.args.targetRNAIndex = os.path.expanduser(args.targetRNAIndex)

        if args.output:
            self.args.output = os.path.expanduser(args.output)

        if type == 'find':
            if args.smallRNAAnnotation:
                self.args.smallRNAAnnotation = os.path.expanduser(
                    args.smallRNAAnnotation)
            if args.targetRNAAnnotation:
                self.args.targetRNAAnnotation = os.path.expanduser(
                    args.targetRNAAnnotation)
            if args.smallRNA:
                self.args.smallRNA = os.path.expanduser(args.smallRNA)
            if args.targetRNA:
                self.args.targetRNA = os.path.expanduser(args.targetRNA)

        self.exit = False

    def hash(self, file):
        sha256 = hashlib.sha256()
        with open(file, 'r+b') as f:
            while True:
                buf = f.read(8192)
                if not buf:
                    break
                sha256.update(buf)
        return sha256.hexdigest()

    def validateFind(self):
        if not os.path.exists(self.args.smallRNA):
            logger.error('{} not found...'.format(self.args.smallRNA))
            self.exit = True
        if not os.path.exists(self.args.targetRNA):
            logger.error('{} not found...'.format(self.args.smallRNA))
            self.exit = True
        if self.args.getGenomicLocationsSmallRNA:
            if not self.args.smallRNAAnnotation:
                logger.error('--smallRNAAnnotation -sa not specified. If you want '
                             'genomic locations for smallRNA, please enable '
                             '--getGenomicLocationsSmallRNA -ggs and specify '
                             'annotation file using --smallRNAAnnotation -sa')
                self.exit = True
            else:
                if not os.path.exists(self.args.smallRNAAnnotation):
                    logger.error('{} not found'.format(self.args.smallRNAAnnotation))
                    self.exit = True
        if self.args.getGenomicLocationsTargetRNA:
            if not self.args.targetRNAAnnotation:
                logger.error('--targetRNAAnnotation -ta not specified. If you want '
                             'genomic locations for targetRNA, please enable '
                             '--getGenomicLocationsTargetRNA -ggt and specify '
                             'annotation file using --targetRNAAnnotation -ta')
                self.exit = True
            else:
                if not os.path.exists(self.args.targetRNAAnnotation):
                    logger.error('{} not found...'.format(self.args.targetRNAAnnotation))
                    self.exit = True
        if self.hash(self.args.smallRNA) == self.hash(self.args.targetRNA):
            logger.error('CLASHChimeras does not detect chimeras between the same '
                         'RNA type yet... Please hang in there, we are planning it '
                         'in the feature')
            self.exit = True

        if self.exit:
            logger.error('Exiting because of the above errors...')
            sys.exit()
        elif self.hash(self.args.smallRNA) == self.hash(self.args.targetRNA):
            logger.error('CLASHChimeras does not detect chimeras between the same '
                         'RNA type yet... Please hang in there, we are planning it '
                         'in the feature')
            sys.exit()

    def validateAlign(self):
        if not os.path.exists(self.args.input):
            logger.error('{} not found...'.format(self.args.input))
            self.exit = True
        if self.args.run == 'bowtie2' and (not self.args.smallRNAIndex or
                                           not self.args.targetRNAIndex):
            logger.error('{}'.format('Please specify --smallRNAIndex -si and ' +
                                     '--targetRNAIndex -ti properly...'))
            self.exit = True
        if self.args.smallRNAIndex:
            bts = glob.glob(self.args.smallRNAIndex + '*.bt2')
            if not len(bts) >= 6:
                logger.error('Something wrong with the {}'.format(
                    self.args.smallRNAIndex))
                logger.error("Can't find all the .bt2 files for that index")
                self.exit = True
        if self.args.targetRNAIndex:
            bts = glob.glob(self.args.targetRNAIndex + '*.bt2')
            if not len(bts) >= 6:
                logger.error('Something wrong with the {}'.format(
                    self.args.targetRNAIndex))
                logger.error("Can't find all the .bt2 files for that index")
                self.exit = True
        if self.args.genomeIndex:
            bts = glob.glob(self.args.genomeIndex.strip() + '*.bt2')
            if not len(bts) >= 6:
                logger.error('Something wrong with the {}'.format(
                    self.args.genomeIndex))
                logger.error("Can't find all the .bt2 files for that index")
                self.exit = True
        if self.args.run == 'tophat' and not self.args.genomeIndex:
            logger.error('Please provide genome index if you want to use tophat..')
            self.exit = True
        if self.exit:
            sys.exit()


class Gencode:
    """This class makes sure the required files (sequences, annotation etc.) are
    present and in working order. It generates sha256 checksums and autodownloads
    the required files from Gencode Genes.
    """

    def __init__(self, release=None, user="anonymous", password=None,
                 path=None, ftpDir=None):

        self.downloaded = []
        self.release = release
        self.reList = ['^GRC.+\.genome\.fa',
                       '^gencode.+chr_patch_hapl_scaff\.annotation\.gtf',
                       '^gencode.+pc_transcripts\.fa',
                       '^gencode.+tRNAs\.gtf',
                       '^gencode.+lncRNA_transcripts\.fa']

        self.host = "ftp.sanger.ac.uk"
        self.user = user
        self.path = path

        self.files = []
        self.otherFiles = {}
        biotypes = ['snRNA', 'snoRNA', 'misc_RNA', 'tRNA']
        for b in biotypes:
            self.otherFiles[b] = 'gencode.v{}.{}_transcripts.fa'.format(self.release,
                                                                        b)
        if self.user == "anonymous":
            self.password = "anonymous"
        else:
            self.password = password

        ftp_host = ftputil.FTPHost(self.host, self.user, self.password)

        self.dir = ftpDir
        fileList = ftp_host.listdir(self.dir)
        for file in fileList:
            for rEx in self.reList:
                if re.match(rEx, file):
                    if not 'primary_assembly' in file:
                        self.files.append(file)
        ftp_host.close()
        _files = self.files.copy()
        _otherFiles = self.otherFiles.copy()

        for file in _files:
            logger.debug('Checking %s' % file)
            _file = file.rpartition(".")[0]
            if os.path.exists(os.path.join(self.path, _file + '.sha256sum')) and \
                    os.path.exists(os.path.join(self.path, _file)):
                with open(os.path.join(self.path, _file + '.sha256sum')) as f:
                    s = f.readline().rstrip()
                _s = self.hash(os.path.join(self.path, _file))
                if s == _s:
                    logger.info("%s is present and verified" % _file)
                    self.files.remove(file)
                else:
                    logger.warn('%s is downloaded but doesnt match the sha256sum' % _file)
                    logger.error('Will be downloaded again')
            else:
                logger.warn('%s will be downloaded' % _file)

        for biotype, file in _otherFiles.items():
            logger.debug('Checking %s' % file)
            _file = file
            if os.path.exists(os.path.join(self.path, _file + '.sha256sum')) and \
                    os.path.exists(os.path.join(self.path, _file)):
                with open(os.path.join(self.path, _file + '.sha256sum')) as f:
                    s = f.readline().rstrip()
                _s = self.hash(os.path.join(self.path, _file))
                if s == _s:
                    logger.info("%s is present and verified" % file)
                    self.otherFiles.pop(biotype)
                else:
                    logger.warn('%s is generated but doesnt match the sha256sum' % file)
                    logger.error('Will be generated again')
            else:
                logger.warn('%s will be generated' % file)

    def download(self):

        if len(self.files) == 0:
            logger.info('Gencode files are downloaded and checksum verified...')
        else:
            with ftputil.FTPHost(self.host, self.user, self.password) as fH:
                for file in self.files:
                    logger.info('Downloading %s from ftp.sanger.ac.uk' % file)
                    fH.download(self.dir + '/' + file, os.path.join(self.path, file),
                                callback=None)
                    self.downloaded.append(file)
                    _file = file.rpartition(".")[0]
                    p = subprocess.Popen(['gunzip', os.path.join(self.path, file)])
                    p.communicate()
                    sha256sum = self.hash(os.path.join(self.path, _file))
                    with open(os.path.join(self.path, _file + '.sha256sum'), 'w') as wH:
                        print(sha256sum, file=wH)
                    logger.info('Downloading, extraction and hashing of %s finished' %
                                file)

        gtfFiles = glob.glob(os.path.join(self.path, "*.annotation.gtf"))

        if len(gtfFiles) == 0:
            logger.warn('This release does not contain annotation file or they are '
                        'not grabbed by regex. Please download it manually from {'
                        '}'.format(self.dir))
            gtfFile = None
        elif len(gtfFiles) == 2:
            gtfFiles.sort()
            gtfFile = gtfFiles[1]
        elif len(gtfFiles) == 1:
            gtfFile = gtfFiles[0]

        genomeFiles = glob.glob(os.path.join(self.path, "*.genome.fa"))

        if len(genomeFiles) == 0:
            logger.warn('This release does not contain genome file or they are not '
                        'grabbed by regex. Please download it manually from {'
                        '}'.format(self.dir))
            genomeFile = None
        else:
            genomeFile = genomeFiles[0]

        tRNAgtfFiles = glob.glob(os.path.join(self.path, "*.tRNAs.gtf"))

        if len(tRNAgtfFiles) == 0:
            logger.warn(('This release does not contain tRNA annotation file or they '
                         'are not grabbed by regex. Please download it manually from {'
                         '}'.format(self.dir)))
            tRNAgtfFile = None
        else:
            tRNAgtfFile = tRNAgtfFiles[0]

        if len(self.otherFiles) == 0:
            logger.info('Other fasta files are generated and checksum verified...')
        else:
            for biotype, file in self.otherFiles.items():
                if biotype == 'tRNA' and (tRNAgtfFile and genomeFile):
                    fasta = Fasta(genome=genomeFile, gtf=tRNAgtfFile)
                    fasta.getBiotype(biotype='tRNA', output=os.path.join(self.path,
                                                                         file))
                    sha256sum = self.hash(os.path.join(self.path, file))
                    with open(os.path.join(self.path, file + '.sha256sum'), 'w') as wH:
                        print(sha256sum, file=wH)
                    logger.info('Extraction and hashing of %s finished' % file)
                elif gtfFile and genomeFile:
                    fasta = Fasta(genome=genomeFile, gtf=gtfFile)
                    fasta.getBiotype(biotype=biotype, output=os.path.join(self.path,
                                                                          file))
                    sha256sum = self.hash(os.path.join(self.path, file))
                    with open(os.path.join(self.path, file + '.sha256sum'), 'w') as wH:
                        print(sha256sum, file=wH)
                    logger.info('Extraction and hashing of %s finished' % file)

    def hash(self, file):
        sha256 = hashlib.sha256()
        with open(file, 'r+b') as f:
            while True:
                buf = f.read(8192)
                if not buf:
                    break
                sha256.update(buf)
        return sha256.hexdigest()


class Mirbase:
    """This class makes sure that requried files (sequences, annotations etc.)
    from Mirbase are downloaded. It generates sha256 checksums and autodownloads
    the required files
    """

    def __init__(self, release=None, user="anonymous", password=None,
                 path=os.path.join(os.path.expanduser("~"), '.CLASHchimeras'),
                 organism='hsa',
                 ftpDir=None):
        self.reList = ['^hairpin\.fa\.gz', '^mature\.fa\.gz']
        organisms = {}
        with ftputil.FTPHost('mirbase.org', 'anonymous', 'anonymous') as fH:
            tH, tP = tempfile.mkstemp()
            fH.download(os.path.join('pub/mirbase', release, 'organisms.txt.gz'), tP)
            for i in gzip.open(tP):
                ii = i.decode('utf-8').strip()
                if not ii.startswith('#'):
                    organisms[ii.split("\t")[0]] = ii.split("\t")[2]

        self.shortOrganism = organism
        self.organism = organisms.get(self.shortOrganism, None)
        if not self.organism:
            logger.error('{} is not a valid organism name in Mirbase'.format(
                self.shortOrganism))
            sys.exit()
        self.host = 'mirbase.org'
        self.release = release
        self.user = user
        self.path = path
        self.files = []

        if self.user == "anonymous":
            self.password = "anonymous"
        else:
            self.password = password

        ftp_host = ftputil.FTPHost(self.host, self.user, self.password)

        self.dir = ftpDir

        fileList = ftp_host.listdir(self.dir)

        self.files.append('genomes/' + self.shortOrganism + '.gff3')
        for file in fileList:
            for rEx in self.reList:
                if re.match(rEx, file):
                    self.files.append(file)
        ftp_host.close()
        _files = self.files.copy()

        for file in _files:
            logger.debug('Checking %s' % file)
            if file.endswith('gz'):
                _file = file.rpartition(".")[0]
            elif '/' in file:
                _file = file.rpartition("/")[2]
            if os.path.exists(os.path.join(self.path, _file + '.sha256sum')) and \
                    os.path.exists(os.path.join(self.path, _file)):
                with open(os.path.join(self.path, _file + '.sha256sum')) as f:
                    s = f.readline().rstrip()
                _s = self.hash(os.path.join(self.path, _file))
                if s == _s:
                    logger.info("%s is present and verified" % file)
                    self.files.remove(file)
                else:
                    logger.warn('%s is downloaded but doesnt match the sha256sum' % file)
                    logger.error('Must run download method')
            else:
                logger.warn('%s will be downloaded when you run download method' % file)

    def download(self):

        if len(self.files) == 0:
            logger.debug('Mirbase files are downloaded and checksum verified...')
        else:
            with ftputil.FTPHost(self.host, self.user, self.password) as fH:
                for file in self.files:
                    logger.info('Downloading %s from mirbase.org' % file)
                    if "/" in file:
                        _file = file.rpartition("/")[2]
                    else:
                        _file = file
                    fH.download(self.dir + '/' + file, os.path.join(self.path, _file),
                                callback=None)
                    if _file.endswith('.gz'):
                        __file = _file.rpartition(".")[0]
                        p = subprocess.Popen(['gunzip', os.path.join(self.path, __file)])
                        p.communicate()
                    else:
                        __file = _file
                    sha256sum = self.hash(os.path.join(self.path, __file))
                    with open(os.path.join(self.path, __file + '.sha256sum'), 'w') as wH:
                        print(sha256sum, file=wH)
                    logger.debug('Downloading, extraction and hashing of %s finished' %
                                 file)
            self.files = []

    def hash(self, file):
        sha256 = hashlib.sha256()
        with open(file, 'r+b') as f:
            while True:
                buf = f.read(8192)
                if not buf:
                    break
                sha256.update(buf)
        return sha256.hexdigest()

    def process(self):
        files = ['hairpin.fa', 'mature.fa']
        for file in files:
            temp = []
            if os.path.exists(os.path.join(self.path, self.shortOrganism + '-' + file)) \
                and os.path.exists(os.path.join(self.path,
                                                self.shortOrganism + '-' + file + '.sha256sum')):
                with open(os.path.join(self.path, file + '.sha256sum')) as f:
                    s = f.readline().rstrip()
                _s = self.hash(os.path.join(self.path, file))
                if s == _s:
                    logger.info("%s is present and verified" % file)
                else:
                    logger.warn('%s is downloaded but doesnt match the sha256sum' % file)
            else:
                logger.debug('Extrating %s sequences from %s' % (self.organism, file))
                with open(os.path.join(self.path, file)) as iH:
                    for rec in SeqIO.parse(iH, 'fasta'):
                        if self.organism in rec.description:
                            _temp_rec = SeqRecord(id=rec.id, description=rec.description,
                                                  seq=rec.seq.back_transcribe())
                            temp.append(_temp_rec)
                SeqIO.write(temp, os.path.join(self.path, self.shortOrganism + '-' + file),
                            'fasta')
                with open(os.path.join(self.path,
                                       self.shortOrganism + '-' + file + '.sha256sum'), 'w') as wH:
                    print(self.hash(os.path.join(self.path,
                                                 self.shortOrganism + '-' + file + '.sha256sum')), file=wH)


class Index:

    def __init__(self, root=os.path.join(os.path.expanduser("~"),
                                         '.CLASHchimeras'), bowtieExecutable=None, tophatExecutable=None):

        self.created = []

        if bowtieExecutable is not None:
            logger.info('Bowtie2 path {} provided and will be used'
                        .format(bowtieExecutable))
            self.bowtieExecutable = bowtieExecutable
        else:
            try:
                bowtieVersion = float(subprocess.check_output(['bowtie2', '--version']
                                                              ).decode('utf-8').split("\n")[0].split()[-1].rpartition(".")[0])
                fullBowtieVersion = subprocess.check_output(['bowtie2', '--version']
                                                            ).decode('utf-8').split("\n")[0].split()[-1]
                if bowtieVersion < 2.2:
                    logger.warn("Please update bowtie2, you can download it from")
                    logger.warn("http://sourceforge.net/projects/bowtie-bio/files/bowtie2/")
                else:
                    logger.info("Bowtie2 version {} found and will be used"
                                .format(fullBowtieVersion))
                    self.bowtieExecutable = 'bowtie2'
            except OSError as e:
                if e.errno == os.errno.ENOENT:
                    logger.error(textwrap.fill("""
            Can't find Bowtie2 in your path. Please install it from
            http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
            """))
                else:
                    logger.error('Something wrong with your bowtie2 command')

        if tophatExecutable is not None:
            logger.info('Tophat2 path {} provided and will be used'
                        .format(tophatExecutable))
            self.tophatExecutable = tophatExecutable
        else:
            try:
                tophatVersion = float(subprocess.check_output(['tophat',
                                                               '--version']).decode('utf-8').split(" v")[-1].rpartition(".")[0])
                fullTophatVersion = subprocess.check_output(['tophat',
                                                             '--version']).decode('utf-8').split(" v")[-1].rstrip()
                if tophatVersion < 2:
                    logger.warn("Please update tophat, you can download it from")
                    logger.warn("http://ccb.jhu.edu/software/tophat/index.shtml")
                else:
                    logger.info("Tophat version {} found and will be used"
                                .format(fullTophatVersion))
                    self.tophatExecutable = 'tophat'
            except OSError as e:
                if e.errno == os.errno.ENOENT:
                    logger.error("Can't find Tophat in your path. Please install it from")
                    logger.error("http://ccb.jhu.edu/software/tophat/index.shtml")
                else:
                    logger.error('Something wrong with your tophat command')

        self.root = root

    def create(self):

        stdoutHandle = open(os.path.join(self.root, 'CLASHChimeras-Index.stdout'),
                            'w')
        stderrHandle = open(os.path.join(self.root, 'CLASHChimeras-Index.stderr'),
                            'w')
        indexCreated = False
        for fa in glob.glob(os.path.join(self.root, '*.fa')):
            indexName = fa.rpartition(".fa")[0]
            self.created.append('{} - {}'.format(indexName, green('Bowtie2',
                                                                  style='bold')))
            if 'genome.fa' in fa:
                genomeIndex = fa.rpartition(".")[0]
            if len(glob.glob(os.path.join(self.root, indexName + '*.bt2'))) == 6:
                logger.info('Bowtie2 index for {} found'.format(fa))
            else:
                try:
                    logger.info('Generating bowtie2 index for {}'.format(fa))
                    indexCreated = True
                    p = subprocess.Popen([self.bowtieExecutable + '-build', fa, indexName],
                                         stdout=stdoutHandle, stderr=stderrHandle)
                    p.communicate()
                    logger.info('Bowtie2 index for {} is generated'.format(fa))
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        logger.error("Can't find bowtie2-build in your path...")
                        logger.error("Please make sure bowtie2 is in your path and")
                        logger.error("bowtie2-build can run from your shell")
                    else:
                        logger.error("Something wrong with your bowtie2-build command")

        if indexCreated:
            logger.warn('The stdout and stderr for bowtie2-build are stored in')
            logger.warn('CLASHChimeras-Index.stdout and CLASHChimeras-Index.stderr')

        for gtf in glob.glob(os.path.join(self.root, '*.annotation.gtf')):
            indexName = gtf.rpartition("/")[-1].rpartition(".")[0]
            self.created.append('{} - {}'.format(gtf.rpartition(".")[0], green(
                                'Tophat2 Transcriptome', style='bold')))
            if len(glob.glob(os.path.join(self.root, indexName + '*.bt2'))) == 6:
                logger.info('Tophat2 transcriptome index found for {}'.format(gtf))
            else:
                try:
                    indexCreated = True
                    logger.info('Generating tophat transcriptome index for {}'
                                .format(gtf))
                    p = subprocess.Popen([self.tophatExecutable, '-G', gtf,
                                          '--transcriptome-index', os.path.join(
                                              self.root, indexName),
                                          genomeIndex], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    p.communicate()
                    logger.info('Tophat2 transcriptome index generated for {}'
                                .format(gtf))
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        logger.error("Can't find tophat in your path...")
                        logger.error("Please make sure tophat is in your path")
                    else:
                        logger.error("Something wrong with your tophat command")

        stdoutHandle.close()
        stderrHandle.close()
