import hashlib
import logging
import os
import subprocess

from clashchimeras.parsers import Fastq

logger = logging.getLogger('root')


class Bowtie:
    """Class for calling bowtie2 aligner
    """

    def __init__(self, args=None):
        self.parameters = []
        self.inputFile = args.input
        self.compressed = args.gzip
        self.bowtieExecutable = args.bowtieExecutable
        self.uniqueInput = self.inputFile.rpartition(".")[0] + '.unique.fasta'
        if os.path.exists(self.uniqueInput) and os.path.exists(
                self.uniqueInput + '.sha256sum'):
            with open(self.uniqueInput + '.sha256sum') as f:
                s = f.readline().rstrip()
            _s = self.hash(self.uniqueInput)
            if s == _s:
                logger.info(
                    '{} is present and verified'.format(self.uniqueInput))
            else:
                logger.warn('{} is present but checksum verification '
                            'failed'.format(self.uniqueInput))
                logger.warn(
                    '{} will be generated again'.format(self.uniqueInput))
        else:
            _f = Fastq(fileName=self.inputFile, compressed=self.compressed)
            _f.createUnique()
            sha256sum = self.hash(self.uniqueInput)
            with open(self.uniqueInput + '.sha256sum', 'w') as wH:
                print(sha256sum, file=wH)
        if not args.reverseComplement:
            self.parameters.append('--norc')
        self.parameters.append('-N')
        self.parameters.append(str(args.mismatch))
        if not args.unaligned:
            self.parameters.append('--no-unal')
        self.parameters.append('-p')
        self.parameters.append(str(args.threads))
        self.parameters.append('--' + args.preset)
        self.smallRNAIndex = args.smallRNAIndex
        self.targetRNAIndex = args.targetRNAIndex
        self.output = args.output
        self.runAligner = True

    def hash(self, file):
        sha256 = hashlib.sha256()
        with open(file, 'r+b') as f:
            while True:
                buf = f.read(8192)
                if not buf:
                    break
                sha256.update(buf)
        return sha256.hexdigest()

    def run(self, output=None, filtered=None):
        self.bowtie2 = []
        if output == 'smallRNA':
            self.outputSam = '{}.{}.sam'.format(self.output, 'smallRNA')
        elif output == 'targetRNA':
            self.outputSam = '{}.{}.sam'.format(self.output, 'targetRNA')
        if not self.bowtieExecutable:
            self.bowtie2.append('bowtie2')
        else:
            self.bowtie2.append(self.bowtieExecutable)
        self.bowtie2.append('-f')
        if filtered:
            self.bowtie2.append(filtered)
        else:
            self.bowtie2.append(self.uniqueInput)
        self.bowtie2.append('-x')
        if output == 'smallRNA':
            self.bowtie2.append('"' + self.smallRNAIndex + '"')
        elif output == 'targetRNA':
            self.bowtie2.append('"' + self.targetRNAIndex + '"')
        self.bowtie2 = self.bowtie2 + self.parameters
        self.bowtie2.append('-S')
        self.bowtie2.append(self.outputSam)

        if os.path.exists(self.outputSam) and os.path.exists(
                self.outputSam + '.sha256sum'):
            with open(self.outputSam + '.sha256sum') as f:
                s = f.readline().rstrip()
            _s = self.hash(self.outputSam)
            if s == _s:
                logger.info(
                    '{} is present and verified'.format(self.outputSam))
                self.runAligner = False
            else:
                logger.warn('{} is present but checksum verification '
                            'failed'.format(self.outputSam))
                logger.warn('{} will be generated again')
                self.runAligner = True

        if self.runAligner:
            logger.info(
                'Executing bowtie2 and generating {}'.format(self.outputSam))
            logger.info('The complete command for bowtie2 is \n{}'.format(
                        " ".join(self.bowtie2)))
            bowtie2Process = subprocess.Popen(self.bowtie2,
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
            bowtie2Stdout, bowtie2Stderr = bowtie2Process.communicate()
            logger.warn(bowtie2Stdout.decode('utf-8'))
            logger.warn(bowtie2Stderr.decode('utf-8'))
            if os.path.exists(self.outputSam):
                logger.info('Bowtie2 execution finished and {} '
                            'generated'.format(self.outputSam))
                sha256sum = self.hash(self.outputSam)
                with open(self.outputSam + '.sha256sum', 'w') as wH:
                    print(sha256sum, file=wH)
            else:
                logger.error('Something went wrong with your bowtie2 command.')
                logger.error(
                    'Please check your bowtie2 stderr and report a bug....')


class Tophat:

    def __init__(self, args=None):
        self.inputFile = args.input
        if args.genomeIndex:
            args.genomeIndex = os.path.expanduser(args.genomeIndex)
        if args.transcriptomeIndex:
            args.transcriptomeIndex = os.path.expanduser(
                args.transcriptomeIndex)
        if args.output:
            args.output = os.path.expanduser(args.output)
        if args.input:
            args.input = os.path.expanduser(args.input)

        self.tophatOutput = args.output
        self.genomeIndex = args.genomeIndex

        self.tophatRun = []
        if args.tophatExecutable is not None:
            self.tophatExecutable = args.tophatExecutable
        else:
            self.tophatExecutable = 'tophat'
        self.tophatRun.append(self.tophatExecutable)
        self.tophatRun.append('--b2-' + args.tophatPreset)
        self.tophatRun.append('--b2-N')
        self.tophatRun.append(str(args.mismatch))

        if args.transcriptomeIndex:
            self.tophatRun.append('--transcriptome-index=' +
                                  args.transcriptomeIndex.strip())
        else:
            logger.warn('Transcriptome index is not provided for tophat')
        self.tophatRun.append('-p')
        self.tophatRun.append(str(args.threads))
        self.tophatRun.append('-o')
        self.tophatRun.append(self.tophatOutput)
        self.tophatRun.append(self.genomeIndex)
        self.tophatRun.append(self.inputFile)

    def run(self):
        logger.info('Executing Tophat2 for {}'.format(self.inputFile))
        logger.info('The complete command for tophat is \n{}'.format(" ".join(
            self.tophatRun)))
        logger.info('The stdout and stderr for tophat run will be '
                    'dumped after its finished')
        tophatProcess = subprocess.Popen(self.tophatRun,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
        tophatStderr, tophatStdout = tophatProcess.communicate()
        logger.warn('Tophat stderr')
        logger.info(tophatStderr.decode("utf-8"))
        logger.warn('Tophat stdout')
        logger.info(tophatStdout.decode("utf-8"))
        logger.info('Execution of Tophat2 finished')
