import argparse
import textwrap
import sys

from colors import yellow

import clashchimeras.log
import clashchimeras.methods
from clashchimeras.parsers import GFF, GTF, Output, SAM
from clashchimeras.initialize import Arguments


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


def parseArguments():
    requiredArgs = getRequiredArgs()
    optionalArgs = getOptionalArgs()
    outputArgs = getOutputArgs()

    parser = argparse.ArgumentParser(
        parents=[requiredArgs, outputArgs, optionalArgs],
        formatter_class=CustomFormatter,

        usage='An example usage is: %(prog)s -s smallRNA.sam -t '
        'targetRNA.sam -o output',

        description=textwrap.dedent('''\
                Given two SAM files, this script tries to find chimeras that
                are observed between a smallRNA and a targetRNA'''),

        add_help=True)

    return parser


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    required.add_argument('--smallRNA', '-s',
                          help='Provide smallRNA alignment SAM file',
                          metavar='SAM file',
                          required=True)
    required.add_argument('--targetRNA', '-t',
                          help='Provide targetRNA alignment SAM file',
                          metavar='SAM file',
                          required=True)

    required.add_argument("--getGenomicLocationsSmallRNA", "-ggs",
                          help="Do you want genomic locations for small RNA?",
                          action="store_true")

    required.add_argument("--getGenomicLocationsTargetRNA", "-ggt",
                          help="Do you want genomic locations for target RNA?",
                          action="store_true")

    required.add_argument("--smallRNAAnnotation", "-sa",
                          action="store",
                          help="""Provide smallRNA annotation gtf(from Gencode)
                        or gff3(from Mirbase). Only provide gtf from Gencode or
                        gff3 from Mirbase. Does not support other gtf files""")

    required.add_argument("--targetRNAAnnotation", "-ta",
                          action="store",
                          help="""Provide targetRNA annotation gtf(from Gencode).
                        Only provide gtf from Gencode. Does not support other
                        gtf files""")

    return parser


def getOutputArgs():
    parser = argparse.ArgumentParser(add_help=False)

    output = parser.add_argument_group('Output')

    output.add_argument("--output", "-o",
                        help="""The output name without extension (.bed .csv will
                      be added)""",
                        metavar='output prefix',
                        required=True)

    return parser


def getOptionalArgs():
    parser = argparse.ArgumentParser(add_help=False)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--overlap", "-ov",
                          help="""Maximum overlap to be set between two
                        molecules when determining chimeras""",
                          default=4)

    optional.add_argument("--gap", "-ga",
                          help="""Maximum gap (number of unaligned nucleotides) to
                        allowed between two molecules within a chimera""",
                          default=9)

    optional.add_argument("--logLevel", "-l",
                          help="Set logging level",
                          default='INFO',
                          choices=['INFO', 'DEBUG', 'WARNING', 'ERROR'])

    return parser


def main():
    """The main function to execute when you run the program as a script
    """
    parser = parseArguments()

    args = parser.parse_args()

    if args.logLevel == 'DEBUG':
        logger = clashchimeras.log.debug_logger('root')
    elif args.logLevel == 'WARNING':
        logger = clashchimeras.log.warning_logger('root')
    elif args.logLevel == 'ERROR':
        logger = clashchimeras.log.error_logger('root')
    else:
        logger = clashchimeras.log.info_logger('root')

    argCheck = Arguments(args, type='find')
    argCheck.validateFind()

    smallRNA = SAM(fileName=args.smallRNA)
    smallRNA.read()

    targetRNA = SAM(fileName=args.targetRNA)
    targetRNA.read()

    if args.getGenomicLocationsSmallRNA:

        if args.smallRNAAnnotation.endswith('gff3'):

            smallRNAAnnotation = GFF(fileName=args.smallRNAAnnotation)
            smallRNAAnnotation.read()

        elif args.smallRNAAnnotation.endswith('gtf'):
            smallRNAAnnotation = GTF(fileName=args.smallRNAAnnotation)
            if 'tRNAs' in args.smallRNAAnnotation:
                smallRNAAnnotation.read(featureType='tRNAscan')
            else:
                smallRNAAnnotation.read()

    if args.getGenomicLocationsTargetRNA:
        if args.targetRNAAnnotation.endswith('gtf'):
            targetRNAAnnotation = GTF(fileName=args.targetRNAAnnotation)
            targetRNAAnnotation.read()

    output = Output(target=args.output,
                    smallRNABed=args.getGenomicLocationsSmallRNA,
                    targetRNABed=args.getGenomicLocationsTargetRNA,
                    overlap=args.overlap,
                    gap=args.gap)

    common = set(smallRNA.records.keys()).intersection(targetRNA.records.keys())

    logger.info('Determining chimeras')

    for index, queryName in enumerate(common, start=1):
        smallRNARecord = smallRNA.access(queryName)
        targetRNARecord = targetRNA.access(queryName)

        chimeraString = clashchimeras.methods.chimeraOrNot(
            smallRNARecord.cigarString, targetRNARecord.cigarString,
            overlap=args.overlap, gap=args.gap)

        if chimeraString:

            if args.getGenomicLocationsSmallRNA:

                try:
                    smallRNALocation = smallRNAAnnotation.coordinates(
                        smallRNARecord.refId,
                        start=smallRNARecord.start,
                        end=smallRNARecord.end)
                except IndexError as e:
                    logger.error(e)
                    logger.error('Are you sure you have the correct annotation file for '
                                 'smallRNA?')
                    logger.error('Please check --smallRNAAnnotation -sa and '
                                 'and make sure they are matching with the alignment '
                                 'file {}'.format(args.smallRNA))
                    sys.exit()
                output.writeSmallRNABed(queryName, smallRNALocation,
                                        targetRNARecord.refId)

            if args.getGenomicLocationsTargetRNA:
                try:
                    targetRNALocation = targetRNAAnnotation.coordinates(
                        targetRNARecord.refId,
                        start=targetRNARecord.start,
                        end=targetRNARecord.end)
                except IndexError as e:
                    logger.error(e.msg)
                    logger.error('Are you sure you have the correct annotation file for '
                                 'targetRNA?')
                    logger.error('Please check --targetRNAAnnotation -ta and '
                                 'and make sure they are matching with the alignment '
                                 'file {}'.format(args.targetRNA))
                    sys.exit()

                output.writeTargetRNABed(queryName, targetRNALocation,
                                         smallRNARecord.refId)

            output.write(queryName, smallRNARecord, targetRNARecord)

    logger.info('Determining chimeras finished')

if __name__ == "__main__":
    main()
