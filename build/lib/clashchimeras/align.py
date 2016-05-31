import argparse
import sys

import clashchimeras.log
from clashchimeras.initialize import Arguments
from clashchimeras.parsers import SAM
from clashchimeras.runners import Bowtie, Tophat


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


def parseArguments():
    requiredArgs = getRequiredArgs()
    bowtieArgs = getBowtie2Args()
    tophatArgs = getTophatArgs()
    optionalArgs = getOptionalArgs()
    outputArgs = getOutputArgs()

    parser = argparse.ArgumentParser(
        parents=[requiredArgs, bowtieArgs, tophatArgs,
                 outputArgs, optionalArgs],
        formatter_class=CustomFormatter,
        description='Given a fastq file, this script executes '
        'bowtie2 and tophat aligners to generate alignment files '
        'necessary for detecting chimeras in the reads',
        usage='\n %(prog)s -i input.fastq -si '
        '/path/to/smallRNA_index -ti /path/to/targetRNA_index -o '
        'output -r bowtie2 \n %(prog)s -i input.fastq -gi '
        '/path/to/genome_index -tri /path/to/transcriptome_index '
        '-o output -r tophat \n \n \n '
        'To see detailed help, please run \n %(prog)s -h',
        add_help=True)

    return parser

    args = parser.parse_args()

    if args.logLevel == 'DEBUG':
        logger = clashchimeras.log.debug_logger('root')
    elif args.logLevel == 'WARNING':
        logger = clashchimeras.log.warning_logger('root')
    elif args.logLevel == 'ERROR':
        logger = clashchimeras.log.error_logger('root')
    else:
        logger = clashchimeras.log.info_logger('root')

    argCheck = Arguments(args, type='align')
    argCheck.validateAlign()

    return args


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Input arguments')

    required.add_argument('--input', '-i',
                          help='Input file containing reads fastq',
                          metavar='raw reads',
                          required=True)

    return parser


def getBowtie2Args():
    parser = argparse.ArgumentParser(add_help=False)

    bowtieArgs = parser.add_argument_group('Bowtie2 arguments')

    bowtieArgs.add_argument("--smallRNAIndex", "-si",
                            help="""Provide the smallRNA bowtie2 index (Usually
                        resides in ~/db/CLASHChimeras or elsewhere if you have
                        specified in --path -pa during initialize)""")

    bowtieArgs.add_argument("--targetRNAIndex", "-ti",
                            help="""Provide the targetRNA bowtie2 index (Usually
                        resides in ~/db/CLASHChimeras or elsewhere if you have
                        specified in --path -pa during initialize)""")

    return parser


def getTophatArgs():
    parser = argparse.ArgumentParser(add_help=False)

    tophatArgs = parser.add_argument_group('Tophat arguments')

    tophatArgs.add_argument("--genomeIndex", "-gi",
                            help="""Provide the genome bowtie2 index (Usually
                        resides in ~/db/CLASHChimeras or elsewhere if you have
                        specified in --path during initialize)""")

    tophatArgs.add_argument("--transcriptomeIndex", "-tri",
                            help="""Provide the transcriptome index as specified
                        in tophat --transcriptome-index""")

    return parser


def getOutputArgs():
    parser = argparse.ArgumentParser(add_help=False)

    output = parser.add_argument_group('Output')

    output.add_argument("--output", "-o",
                        help="""The output name without extension (.sam .bam will
                      be added)""",
                        metavar='output prefix',
                        required=True)

    return parser


def getOptionalArgs():
    parser = argparse.ArgumentParser(add_help=False)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--run", "-r",
                          help="Run the following aligner for raw reads",
                          default='bowtie2',
                          choices=['bowtie2', 'tophat'])

    optional.add_argument("--logLevel", "-l",
                          help="Set logging level",
                          default='INFO',
                          choices=['INFO', 'DEBUG', 'WARNING', 'ERROR'])

    optional.add_argument("--gzip", "-gz", action="store_true",
                          help="Whether your input file is gzipped")

    optional.add_argument("--bowtieExecutable", "-be",
                          help="""Provide bowtie2 executable if it's not present
                        in your path""")

    optional.add_argument("--tophatExecutable", "-te",
                          help="""Provide Tophat executable if it's not present in
                        your path""")

    optional.add_argument("--preset", "-p", default='sensitive-local',
                          choices=['very-fast', 'fast', 'sensitive',
                                   'very-sensitive', 'very-fast-local', 'fast-local',
                                   'sensitive-local', 'very-sensitive-local'],
                          help="Provide preset for bowtie2")

    optional.add_argument("--tophatPreset", "-tp",
                          choices=['very-fast', 'fast', 'sensitive',
                                   'very-sensitive'],
                          default='very-sensitive',
                          help="Provide preset for Tophat")

    optional.add_argument("--mismatch", "-m", type=int,
                          choices=[0, 1], default=1,
                          help="""Number of seed mismatches as represented in
                        bowtie2 as -N""")

    optional.add_argument("--reverseComplement", "-rc",
                          action="store_true",
                          help="""Align to reverse complement of reference as
                        represented in bowtie2 as --norc""")

    optional.add_argument("--unaligned", "-un",
                          action="store_true",
                          help="""Whether to keep unaligned reads in the output
                        sam file. Represented in bowtie2 as --no-unal""")

    optional.add_argument("--threads", "-n", default=1,
                          help="Specify the number of threads")

    return parser


def main():

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

    argCheck = Arguments(args)
    argCheck.validateAlign()

    if args.run == 'bowtie2':
        b = Bowtie(args=args)
        b.run(output='smallRNA')
        s = SAM(fileName=b.outputSam)
        filteredFasta = s.filterPotentialChimeras(target=b.outputSam)
        b.run(output='targetRNA', filtered=filteredFasta)

    if args.run == 'tophat':
        t = Tophat(args=args)
        t.run()

if __name__ == "__main__":
    main()
