import argparse
import os
import textwrap

import clashchimeras.log
from clashchimeras.initialize import Releases, Index


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
        description=textwrap.dedent("""\
              Downloads required sequences and create bowtie2 indexes
              required for alignment"""),
        usage='An example usage is: %(prog)s -gor "H.sapiens" -mor hsa',
        add_help=True)

    return parser


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    required.add_argument("--gencodeOrganism", "-gor",
                          default="H.sapiens", choices=["H.sapiens",
                                                        "M.musculus"],
                          help="""Select model organism""",
                          required=True)

    required.add_argument("--mirbaseOrganism", "-mor",
                          default='hsa',
                          help="""Select organism to download microRNAs for""",
                          required=True)

    return parser


def getOutputArgs():
    parser = argparse.ArgumentParser(add_help=False)

    output = parser.add_argument_group('Output')

    output.add_argument("--path", "-pa",
                        help="""Location where all the database files and
                      indexes will be downloaded""",
                        default='~/db/CLASHChimeras',
                        metavar='path')

    return parser


def getOptionalArgs():
    parser = argparse.ArgumentParser(add_help=False)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--logLevel", "-l",
                          help="Set logging level",
                          default='INFO',
                          choices=['INFO', 'DEBUG', 'WARNING', 'ERROR'])

    optional.add_argument("--bowtieExecutable", "-be",
                          help="""Provide bowtie2 executable if it's not present
                        in your path""")

    optional.add_argument("--tophatExecutable", "-te",
                          help="""Provide Tophat executable if it's not present in
                        your path""")

    optional.add_argument("--miRNA", "-mi",
                          choices=['mature', 'hairpin'], default='hairpin',
                          help="Which miRNA sequences to align")

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

    r = Releases(gencodeOrganism=args.gencodeOrganism,
                 mirbaseOrganism=args.mirbaseOrganism,
                 path=os.path.expanduser(args.path))

    i = Index(root=r.gencodePath)
    i.create()

    logger.warn('Please find the indexes created listed below..')
    logger.warn('Use them when you run align-for-chimeras')

    for _i in i.created:
        logger.warn(_i)

    i = Index(root=r.mirbasePath)
    i.create()

    logger.warn('Please find the indexes created listed below..')
    logger.warn('Use them when you run align-for-chimeras')

    for _i in i.created:
        logger.warn(_i)

if __name__ == "__main__":
    main()
