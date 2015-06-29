import argparse
import textwrap

from colors import yellow
from progress.bar import Bar

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
              #
              # epilog=textwrap.dedent('''\
              #
              #
              #   Additional help for setting annotation files
              #   ============================================
              #
              #    XX represents the Gencode Release. The Gencode database files
              #    are present in your --path folder when you ran
              #    download-for-chimeras
              #
              #     --smallRNAAnnotation -sa
              #       * If your biotype is tRNA, set the gtf file as
              #         gencode.vXX.tRNAs.gtf found in Gencode folder
              #       * If your biotype is miRNA, set the gtf file as hsa.gff3
              #         found in Mirbase folder
              #       * If your biotype is snRNA, snoRNA, misc_RNA, set the gtf
              #         as gencode.vXX.chr_patch_hapl_scaff.annotation.gtf found
              #         in Gencode folder
              #
              #     --targetRNAAnnotation -ta
              #       * If your biotype is protein_coding, snRNA, snoRNA,
              #         misc_RNA, set the gtf file as
              #         gencode.vXX.chr_patch_hapl_scaff.annotation.gtf found
              #         in Gencode folder
              #       * If your biotype is tRNA, set the gtf file as
              #         gencode.vXX.tRNAs.gtf found in Gencode folder
              # '''),

              add_help=False)

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

  optional.add_argument("--help", "-h", action="help",
                        help="Show this help message and exit")

  optional.add_argument("--smallRNAAnnotation", "-sa",
                        action="store",
                        help="""Provide smallRNA annotation gtf(from Gencode)
                        or gff3(from Mirbase). Only provide gtf from Gencode or
                        gff3 from Mirbase. Does not support other gtf files""")

  optional.add_argument("--targetRNAAnnotation", "-ta",
                        action="store",
                        help="""Provide targetRNA annotation gtf(from Gencode).
                        Only provide gtf from Gencode. Does not support other
                        gtf files""")

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

  argCheck = Arguments(args)
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

  updateInterval = int(len(common)/100)
  pbar = Bar(width=80, max=100, fill=yellow(u"█"),
    suffix="%(percent)d%% - %(eta)ds")

  for index, queryName in enumerate(common, start=1):
    smallRNARecord = smallRNA.access(queryName)
    targetRNARecord = targetRNA.access(queryName)

    chimeraString = clashchimeras.methods.chimeraOrNot(
      smallRNARecord.cigarString, targetRNARecord.cigarString,
      overlap=args.overlap, gap=args.gap)

    if chimeraString:

      if args.getGenomicLocationsSmallRNA:

        smallRNALocation = smallRNAAnnotation.coordinates(
                              smallRNARecord.refId,
                              start=smallRNARecord.start,
                              end=smallRNARecord.end)
        output.writeSmallRNABed(queryName, smallRNALocation,
                                targetRNARecord.refId)

      if args.getGenomicLocationsTargetRNA:

        targetRNALocation = targetRNAAnnotation.coordinates(
                              targetRNARecord.refId,
                              start=targetRNARecord.start,
                              end=targetRNARecord.end)

        output.writeTargetRNABed(queryName, targetRNALocation,
                                 smallRNARecord.refId)

      output.write(queryName, smallRNARecord, targetRNARecord)


    if index % updateInterval == 0:
      pbar.next()
  pbar.finish()
  logger.info('Determining chimeras finished')

if __name__ == "__main__":
   main()
