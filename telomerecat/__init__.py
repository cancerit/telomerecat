import sys
from .version import __version__

# oh 42, the answer to everything.
RANDOM_SEED = 42

def add_input_telbam_arg(parser):
  parser.add_argument(
    "input", metavar="TELBAM(S)", nargs="+", help="The TELBAM(s) that we wish to analyse"
  )


def add_input_bam_arg(parser):
  parser.add_argument(
    "input",
    metavar="BAM(S)",
    nargs="+",
    help=("BAM file(s) for which we wish to\n" "generate telomere length estimates"),
  )


def add_trim_arg(parser):
  parser.add_argument(
    "-t",
    "--trim",
    type=int,
    metavar="INT",
    default=0,
    help="Use only the amount of sequence specified by this  \n"
    "option (i.e if the value 90 is supplied\n"
    "then only the first 90 bases are\n"
    "considered) [Default: Whole read]",
  )


def add_outbam_dir_arg(parser):
  parser.add_argument(
    "--outbam_dir",
    metavar="DIR",
    type=str,
    default=None,
    help=(
      "Output folder for bams which have all reads with at least 2 occurences\n"
      "of the telomeric hexamer. [Default: None, telbams will be discarded]"
    ),
  )


def add_output_csv_arg(parser):
  parser.add_argument(
    "--output",
    metavar="PATH",
    type=str,
    default="./telomerecat_length.csv",
    help=(
      "Specify output path for length estimation CSV.\n"
      "[Default: './telomerecat_length.csv']"
    ),
  )


def add_nreads_for_each_task_arg(parser):
  parser.add_argument(
    "-s",
    type=int,
    metavar="INT",
    default=250000,
    help=("The amount of reads considered by each\n" "distributed task. [Default: 250000]"),
  )


def add_file_input_arg(parser):
  parser.add_argument(
    '-i', '--file_input',
    action="store_true",
    default=False,
    help="Specify whether the input file is a telbam or a txt file\n"
          "that contains one telbam file per row. If the -cnt flag is\n"
          "also used then the input txt file also contains coverage\n"
          "and number of chromosomes on each row (separated by comma)."
  )


def add_pseudobulk_arg(parser):
  parser.add_argument(
    '-b', '--pseudobulk',
    metavar='TELBAM',
    type=str,
    nargs='?',
    default=None,
    help="Path to pseudobulk telbam that gets used to create a bulk error\n"
          "profile and sample variance that is used to categorize\n"
          "read types for all cells. A telomere length estimate won't\n"
          "be given for this pseudobulk telbam. [Default: None]"
  )


def add_error_path_arg(parser):
  parser.add_argument(
    '-ep', '--error_path',
    metavar='PATH',
    type=str,
    nargs='?',
    default=None,
    help="Specify a path to save error_profile arrays to.\n"
          "[Default: None]"
  )


def add_error_list_arg(parser):
  parser.add_argument(
    '-el', '--error_list',
    action="store_true",
    default=False,
    help="Specify whether the non-global error profiles should be\n"
          "recorded in a list. [Default: None]"
  )


def add_cov_ntel_arg(parser):
  parser.add_argument(
    '-cnt', '--cov_ntel',
    action="store_true",
    default=False,
    help="Use this option when input txt file also has coverage and number\n"
          "of telomeres. Telomere length will be calculated using F1, F2a_c, coverage,\n"
          "and number of telomeres in this case."
  )

add_arg = {
  'input_telbam': add_input_telbam_arg,
  'input_bam': add_input_bam_arg,
  'trim': add_trim_arg,
  'outbam_dir': add_outbam_dir_arg,
  'output_csv': add_output_csv_arg,
  'nreads_for_task': add_nreads_for_each_task_arg,
  'file_input': add_file_input_arg,
  'pseudobulk': add_pseudobulk_arg,
  'error_path': add_error_path_arg,
  'error_list': add_error_list_arg,
  'cov_ntel': add_cov_ntel_arg
}

def exit_with_msg(message):
  sys.stderr.write(message)
  sys.stderr.flush()
  sys.exit(1)
