"""
Create TELBAMS given a set of BAM files.

TELBAMS are comprised of all the reads in a BAM file that include the
canonic TTAGGG telomere sequence

Author: jhrf
"""

import os
import textwrap
import gc

import parabam

from . import add_arg, exit_with_msg

# Here we define what constitutes a telomereic read.
# If both reads in a pair match the critea we store
# them to the telbam
def rule(reads, constants, master):
  tel_pats = constants["tel_pats"]
  telomere_status = False

  for read in iter(reads):
    for pattern in tel_pats:
      if read.seq is not None and pattern in read.seq:
        telomere_status = True
        break

  results = []
  if telomere_status:
    results = [("telbam", reads[0]), ("telbam", reads[1])]
  return results


# This whole package is essentially just a wrapper for a call to parabam subset
class Bam2Telbam(parabam.core.Interface):
  """Interface to interact with telomerecat programatically. This interface
  is also called by the `telomerecat` script"""

  def __init__(
    self,
    temp_dir=None,
    task_size=250000,
    total_procs=4,
    reader_n=2,
    verbose=False,
    announce=True,
    cmd_run=False,
  ):

    super(Bam2Telbam, self).__init__(
      instance_name="telomerecat bam2telbam",
      temp_dir=temp_dir,
      task_size=task_size,
      total_procs=total_procs,
      reader_n=reader_n,
      verbose=verbose,
      announce=announce,
      cmd_run=cmd_run,
    )

  def run_cmd(self):
    """Called from the master `telomerecat` script which handels the
    cmd line interface. Requires an argparse parser. Users should call
    the `run` function."""

    self.run(input_paths=self.cmd_args.input, outbam_dir=self.cmd_args.outbam_dir)

  def run(self, input_paths, outbam_dir=None):
    """The main function for invoking the part of the
       program which creates a telbam from a bam

    Arguments:
      bams (list): The BAM files we wish to run telomerecat telbam on
      total_procs (int): The maximum numbers of task that will be run at one time
      task_size (int): The amount of reads that any one task will process concurrently
      verbose (int): Expects an int from 0 to 2.
               The level of output produced by telomerecat
      keep_in_temp (bool): Files will be kept in temp file after processing.
                 Useful for incorporation into pipelines
      announce (bool): Specify whether the program returns a welcome string."""

    self.__introduce__()

    subset_types = ["telbam"]
    tel_pats = ["TTAGGGTTAGGG", "CCCTAACCCTAA"]

    # need to define my constants and engine here:
    telbam_constants = {"thresh": 1, "tel_pats": tel_pats}

    final_output_paths = {}

    keep_in_temp = outbam_dir is None
    if not keep_in_temp:
      # check if the folder is writable
      if not os.path.exists(outbam_dir):
        exit_with_msg(f"Error: can not find outbam_dir path: '{outbam_dir}'")
      if not os.access(outbam_dir, os.W_OK | os.X_OK):
        exit_with_msg(f"Error: do not have right permission to write into outbam_dir path: '{outbam_dir}'")

    for input_path in input_paths:

      if self.verbose:
        print(f" Generating TELBAM from: {input_path}")
        print(f"\t- TELBAM generation started {self.__get_date_time__()}")

      subset_interface = parabam.Subset(
        temp_dir=self.temp_dir,
        total_procs=self.total_procs,
        task_size=self.task_size,
        reader_n=self.reader_n,
        verbose=self.verbose,
        pair_process=True,
        include_duplicates=True,
        keep_in_temp=keep_in_temp,
      )
      # call to parabam subset
      telbam_paths = subset_interface.run(
        input_paths=[input_path],
        subsets=subset_types,
        constants=telbam_constants,
        rule=rule,
        outbam_dir=outbam_dir
      )
      if self.verbose:
        print(f"\t- TELBAM generation finished {self.__get_date_time__()}")

      gc.collect()
      final_output_paths.update(telbam_paths)

    self.__goodbye__()
    return final_output_paths

  def get_parser(self):
    parser = self.default_parser()
    parser.description = textwrap.dedent(
      """\
    %s
    %s

       The bam2telbam command allows you to generate a TELBAM
       from a parent BAM file. A TELBAM is a file including all of the
       reads with at least 2 occurences of the telomeric hexamer.

       Once you have genereated a TELBAM you may then generate length
       estimates more quickly, when compared to running the `bam2length`
       command on a full BAM file. This is helpful if you intend to
       generate TL estimates more than once or if you require a
       collection of telomere reads.

       Given the following BAM file:

         example_bam_name.bam

       `telomerecat bam2telbam` will create the following TELBAM in the
       directory which it is run:

         example_bam_name_telbam.bam

       To find out how to generate a length estimate from a TELBAM
       type `telomerecat telbam2length` into your terminal

    %s
     """
      % (self.instance_name, self.header_line, self.header_line,)
    )

    for arg_name in ['input_bam', 'outbam_dir']:
      add_arg[arg_name](parser)

    return parser


if __name__ == "__main__":
  print("Type telomerecat -h for help!")

# ....happily ever after.
