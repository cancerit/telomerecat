"""
Genereate length estimates for given BAM files.

This simply strings together bam2telbam
and telbam2length packages

Author: jhrf
"""

import textwrap
from argparse import SUPPRESS
from telomerecat.core import TelomerecatInterface

# import args
from . import add_arg

class Bam2Length(TelomerecatInterface):
  def __init__(
    self,
    temp_dir=None,
    task_size=250000,
    total_procs=4,
    reader_n=1,
    verbose=False,
    announce=True,
    cmd_run=False,
  ):

    super(Bam2Length, self).__init__(
      instance_name="telomerecat bam2length",
      temp_dir=temp_dir,
      task_size=task_size,
      total_procs=total_procs,
      reader_n=reader_n,
      verbose=verbose,
      announce=announce,
      cmd_run=cmd_run,
    )

  def run_cmd(self):
    # TODO: existence of self.cmd_args.outbam_dir should be checked before proceeding
    self.run(
      input_paths=self.cmd_args.input,
      output_path=self.cmd_args.output,
      inserts_path=self.cmd_args.insert,
      outbam_dir=self.cmd_args.outbam_dir,
      correct_f2a=self.cmd_args.enable_correction,
      simulator_n=self.cmd_args.simulator_runs,
      seed_randomness=self.cmd_args.seed_randomness
    )

  def run(
    self,
    input_paths,
    output_path=None,
    outbam_dir=None,
    inserts_path=None,
    correct_f2a=False,
    simulator_n=10,
    seed_randomness=False
  ):

    # Import here to avoid infinite loop on import
    from telomerecat.bam2telbam import Bam2Telbam
    from telomerecat.telbam2length import Telbam2Length

    self.__introduce__()

    telbam_interface = Bam2Telbam(
      temp_dir=self.temp_dir,
      total_procs=self.total_procs,
      task_size=self.task_size,
      reader_n=self.reader_n,
      verbose=self.verbose,
      announce=False,
    )

    out_files = telbam_interface.run(input_paths=input_paths, outbam_dir=outbam_dir)

    length_paths = self.collapse_out_files(out_files)
    length_interface = Telbam2Length(
      temp_dir=self.temp_dir,
      total_procs=self.total_procs,
      task_size=self.task_size,
      reader_n=self.reader_n,
      announce=False,
      verbose=self.verbose,
    )

    length_interface.run(
      input_paths=length_paths,
      output_path=output_path,
      inserts_path=inserts_path,
      simulator_n=simulator_n,
      correct_f2a=correct_f2a,
      seed_randomness=seed_randomness
    )

    self.__goodbye__()
    self.interface_exit()

  def collapse_out_files(self, out_files):
    length_paths = []
    for path in out_files.keys():
      length_paths.append(out_files[path]["telbam"])
    return length_paths

  def get_parser(self):
    parser = self.default_parser()
    parser.description = textwrap.dedent(
      """\
    %s
    %s

      The bam2length command allows the user to genereate a telomere
      length estimate from a BAM file.

      The majority of the time taken running the bam2length script is
      spent collecting telomeric reads from the BAM file. By default
      this command will retain the TELBAM generted as part
      of the analysis.

      If you wish to generate TELBAMS seperately from length estimation
      you should use the bam2telbam command.

      Type `telomerecat bam2telbam` to find out more.

    %s
    """
      % (self.instance_name, self.header_line, self.header_line,)
    )

    for arg_name in ['input_bam', 'outbam_dir', 'output_csv', 'trim', 'nreads_for_task', 'reference']:
      add_arg[arg_name](parser)

    parser.add_argument("-s", help=SUPPRESS)
    parser.add_argument("-f", help=SUPPRESS)

    return parser


if __name__ == "__main__":
  print("Please do not run this script directly. Type telomerecat -h for more information.")


# ....happily ever after.
