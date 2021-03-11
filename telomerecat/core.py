"""
This module defines a class that all other commands inherit from.

It's primary responsibility is to define the default parameters for argparse

Author: jhrf
"""

from argparse import SUPPRESS
import parabam
from abc import ABCMeta


class TelomerecatInterface(parabam.core.Interface, metaclass=ABCMeta):

  def __init__(
    self,
    instance_name,
    temp_dir=None,
    task_size=10000,
    total_procs=4,
    reader_n=2,
    verbose=False,
    announce=True,
    cmd_run=False,
  ):

    super(TelomerecatInterface, self).__init__(
      instance_name=instance_name,
      temp_dir=temp_dir,
      task_size=task_size,
      total_procs=total_procs,
      reader_n=reader_n,
      verbose=verbose,
      announce=announce,
      cmd_run=cmd_run,
    )

  def __output__(self, outstr, level=-1):
    if self.verbose and (self.verbose >= level or level == -1):
      print(outstr)

  def default_parser(self):
    parser = super(TelomerecatInterface, self).default_parser()

    parser.add_argument(
      "--insert",
      metavar="CSV",
      type=str,
      default=None,
      help=(
        "A file specifying the insert length mean and\n"
        "std for each input sample. If not present\n"
        "telomerecat will automatically estimate\n"
        "insert length of sample [Default: None]"
      ),
    )
    parser.add_argument(
      "-N",
      "--simulator_runs",
      metavar="INT",
      type=int,
      default=10,
      help=(
        "The amount of times to run the length simulator.\n"
        "A higher number better captures the uncertainty \n"
        "produced by the insert length\n"
        "distribution [Default 10]"
      ),
    )
    parser.add_argument(
      "-e",
      "--enable_correction",
      action="store_true",
      default=False,
      help="Correction will be applied to F2a values",
    )

    parser.add_argument(
      "--seed_randomness",
      action="store_true",
      default=False,
      # Use the seed (42) for all randomnesses in telomerecat to produce stable results.
      help=SUPPRESS  # this option is hidden from users.
    )
    parser.add_argument(
        '-w', '--prior_weight',
        metavar='INT', type=int, default=3,
        help=('The weight given to the prior expectation\n'
              'in F2a correction [Default 3]'))

    return parser
