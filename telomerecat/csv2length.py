"""
Generate TL estimates given a CSV file.

Run the telomerecat length estimate on a CSV containing the
the following fields:
  F1, F2, F4, Psi, Insert_mean, Insert_sd

Author: jhrf
"""

import sys
import textwrap
import time
import random
import math

import numpy as np
import pandas as pd

from argparse import SUPPRESS
from functools import partial
from shutil import copy
from multiprocessing import Pool, freeze_support

from telomerecat import core

from . import RANDOM_SEED

class LengthSimulator(object):
  def __init__(self, insert_mu, insert_sigma, complete, boundary, read_len, seed_randomness=False):

    self._insert_mu = insert_mu
    self._insert_sigma = insert_sigma

    self._obs_total = complete + boundary  # observed total
    self._read_thresh = int(self._obs_total * 0.0001)

    self._complete = complete
    self._boundary = boundary
    self._read_len = read_len
    self._read_sim = self.__simulate_reads__

    self._get_read_len = lambda: self._read_len
    self.seed_randomness = seed_randomness

  def start(self):
    if self._boundary <= 0 or self._complete <= 0:
      return 0
    else:
      return self.__run_sim__()

  def __run_sim__(self):
    found = False

    tel_len = int((float(self._complete) / self._boundary) * (self._insert_mu - 100))
    max_its = 5000
    its = 0

    best_result = float("inf")
    best_tel_len = tel_len

    read_simulator = self._read_sim
    get_factor = self.__get_factor__

    while not found and its < max_its:
      sim_comp, sim_boun, invalid_reads = read_simulator(tel_len)
      #   result can be positive or negative thus
      #   influencing elongation and shortening
      #   A negative diference (overestimate) results
      #   in a shorter telomere length next iteration
      result = self.is_same_ratio(sim_comp, sim_boun)
      abs_result = abs(result)

      if abs_result < best_result:
        best_result = abs_result
        best_tel_len = tel_len

      if result == 0:
        found = True
      else:
        mult = -1 if result < 0 else 1
        factor = get_factor(abs_result)
        tel_len += factor * mult
      its += 1

    return best_tel_len

  def __get_factor__(self, abs_result):
    total = self._obs_total
    if self.seed_randomness:
      random.seed(RANDOM_SEED)
    if abs_result < (total * 0.001):
      factor = 1
    elif abs_result < (total * 0.2):
      factor = random.sample([1, 2, 2, 4, 25, 200, 500], 1)[0]
    elif abs_result < (total * 0.5):
      factor = random.sample([20, 50, 100, 200, 1000, 2000], 1)[0]
    else:
      factor = int(np.log(abs_result) * 300) / random.sample([1, 2], 1)[0]
    return factor

  def is_same_ratio(self, comp, boun):
    difference = abs(comp - self._complete)
    if difference <= self._read_thresh:
      return 0
    else:
      return self._complete - comp

  def __simulate_reads__(self, tel_len):
    if tel_len < self._insert_mu:
      return int(self._complete * 0.50), self._boundary, 0

    #  speedup
    is_complete = self.is_complete
    is_boundary = self.is_boundary

    insert_mean = self._insert_mu
    insert_sigma = self._insert_sigma

    obs_total = self._obs_total
    #  speedup

    est_complete = 0
    est_boundary = 0
    invalid_read = 0
    total = 0

    while total < (obs_total):
      if self.seed_randomness:
        random.seed(RANDOM_SEED)
      insert_size = random.gauss(insert_mean, insert_sigma)
      if self.seed_randomness:
        random.seed(RANDOM_SEED)
      location = random.randint(0, int(math.floor(tel_len)))
      if is_complete(location, insert_size):
        est_complete += 1
        total += 1
      elif is_boundary(location, insert_size):
        est_boundary += 1
        total += 1
      else:
        invalid_read += 1

    return (est_complete, est_boundary, invalid_read)

  def is_complete(self, location, insert_size):
    return (location - insert_size) > 0

  def is_boundary(self, location, insert_size):
    return (location - self._get_read_len()) > 0


class Csv2Length(core.TelomerecatInterface):
  def __init__(self, temp_dir=None, total_procs=8, verbose=False, announce=True, cmd_run=False):

    super(Csv2Length, self).__init__(
      instance_name="telomerecat csv2length",
      temp_dir=temp_dir,
      total_procs=total_procs,
      verbose=verbose,
      announce=announce,
      cmd_run=cmd_run,
    )

  def run_cmd(self):
    output_paths = self.__handle_cmd_outpaths__()

    self.run(
      input_paths=self.cmd_args.input,
      correct_f2a=not self.cmd_args.disable_correction,
      output_paths=output_paths,
      prior_weight=self.cmd_args.prior_weight,
      seed_randomness=self.cmd_args.seed_randomness,
      simulator_n=self.cmd_args.simulator_runs
    )

  def run(self, input_paths, output_paths=[], correct_f2a=True, prior_weight=3, seed_randomness=False, simulator_n=10):

    self.__introduce__()
    self.__generate_output_paths__(input_paths, output_paths)

    self.__output__(" Commencing length estimation | %s\n" % (self.__get_date_time__(),), 1)

    for input_path, output_path in zip(input_paths, output_paths):
      if self.announce:
        self.__output__("\tInput: %s\n" % (input_path,), 1)
        self.__output__(" \tOutput: %s\n" % (output_path,), 1)

      counts = pd.read_csv(input_path)
      counts = self.__get_length_from_dataframe__(
        counts, simulator_n, correct_f2a, prior_weight, seed_randomness
      )

      self.__output_length_results__(counts, output_path)
      self.__output__("\n", 1)

    self.__goodbye__()
    return output_paths

  def __handle_cmd_outpaths__(self):
    output_paths = None
    if self.cmd_args.output is None:
      output_paths = []
    else:
      output_paths = self.cmd_args.output.split(",")
      output_paths = [p.strip() for p in output_paths]

    return output_paths

  def __get_length_from_dataframe__(
    self, counts, simulator_n, correct_f2a=True, prior_weight=3, seed_randomness=False, simulate_lengths=True,
  ):

    counts["F2a"] = counts["F2"] - counts["F4"]

    if correct_f2a:
      counts["F2a_c"] = self.__get_corrected_f2a__(counts, prior_weight)
    else:
      counts["F2a_c"] = counts["F2a"]

    if "coverage" in counts.columns and "num_tel" in counts.columns:
      # use coverage & number of telomeres to estimate mean telomere length
      lengths = self.__get_cov_ntel_lengths__(counts)
      counts["Length"] = lengths
      # use simulator to generate length standard deviation (uncertainty based on coverage)
      __, counts["Length_std"] = self.__get_lengths__(counts, seed_randomness, simulator_n)
    elif simulate_lengths:
      counts["Length"], counts["Length_std"] = self.__get_lengths__(counts, seed_randomness, simulator_n)
    else:
      lengths = self.__quick_length__(counts)
      counts["Length"] = lengths
      counts["Length_std"] = [0.000] * len(lengths)  # stdev is always 0 with quick_length
    return counts

  def __get_corrected_f2a__(self, counts, prior_weight=3):
    # include very small float eps so that 0/0 = NaN turns into 0/eps = 0
    theta_observed = counts["F2a"] / (counts["F2"] + counts["F4"] + np.finfo(float).eps)

    prior_weight = 3
    theta_expected = sum(theta_observed * counts["F2"]) / (sum(counts["F2"]) + np.finfo(float).eps)

    theta_corrected = (
      (theta_observed * (counts["Psi"])) + (theta_expected * prior_weight)
    ) / ((counts["Psi"]) + (prior_weight))

    corrected_f2_counts = (counts["F2"] + counts["F4"]) * theta_corrected

    return corrected_f2_counts.round(3)

  def __get_cov_ntel_lengths__(self, counts):
    """ Calculate telomere length based on coverage, number of telomere, and read counts. """
    # TODO: estimate length uncertainty for each sample (row) based on 
    # missingness of coverage, return these uncertainties as length_stds
    lengths = []
    for i, sample in counts.iterrows():
      length = ((sample["F1"] * 2 + sample["F2a_c"]) * sample["Read_length"]) / \
          (sample["coverage"] * sample["num_tel"])
      lengths.append(round(length, 3))
    return lengths

  def __quick_length__(self, counts):
    lengths = []
    for i, sample in counts.iterrows():
      factor = sample["Read_length"] / (sample["F2a_c"] - sample["Read_length"])
      scale = sample["F2a_c"] + (sample["F2a_c"] * factor)
      length = ((sample["F1"] / scale) * sample["Insert_mean"]) + sample["Insert_mean"]
      lengths.append(round(length, 3))
    return lengths

  def __get_lengths__(self, counts, seed_randomness, simulator_n):
    length_means = []
    length_stds = []
    for i, sample in counts.iterrows():
      sample_intro = "\t- %s | %s\n" % (sample["Sample"], self.__get_date_time__())
      self.__output__(sample_intro, 2)

      # just say length is NA for cases that cause errors
      if sample["Insert_sd"] <= 0 or np.isnan(sample["F2a_c"]):
        length_mean = "NA"
        length_std = "NA"
      else:
        length_mean, length_std = run_simulator_par(
              sample["Insert_mean"],
              sample["Insert_sd"],
              sample["F1"],
              sample["F2a_c"],
              self.total_procs,
              sample["Read_length"],
              seed_randomness,
              simulator_n)

      length_means.append(length_mean)
      length_stds.append(length_std)
    return length_means, length_stds

  def __generate_output_paths__(self, input_paths, output_paths):
    if len(output_paths) < len(input_paths):
      for i in range(len(input_paths)):
        path = "./telomerecat_length_%d_%d.csv" % (time.time(), i)
        output_paths.append(path)

  def __output_length_results__(self, counts, count_path):
    counts.to_csv(count_path, index=None)

  def __copy_out_of_temp__(self, file_paths, copy_path="."):
    for fil in file_paths:
      copy(fil, copy_path)

  def get_parser(self):
    parser = self.default_parser()
    parser.description = textwrap.dedent(
      """\
    %s
    %s

      The csv2length command allows the user to genereate new TL
      estimates from a previously generated CSV. This could be useful
      if the user wishes to re-run analysis with different options 
      (such as F2a correction) and does not wish to reanalyse
      the TELBAM files.

      Example useage:

      telomerecat csv2length /path/to/telomerecat_length_123.csv

      This will generate a new .csv file with revised telomere length
      estimates for the `telomerecat_length_123.csv` file.

    %s
    """
      % (self.instance_name, self.header_line, self.header_line,)
    )

    parser.add_argument(
      "input", metavar="CSV(s)", nargs="+", help="The CSV(s) that we wish to reanalyse"
    )
    parser.add_argument(
      "--output",
      metavar="CSV",
      type=str,
      nargs="?",
      default=None,
      help=(
        "Specify output path for length estimation CSV.\n"
        "If more than one input CSV is provided\n"
        "the user will need to provide a comma seperated list\n"
        "with one entry per input. For example:\n\t"
        '--output="zzz.csv, yyy.csv". [Default: None]'
      ),
    )
    parser.add_argument("-s", help=SUPPRESS)
    parser.add_argument("-f", help=SUPPRESS)

    return parser


"""
 These functions are implemented as global
 beacuse multiprocessing does not work with
 class bound methods
"""


def check_results(sim_results):
  if 0 in sim_results:
    sys.stderr.write(
      "[WARNING] Telomere length reported zero. This means telomercat\n"
      + "\tfailed to identify enough complete or boundary reads.\n"
      + "\tThis  may mean your original sample was preprocessed to remove \n"
      + "\ttelomere reads. Alternatively this sample could have \n"
      + "\tvery short average TL.\n"
    )


def run_simulator(insert_mu, insert_sigma, complete, boundary, proc, read_len, n=10):

  simmer = LengthSimulator(insert_mu, insert_sigma, complete, boundary, read_len)
  res = []
  invalid_count = []
  for i in range(n):
    length, invalid = simmer.start()
    res.append(length)
    invalid_count.append(invalid)

  check_results(res)
  return (np.mean(res), np.std(res))


def estimator_process(job, insert_mu, insert_sigma, complete, boundary, read_len, seed_randomness):

  length_estimator = LengthSimulator(insert_mu, insert_sigma, complete, boundary, read_len, seed_randomness)
  results = length_estimator.start()
  return results


def run_simulator_par(insert_mu, insert_sigma, complete, boundary, proc, read_len, seed_randomness, simulator_n):

  freeze_support()
  p = Pool(proc)
  sim_partial = partial(
    estimator_process,
    insert_mu=insert_mu,
    insert_sigma=insert_sigma,
    complete=complete,
    boundary=boundary,
    read_len=read_len,
    seed_randomness=seed_randomness
  )

  results = p.map(sim_partial, range(simulator_n))
  p.close()
  check_results(results)
  return (np.mean(results), np.std(results))


if __name__ == "__main__":
  print ("Do not run this script directly. Type `telomerecat` for help.")
