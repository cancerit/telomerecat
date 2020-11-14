"""
Create length estimates given a set of TELBAMS.

Author: jhrf
"""

import textwrap
import time
import os
import re
import parabam

import numpy as np
import random
import pandas as pd

# from itertools import izip
from collections import namedtuple

from telomerecat.csv2length import Csv2Length
from telomerecat.core import TelomerecatInterface

# import args
from . import add_arg
from functools import partial

from . import RANDOM_SEED

class SimpleReadFactory(object):
  def __init__(self, vital_stats=None, trim_reads=0):
    self._SimpleRead = namedtuple(
      "SimpleRead", "seq qual" + " five_prime pattern mima_loci n_loci" + " avg_qual"
    )

    if vital_stats:
      self._read_len = vital_stats["read_len"]
      self._phred_offset = vital_stats["phred_offset"]
    else:
      self._read_len = 100
      self._phred_offset = 33

    self._trim_reads = trim_reads
    self._compliments = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    self.mima_logic = MismatchingLociLogic()

  def get_simple_read(self, read):
    seq, qual = self.__flip_and_compliment__(read)
    if self._trim_reads > 0:
      seq, qual = (seq[: self._trim_reads], qual[: self._trim_reads])

    (mima_loci, frameshift_loci), pattern = self.mima_logic.get_telo_mismatch(seq, qual)

    avg_qual, n_loci = self.__get_phred_score__(seq, qual, mima_loci, frameshift_loci)

    simple_read = self._SimpleRead(
      seq, qual, self.__get_five_prime__(pattern), pattern, mima_loci, n_loci, avg_qual
    )

    # self.mima_logic.print_mima(seq, qual, pattern)

    return simple_read

  def __get_phred_score__(self, seq, qual, mima_loci, frameshift_loci):
    if len(mima_loci) + len(frameshift_loci) == 0:
      return 0, 0

    remove_mimas = []
    phreds = []

    for start, end in frameshift_loci:
      fuse_phreds = []

      for i in range(start, end):
        if i in mima_loci:
          remove_mimas.append(i)
        fuse_phreds.append(ord(qual[i]) - self._phred_offset)

      phreds.append(min(fuse_phreds))

    for loci in mima_loci:
      if loci not in remove_mimas:
        phreds.append(ord(qual[loci]) - self._phred_offset)

    return np.mean(phreds), len(phreds)

  def __get_average_qual__(self, qual, mima_loci):
    if len(mima_loci) == 0:
      return 0
    phreds = np.array([ord(qual[i]) - self._phred_offset for i in mima_loci])
    return np.mean(phreds)

  def __trim_seq__(self, seq, qual):
    cutoff = 0
    min_sequence = 0
    for q in qual:
      qual_byte = ord(q) - self._phred_offset
      if qual_byte == 0:
        min_sequence += 1
      else:
        min_sequence = 0

      if min_sequence == 5:
        cutoff = cutoff - 5
        break

      cutoff += 1

    return seq[:cutoff], qual[:cutoff]

  def __get_five_prime__(self, pattern):
    if pattern is None:
      return None
    else:
      return pattern == "CCCTAA"

  def __get_pattern__(self, seq):
    cta, tag = "CCCTAA", "TTAGGG"
    pattern = None
    if cta in seq or tag in seq:
      if seq.count(cta) > seq.count(tag):
        pattern = cta
      else:
        pattern = tag
    return pattern

  def __flip_and_compliment__(self, read):
    if read.is_reverse:
      compliments = self._compliments
      seq_compliment = [compliments[base] for base in read.seq]
      seq_compliment = "".join(seq_compliment)
      return (seq_compliment[::-1], read.qual[::-1])
    else:
      return (read.seq, read.qual)


class MismatchingLociLogic(object):
  def get_telo_mismatch(self, seq, qual):

    c_rich_count = seq.count("CCCTAA")
    g_rich_count = seq.count("TTAGGG")

    if c_rich_count > g_rich_count:
      return self.get_mismatching_loci(seq, qual, "CCCTAA"), "CCCTAA"
    elif g_rich_count > c_rich_count:
      return self.get_mismatching_loci(seq, qual, "TTAGGG"), "TTAGGG"
    else:
      c_mima, c_fuse = self.get_mismatching_loci(seq, qual, "CCCTAA")
      g_mima, g_fuse = self.get_mismatching_loci(seq, qual, "TTAGGG")

      c_score = len(c_mima) + len(c_fuse)
      g_score = len(g_mima) + len(g_fuse)

      if c_score < g_score:
        return (c_mima, c_fuse), "CCCTAA"
      else:
        return (g_mima, g_fuse), "TTAGGG"

  def get_mismatching_loci(self, seq, qual, pattern):
    segments = re.split("(%s)" % (pattern,), seq)

    segments = self.join_complete_segments(segments, pattern)

    segment_offsets = self.get_segment_offsets(segments)
    self.extend_offsets(segments, pattern, segment_offsets)

    mima_loci, fuse_loci = self.offsets_to_loci(seq, qual, pattern, segment_offsets)

    return mima_loci, fuse_loci

  def join_complete_segments(self, segments, pattern):

    new_segements = []
    current_segment = ""

    for segment in segments:
      if segment == "":
        continue
      elif segment == pattern:
        current_segment += pattern
      else:
        if len(current_segment) > 0:
          new_segements.append(current_segment)
          current_segment = ""
        new_segements.append(segment)

    if len(current_segment) > 0:
      new_segements.append(current_segment)
    return new_segements

  def offsets_to_loci(self, seq, qual, pattern, segment_offsets):

    mima_loci = []
    fuse_loci = []
    seq_len = len(seq)

    # print segment_offsets
    # print [seq[s:e] for s, e in segment_offsets]

    for start, end in segment_offsets:
      segment = seq[start:end]
      segment_qual = qual[start:end]

      if start > 0 and start < end:
        # deletion event
        # +1 for exclusive upperrange
        end_of_range = min((start + 1, seq_len,))
        fuse_loci.append((start - 1, end_of_range))

      if pattern not in segment:
        if start > end:
          # These are fused loci
          fuse_loci.append((end, start))

        elif start < end:
          if (end - start) > 1:
            segment_mima = self.compare_to_telo(segment, segment_qual, pattern)

            segment_mima = [s + start for s in segment_mima]
            mima_loci.extend(segment_mima)
          else:
            mima_loci.append(start)

        elif start == end:
          # complete merge
          pass

    fuse_loci = self.filter_fuse_loci(mima_loci, fuse_loci)
    return mima_loci, fuse_loci

  def filter_fuse_loci(self, mima_loci, fuse_loci):
    remove_candidates = []
    offset = 0
    for loci in mima_loci:
      add_to_offset = 0
      for fuse_id, (start, end) in enumerate(fuse_loci[offset:]):
        if start <= loci and loci < end:
          remove_candidates.append(offset + fuse_id)
        elif loci > end:
          add_to_offset += 1
        elif loci < start:
          offset += add_to_offset
          break

    filtered_fuse = []
    for fuse_id, fuse in enumerate(fuse_loci):
      if fuse_id not in remove_candidates:
        filtered_fuse.append(fuse)

    all_fuse_loci = []
    for start, end in fuse_loci:
      all_fuse_loci.extend(range(start, end))

    merged_fuse = []
    current_fuse = []
    for i in set(all_fuse_loci):
      current_fuse.append(i)
      if len(current_fuse) > 1 and (current_fuse[-1] - current_fuse[-2]) > 1:
        merged_fuse.append((current_fuse[0], current_fuse[-1] + 1,))
        current_fuse = []

    return filtered_fuse

  def compare_to_telo(self, seq, qual, pattern):
    comparisons = []
    best_score = float("inf")

    for i in range(len(pattern)):
      comparison_seq = self.telo_sequence_generator(pattern, i)
      mima_loci = []
      qual_bytes = []
      score = 0
      for s, (l, c) in zip(seq, comparison_seq):
        if s != c:
          mima_loci.append(l)
          qual_bytes.append(ord(qual[l]))

          score += 1
          if score > best_score:
            break

      if score <= best_score:
        best_score = score

        if len(mima_loci) == 0:
          avg_phred = 0
        else:
          avg_phred = np.mean(qual_bytes)
        comparisons.append((score, list(mima_loci), avg_phred,))

    return self.get_best_offset(best_score, comparisons)

  def get_best_offset(self, best_score, comparisons):
    best_comparisons = []
    for score, loci, avg_phred in comparisons:
      if score == best_score:
        best_comparisons.append((avg_phred, loci))
    best_comparisons.sort(key=lambda x: x[0])
    return best_comparisons[0][1]

  def extend_offsets(self, segments, pattern, segment_offsets):
    for seg_id, segment in enumerate(segments):
      if pattern in segment:

        cur_seg_offsets = segment_offsets[seg_id]

        if seg_id > 0:
          # extend_backwards
          prev_seg_offsets = segment_offsets[seg_id - 1]

          new_offset = self.compare_to_pattern(
            segments[seg_id - 1], pattern, reverse=True
          )

          prev_seg_offsets[1] = prev_seg_offsets[1] - new_offset
          cur_seg_offsets[0] = cur_seg_offsets[0] - new_offset

        if seg_id < (len(segments) - 1):
          # extend_forwards
          next_seg_offsets = segment_offsets[seg_id + 1]
          new_offset = self.compare_to_pattern(segments[seg_id + 1], pattern)

          next_seg_offsets[0] = next_seg_offsets[0] + new_offset
          cur_seg_offsets[1] = cur_seg_offsets[1] + new_offset

  def get_segment_offsets(self, segments):
    segment_offsets = []
    offset = 0

    for segment in segments:
      segment_offsets.append(
        [offset, offset + len(segment),]
      )
      offset += len(segment)

    return segment_offsets

  def telo_sequence_generator(self, pattern, offset):

    offset_pattern = (pattern[offset:] + pattern)[: len(pattern)]
    i = 0
    while True:
      yield i, offset_pattern[i % len(pattern)]
      i += 1

  def compare_to_pattern(self, seq, pattern, reverse=False):
    i = 0
    generator = self.get_sequence_generator(seq, pattern, reverse)

    for c, s in generator:
      if c == s:
        i += 1
      else:
        break

    return i

  def get_sequence_generator(self, seq, pattern, reverse):
    def forward_gen(seq, pattern):
      for c, s in zip(seq, pattern):
        yield c, s
      return

    def reverse_gen(seq, pattern):
      for c, s in zip(seq[::-1], pattern[::-1]):
        yield c, s
      return

    if reverse:
      generator = reverse_gen(seq, pattern)
    else:
      generator = forward_gen(seq, pattern)

    return generator

  def print_mima(self, seq, qual, pat):
    print("-")
    loci_status, mima_loci, fuse_loci = self.get_loci_status(seq, qual, pat)

    print("Mima: %s" % str(mima_loci))
    print("Fuse: %s", str(fuse_loci))
    print(seq)
    print(loci_status)
    print(qual)

  def get_loci_status(self, seq, qual, pat):
    mima_loci, fuse_loci = self.get_mismatching_loci(seq, qual, pat)
    loci_status = []

    for i in range(len(seq)):
      if i in mima_loci:
        loci_status.append("X")
      else:
        loci_status.append("_")

    for start, end in fuse_loci:
      loci_status[start:end] = ["F"] * (end - start)

    return "".join(loci_status), mima_loci, fuse_loci


class VitalStatsFinder(object):
  def __init__(self, temp_dir, total_procs, task_size, trim_length=0):
    self.temp_dir = temp_dir
    self._total_procs = total_procs
    self._task_size = task_size
    self._trim_length = trim_length

  def __csv_to_dict__(self, stats_path):
    insert_dat = pd.read_csv(stats_path).to_dict(orient="record")[0]

    ins_N = int(insert_dat["N"])
    if ins_N == 0:
      insert_mean = -1
      insert_sd = -1
    else:
      ins_sum = int(insert_dat["sum"])
      ins_power_2 = int(insert_dat["power_2"])

      insert_mean, insert_sd = self.__get_mean_and_sd__(ins_sum, ins_power_2, ins_N)

    min_qual = int(insert_dat["min_qual"])
    qual_mean, qual_sd = self.__get_mean_and_sd__(
      insert_dat["qual_sum"], insert_dat["qual_power_2"], insert_dat["qual_N"]
    )

    return {
      "insert_mean": insert_mean,
      "insert_sd": insert_sd,
      "min_qual": min_qual,
      "max_qual": int(insert_dat["max_qual"]),
      "read_len": int(insert_dat["read_len"]),
      "qual_mean": qual_mean,
      "qual_sd": qual_sd,
    }

  def __get_mean_and_sd__(self, x_sum, x_power_2, x_N):

    x_mean = x_sum / x_N
    x_sd = np.sqrt(float((x_N * x_power_2)) - float(x_sum ** 2)) / x_N

    return x_mean, x_sd

  def get_vital_stats(self, sample_path):

    vital_stats_csv = self.__run_vital_rule__(sample_path)
    vital_stats = self.__csv_to_dict__(vital_stats_csv)
    vital_stats["phred_offset"] = vital_stats["min_qual"]
    vital_stats["initial_read_len"] = vital_stats["read_len"]

    if self._trim_length > 0:
      vital_stats["read_len"] = self._trim_length

    return vital_stats

  def __run_vital_rule__(self, sample_path, keep_in_temp=True):
    def rule(read, constants, master):
      stats = {}

      if read.is_read1 and read.is_proper_pair and read.mapq > 38:
        insert_size = abs(read.template_length)
        stats["sum"] = insert_size
        stats["power_2"] = insert_size ** 2
        stats["N"] = 1

      stats["read_len"] = len(read.seq)
      byte_vals = [ ord(char) for char in read.qual ]

      min_qual = min(byte_vals)
      max_qual = max(byte_vals)

      qual_mean = np.mean(byte_vals)
      stats["qual_sum"] = qual_mean
      stats["qual_power_2"] = qual_mean ** 2
      stats["qual_N"] = 1

      stats["min_qual"] = min_qual
      stats["max_qual"] = max_qual

      return stats

    structures = {}

    structures["sum"] = {"data": 0, "store_method": "cumu"}
    structures["power_2"] = {"data": 0, "store_method": "cumu"}
    structures["N"] = {"data": 0, "store_method": "cumu"}
    structures["read_len"] = {"data": 0, "store_method": "max"}

    structures["min_qual"] = {"data": 999, "store_method": "min"}
    structures["max_qual"] = {"data": 0, "store_method": "max"}

    structures["qual_sum"] = {"data": 0, "store_method": "cumu"}
    structures["qual_power_2"] = {"data": 0, "store_method": "cumu"}
    structures["qual_N"] = {"data": 0, "store_method": "cumu"}

    stat_interface = parabam.Stat(
      temp_dir=self.temp_dir,
      total_procs=self._total_procs,
      task_size=10000,
      keep_in_temp=keep_in_temp,
    )

    out_paths = stat_interface.run(
      input_paths=[sample_path], constants={}, rule=rule, struc_blueprint=structures
    )

    return out_paths["global"]["stats"]


class ReadStatsFactory(object):
  def __init__(self, temp_dir, total_procs=4, task_size=5000, trim_reads=0, seed_randomness=False, debug_print=False):

    self.temp_dir = temp_dir
    self._total_procs = total_procs
    self._task_size = task_size

    self._debug_print = debug_print
    self._trim = trim_reads
    self.seed_randomness = seed_randomness

  def get_read_counts(self, path, vital_stats, sample_name=None, error_profile=None, error_path=None, error_list=False):
    read_stat_paths = self.run_read_stat_rule(path, vital_stats)

    read_array = self.__path_to_read_array__(read_stat_paths["read_array"])

    # only build error profile if global error isn't passed into this function
    if error_profile is None:
      error_profile = self.__paths_to_error_profile__(read_stat_paths)
      # save error profile if error_path is specified
      if error_path is not None:
        if error_list:
          # use list format if more than one error profile needs to be saved
          temp_error_path = os.path.dirname(error_path) + "/" + str(sample_name)
          temp_error_path = temp_error_path.replace(".bam", "_error_profile.txt")
          self.__save_error_profile__(error_profile, temp_error_path)
          with open(error_path, "a") as f:
            # add error path name to list txt file
            f.write(str(os.path.basename(temp_error_path)) + '\n')
        else:
          self.__save_error_profile__(error_profile, error_path)

    # always build vairance
    sample_variance = self.__paths_to_sample_variance__(read_stat_paths)

    read_counts = self.read_array_to_counts(read_array, error_profile, sample_variance)

    self.__delete_analysis_paths__(read_stat_paths)

    return read_counts

  def get_global_error(self, bulk_path, vital_stats, error_path=None):
    read_stat_paths = self.run_read_stat_rule(bulk_path, vital_stats)
    read_array = self.__path_to_read_array__(read_stat_paths["read_array"])
    error_profile = self.__paths_to_error_profile__(read_stat_paths)

    # save error profile if error_path is specified
    if error_path is not None:
      self.__save_error_profile__(error_profile, error_path)

    self.__delete_analysis_paths__(read_stat_paths)

    return error_profile

  def __delete_analysis_paths__(self, read_stat_paths):
    for analysis, path in read_stat_paths.items():
      os.remove(path)

  def __paths_to_error_profile__(self, read_stat_paths):
    random_counts = pd.read_csv(read_stat_paths["random_counts"], header=None).values
    read_counts = pd.read_csv(read_stat_paths["mima_counts"], header=None).values
    error_profile = self.__get_error_profile__(read_counts, random_counts)
    return error_profile

  def __save_error_profile__(self, error_profile, error_path):
      np.savetxt(error_path, error_profile.astype(int), delimiter=',')

  def __paths_to_sample_variance__(self, read_stat_paths):
    read_counts = pd.read_csv(read_stat_paths["mima_counts"], header=None).values
    sample_variance = self.__get_sample_variance__(read_counts)
    return sample_variance

  def __get_sample_variance__(self, read_counts):
    read_counts[0, :] = 0
    mask = read_counts > 0
    mask[40:, :] = False

    read_variance = read_counts[mask].std() / read_counts[mask].mean()
    return read_variance

  def __get_error_profile__(self, read_counts, random_counts, thresh=None):
    error_profile = self.__get_significantly_enriched__(read_counts, random_counts, thresh)

    error_profile = self.__remove_noise__(error_profile)
    error_profile = self.__prune_error_profile__(error_profile)
    error_profile = self.__rationalise_error_profile__(error_profile)

    error_profile[: int(read_counts.shape[0] * 0.1), :] = 1

    return error_profile

  def __get_significantly_enriched__(self, read_counts, random_counts, thresh=None):
    dif_counts = read_counts - random_counts
    ten_percent = int(read_counts.shape[0] * 0.1)

    if thresh is None:
      mask = self.__get_exclusion_mask__(read_counts)
      arg_max_index = (read_counts * mask).argmax()
      dif_loci_x, dif_loci_y = np.unravel_index(arg_max_index, dif_counts.shape)

      hi_thresh = dif_counts[
        int(dif_loci_x - ten_percent) : int(dif_loci_x + ten_percent),
        dif_loci_y - 15 : dif_loci_y + 1,
      ]
      hi_thresh = hi_thresh.flatten()

      # return a matrix of 0's when hi_thresh is invalid
      if hi_thresh.shape[0] <= 0:
        return np.zeros(dif_counts.shape)

      thresh = np.percentile(hi_thresh, 95)

    if self._debug_print:
      print("Thresh:", thresh)

    error_profile = (dif_counts * (dif_counts > 0)) > thresh
    error_profile = error_profile * 1

    return error_profile

  def __remove_noise__(self, error_profile):
    row_max, col_max = error_profile.shape
    error_profile[int(row_max * 0.11) :, int(col_max * 0.7) :] = 0
    error_profile[int(row_max * 0.40) :, 0] = 0
    error_profile[int(row_max * 0.55) :, :] = 0

    return error_profile

  def __get_exclusion_mask__(self, read_counts):
    mask = np.zeros(read_counts.shape)
    x_start = int(read_counts.shape[0] * 0.2)
    y_start = int(read_counts.shape[1] * 0.5)
    mask[x_start:, y_start:] = 1
    return mask

  def __prune_error_profile__(self, error_profile):

    isolated_mask = self.__get_isolated_mask__(error_profile)
    continuous_mask = self.__get_continuous_mask__(error_profile)

    combined_mask = (isolated_mask + continuous_mask) > 0

    return error_profile * combined_mask

  def __get_continuous_mask__(self, error_profile):
    row_continuous = self.__transform_continous_matrix__(error_profile)
    col_continuous = self.__transform_continous_matrix__(error_profile.transpose())
    col_continuous = col_continuous.transpose()

    max_continuous = row_continuous * (row_continuous > col_continuous)
    max_continuous = max_continuous + (col_continuous * (col_continuous >= row_continuous))

    continusous_mask = max_continuous >= 4
    return continusous_mask

  def __transform_continous_matrix__(self, error_profile):
    continuous = np.zeros(error_profile.shape)

    for row_i in range(0, error_profile.shape[0]):
      sequence = 0
      start_index = -1
      for col_i in range(0, error_profile.shape[1]):
        value = error_profile[row_i, col_i]
        if value == 0:
          if sequence > 0:
            continuous[row_i, start_index:col_i] = sequence
            sequence = 0
            start_index = -1
        elif value > 0:
          if start_index == -1:
            start_index = col_i
          sequence += 1
    return continuous

  def __get_isolated_mask__(self, error_profile):
    isolated_mask = np.ones(error_profile.shape)
    first_locis = self.__get_first_loci__(error_profile)
    for row_i in range(1, error_profile.shape[0]):
      if first_locis[row_i] == -1:
        # skip rows with no entries
        continue
      else:
        for col_i in range(error_profile.shape[1]):
          if (not error_profile[row_i, col_i]) or col_i == 0:
            continue
          elif self.__prune_decision__(row_i, col_i, error_profile):
            isolated_mask[row_i, col_i] = 0
    return isolated_mask

  def __prune_decision__(self, row_i, col_i, error_profile):
    neighbours = [
      (row_i - 1, col_i + 1),
      (row_i - 1, col_i),
      (row_i - 1, col_i - 1),
      (row_i, col_i + 1),
      (row_i, col_i - 1),
      (row_i + 1, col_i + 1),
      (row_i + 1, col_i),
      (row_i + 1, col_i - 1),
    ]

    try:
      return self.__get_neighbor_sum__(row_i, col_i, error_profile, neighbours) < 4
    except IndexError:
      return False

  def __get_neighbor_sum__(self, row_i, col_i, error_profile, neighbours):
    neighbours_sum = sum([error_profile[r, c] for (r, c) in neighbours])
    return neighbours_sum

  def __get_first_loci__(self, error_profile):
    first_loci = []
    for row_i in range(error_profile.shape[0]):
      if any(error_profile[row_i, :]):
        for col_i in range(error_profile.shape[1]):
          if error_profile[row_i, col_i]:
            first_loci.append(col_i)
            break
      else:
        first_loci.append(-1)
    return first_loci

  def __rationalise_error_profile__(self, error_profile):
    if error_profile.sum() > 0:
      # ensure each row makes sense
      start_row = np.where(error_profile)[0].max()
      global_loci = 0
      for i in reversed(range(0, start_row + 1)):
        error_bins_in_row = np.where(error_profile[i, :])[0]
        if len(error_bins_in_row) > 0:
          cur_loci = error_bins_in_row.max()
        else:
          cur_loci = 0

        global_loci = cur_loci
        error_profile[i, : global_loci + 1] = True

      # ensure each column makes sense
      start_col = np.where(error_profile)[1].max()
      global_loci = 0
      for j in reversed(range(0, start_col + 1)):
        error_bins_in_col = np.where(error_profile[:, j])[0]
        if len(error_bins_in_col) > 0:
          cur_loci = error_bins_in_col.max()
        else:
          cur_loci = 0

        global_loci = cur_loci
        error_profile[:global_loci + 1, j] = True

    return error_profile

  def __array_to_file__(self, array, unique):
    df = pd.DataFrame(array)
    out_path = "./%s-tmctout.csv" % (unique)
    df.to_csv(out_path, index=False, header=False)
    return out_path

  def __path_to_read_array__(self, read_array_path):
    return pd.read_csv(read_array_path, header=None).values

  def read_array_to_counts(self, read_array, error_profile, sample_variance):
    complete_reads, boundary_reads = self.__get_complete_status__(read_array, error_profile)

    f2_count, f4_count = self.__get_boundary_counts__(boundary_reads)
    f1_count = self.__get_f1_count__(complete_reads)

    return_dat = {
      "F2": int(f2_count),
      "F1": int(f1_count),
      "F4": f4_count,
      "sample_variance": sample_variance,
    }

    return return_dat

  def __get_f1_count__(self, complete_reads):
    return float(complete_reads.shape[0]) / 2

  def __get_boundary_counts__(self, boundary_reads):
    f2_count, f4_count, total_reads = self.__get_read_counts__(boundary_reads)
    return f2_count, f4_count

  def __get_read_counts__(self, boundary_reads):
    f2_count = sum(boundary_reads[:, 3] == 1)
    f4_count = sum(boundary_reads[:, 3] == 0)
    total_reads = boundary_reads.shape[0]
    return f2_count, f4_count, total_reads

  def __get_complete_status__(self, read_array, error_profile):
    boundary_indicies = []
    complete_indicies = []

    for i in range(int(read_array.shape[0])):
      read_info = list(map(int, read_array[i, [0, -2]]))
      pair_info = list(map(int, read_array[i, [2, -1]]))

      read = error_profile[read_info[0], read_info[1]]
      pair = error_profile[pair_info[0], pair_info[1]]

      if read and pair:
        complete_indicies.append(i)
      elif (not read) and pair:
        boundary_indicies.append(i)

    return read_array[complete_indicies, :], read_array[boundary_indicies, :]

  def run_read_stat_rule(self, path, vital_stats, keep_in_temp=True):

    simple_read_factory = SimpleReadFactory(vital_stats, trim_reads=self._trim)
    phred_offset = vital_stats["phred_offset"]

    maxtrix_max = (vital_stats["max_qual"] - phred_offset) + 1
    matrix_shape = (vital_stats["read_len"] + 1, maxtrix_max)

    def get_return_stats(reads):

      return_stats = [
        len(reads[0].mima_loci),
        int(reads[0].five_prime),
        len(reads[1].mima_loci),
        int(reads[1].five_prime),
        reads[0].avg_qual,
        reads[1].avg_qual,
      ]

      return return_stats

    def rule(reads, constants, master, seed_randomness):
      simple_reads = [simple_read_factory.get_simple_read(read) for read in reads]
      return_dat = np.zeros((2, 6))
      return_dat[0, :] = get_return_stats(simple_reads)
      return_dat[1, :] = get_return_stats(simple_reads[::-1])

      random_counts = np.zeros(matrix_shape)
      mima_counts = np.zeros(matrix_shape)

      for read in simple_reads:
        mima_counts[read.n_loci, int(read.avg_qual)] += 1

        sample_size = len(read.mima_loci)
        if sample_size > 0:
          if seed_randomness:
            random.seed(RANDOM_SEED)
          rand_quals = random.sample(list(read.qual), sample_size)
          qual_bytes = [ord(q) - phred_offset for q in rand_quals]
          rand_avg = np.mean(qual_bytes)

          random_counts[int(sample_size), int(rand_avg)] += 1

      results = {
        "read_array": np.array(return_dat),
        "random_counts": random_counts,
        "mima_counts": mima_counts,
      }

      return results

    structures = {
      "read_array": {"data": np.zeros((2, 6)), "store_method": "vstack"},
      "mima_counts": {"data": np.zeros(matrix_shape), "store_method": "cumu"},
      "random_counts": {"data": np.zeros(matrix_shape), "store_method": "cumu"},
    }

    stat_interface = parabam.Stat(
      temp_dir=self.temp_dir,
      pair_process=True,
      total_procs=self._total_procs,
      task_size=self._task_size,
      keep_in_temp=keep_in_temp,
      verbose=0,
    )

    out_paths = stat_interface.run(
      input_paths=[path], constants={}, rule=partial(rule, seed_randomness=self.seed_randomness), struc_blueprint=structures
    )

    return out_paths[path]


class Telbam2Length(TelomerecatInterface):
  def __init__(
    self,
    temp_dir=None,
    task_size=10000,
    total_procs=8,
    reader_n=2,
    verbose=False,
    announce=True,
    cmd_run=False,
  ):

    super(Telbam2Length, self).__init__(
      instance_name="telomerecat telbam2length",
      temp_dir=temp_dir,
      task_size=task_size,
      total_procs=total_procs,
      reader_n=reader_n,
      verbose=verbose,
      announce=announce,
      cmd_run=cmd_run,
    )

  def run_cmd(self):
    coverages = None
    num_tels = None

    if self.cmd_args.file_input:  # input is txt files containing telbam paths
      if self.cmd_args.cov_ntel:  # file input contains coverage and number of telomeres too
        telbams_paths = []
        coverages = []
        num_tels = []
        for file in self.cmd_args.input:
          with open(file, 'r') as f:
            for line in f:
              items = line.split(",")
              telbams_paths.append(items[0])
              coverages.append(items[1])
              num_tels.append(items[2])
      else:  # txt file input just has telbam paths
        telbams_paths = []
        for file in self.cmd_args.input:
          with open(file, 'r') as f:
            lines = [line.strip() for line in f]
            telbams_paths.extend(lines)

    else:  # input is already a list of telbam paths
      telbams_paths = self.cmd_args.input

    self.run(input_paths=telbams_paths,
            trim=self.cmd_args.trim,
            output_path=self.cmd_args.output,
            simulator_n=self.cmd_args.simulator_runs,
            correct_f2a=self.cmd_args.enable_correction,
            inserts_path=self.cmd_args.insert,
            bulk_path=self.cmd_args.pseudobulk,
            error_path=self.cmd_args.error_path,
            error_list=self.cmd_args.error_list,
            cov_ntel=self.cmd_args.cov_ntel,
            coverages=coverages,
            num_tels=num_tels,
            seed_randomness=self.cmd_args.seed_randomness
    )

  def run(
    self,
    input_paths,
    trim=0,
    output_path=None,
    correct_f2a=False,
    simulator_n=10,
    inserts_path=None,
    bulk_path=None,
    error_path=None,
    error_list=False,
    cov_ntel=False,
    coverages=None,
    num_tels=None,
    seed_randomness=False
  ):

    """The main function for invoking the part of the
       program which creates a telbam from a bam

    Arguments:
      inputs_paths (list): The TELBAMs that we wish to estimate TL estimates for
      output_path (string): Specify a path to output results to (optional)
      inserts_path (string): A path to a file containing insert length estimates
                   for each TELBAM. Formatted as follows:
                    example_telbam.bam, insert_mean, insert_sd
      bulk_path (string): A path to a specified pseudobulk telbam (optional)
      error_path (string): A directory path used to store error profiles as CSVs (optional)
      error_list (bool): Denote whether to store non-global error profiles to a list (optional)
      cov_ntel (bool): Denote whether coverages and num_tels have been provided 
                    within the input_paths txt file (optional)
      coverages (list): Coverages extracted from the original input_paths txt file
      num_tels (list): Number of telomeres per sample extracted from the original input_paths txt file
    """

    self.__introduce__()
    names = [ os.path.basename(path).replace("_telbam", "") for path in input_paths ]

    output_csv_path = output_path
    temp_csv_path = self.__get_temp_path__(cov_ntel)

    insert_length_generator = self.__get_insert_generator__(inserts_path)

    self.__output__(
      " Collecting meta-data for all samples | %s\n" % (self.__get_date_time__(),), 1
    )

    vital_stats_finder = VitalStatsFinder(
      self.temp_dir, self.total_procs, self.task_size, trim
    )

    # find global error profile when pseudobulk telbam path is provided
    if bulk_path is not None:
      vital_stats = vital_stats_finder.get_vital_stats(bulk_path)

      self.__check_vital_stats_insert_size__(inserts_path,
                                             insert_length_generator,
                                             vital_stats)

      # specify bulk sample name in error path
      if error_path is not None:
        current_error_path = error_path
      else:
        current_error_path = None

      global_error_profile = self.__get_global_error__(bulk_path,
                                                       vital_stats,
                                                       self.total_procs,
                                                       trim,
                                                       error_path=current_error_path)
    # set global_error_profile to None so its get calculated on a per-telbam basis in __get_read_types__()
    else:
      global_error_profile = None


    # write read_type_counts to temp csv for each telbam
    for i, sample_path in enumerate(input_paths):
      sample_intro = "\t- %s | %s\n" % (names[i], self.__get_date_time__())

      self.__output__(sample_intro, 2)

      vital_stats = vital_stats_finder.get_vital_stats(sample_path)

      self.__check_vital_stats_insert_size__(
        inserts_path, insert_length_generator, vital_stats
      )

      # create file name for current error profile
      if error_path is not None and global_error_profile is None:
        current_error_path = error_path
      else:
        current_error_path = None

      read_type_counts = self.__get_read_types__(
        sample_path, vital_stats, self.total_procs, trim, \
        seed_randomness, names[i], error_profile=global_error_profile, \
        error_path=current_error_path, error_list=error_list
      )

      # only include coverages and num_tels if they are non-empty lists
      if cov_ntel:
        self.__write_to_csv__(read_type_counts,
                              vital_stats,
                              temp_csv_path,
                              names[i],
                              float(coverages[i]),
                              float(num_tels[i]))
      else:
        self.__write_to_csv__(read_type_counts,
                              vital_stats,
                              temp_csv_path,
                              names[i])

    self.__output__("\n", 1)
    length_interface = Csv2Length(
      temp_dir=self.temp_dir,
      total_procs=self.total_procs,
      verbose=self.verbose,
      announce=False,
      cmd_run=False,
    )

    length_interface.run(
      input_paths=[temp_csv_path],
      output_paths=[output_csv_path],
      correct_f2a=correct_f2a,
      seed_randomness=seed_randomness,
      simulator_n=simulator_n
    )

    self.__print_output_information__(output_csv_path)
    self.__goodbye__()

    return output_csv_path

  def __print_output_information__(self, output_csv_path):
    self.__output__((" Length estimation results " "written to the following file:\n"), 1)
    self.__output__("\t./%s\n\n" % (os.path.basename(output_csv_path,)))

  def __get_insert_generator__(self, inserts_path):
    if inserts_path:
      with open(inserts_path, "r") as inserts_file:
        for line in inserts_file:
          yield map(float, line.split(","))

  def __check_vital_stats_insert_size__(
    self, inserts_path, insert_length_generator, vital_stats
  ):
    if inserts_path:
      insert_mean, insert_sd = insert_length_generator.__next__()
      vital_stats["insert_mean"] = insert_mean
      vital_stats["insert_sd"] = insert_sd
      self.__output__(
        "\t\t+ Using user defined insert size: %d,%d\n" % (insert_mean, insert_sd), 2
      )
    elif vital_stats["insert_mean"] == -1:
      default_mean, default_sd = 350, 25
      vital_stats["insert_mean"] = default_mean
      vital_stats["insert_sd"] = default_sd
      self.__output__(
        "\t\t+ Failed to estimate insert size. Using default: %d,%d\n"
        % (default_mean, default_sd),
        2,
      )

  def __get_read_types__(
    self, sample_path, vital_stats, total_procs, trim, seed_randomness, sample_name,
    read_stats_factory=None, error_profile=None, error_path=None, error_list=False
  ):

    if read_stats_factory is None:
      read_stats_factory = ReadStatsFactory(
        temp_dir=self.temp_dir, total_procs=total_procs, trim_reads=trim,
        seed_randomness=seed_randomness, debug_print=False
      )

    read_type_counts = read_stats_factory.get_read_counts(sample_path,
                                                          vital_stats,
                                                          sample_name=sample_name,
                                                          error_profile=error_profile,
                                                          error_path=error_path,
                                                          error_list=error_list)
    return read_type_counts

  def __get_global_error__(self, sample_path,
                                 vital_stats,
                                 total_procs,
                                 trim,
                                 read_stats_factory=None,
                                 error_path=None):
    if read_stats_factory is None:
      read_stats_factory = ReadStatsFactory(temp_dir=self.temp_dir,
                                            total_procs=total_procs,
                                            trim_reads=trim,
                                            debug_print=False)

    # build error profile
    error_profile = read_stats_factory.get_global_error(sample_path,
                                                        vital_stats,
                                                        error_path=error_path)

    return error_profile

  def __get_temp_path__(self, cov_ntel=False):
    temp_path = os.path.join(self.temp_dir, "telomerecat_temp_%d.csv" % (time.time()))
    self.__create_output_file__(temp_path, cov_ntel)
    return temp_path

  def __create_output_file__(self, output_csv_path, cov_ntel=False):
    with open(output_csv_path, "w") as total:
      # use if-statement to know whether coverage & num_tel belong in csv header
      if cov_ntel:
        header = ("Sample,F1,F2,F4,Psi,Insert_mean,Insert_sd,"
                  "Read_length,Initial_read_length,coverage,num_tel\n")
      else:
        header = ("Sample,F1,F2,F4,Psi,Insert_mean,Insert_sd,"
                  "Read_length,Initial_read_length\n")
      total.write(header)
    return output_csv_path

  def __write_to_csv__(self, read_type_counts, vital_stats, output_csv_path, name, coverage=None, num_tel=None):
    with open(output_csv_path, "a") as counts:
      # use if-statement to see whether we have a coverage or num_tel to write on this line
      if coverage is None or num_tel is None:
        counts.write("%s,%d,%d,%d,%.3f,%.3f,%.3f,%d,%d\n" % \
                                    (name,
                                     read_type_counts["F1"],
                                     read_type_counts["F2"],
                                     read_type_counts["F4"],
                                     read_type_counts["sample_variance"],
                                     vital_stats["insert_mean"],
                                     vital_stats["insert_sd"],
                                     vital_stats["read_len"],
                                     vital_stats["initial_read_len"]))
      else:
        counts.write("%s,%d,%d,%d,%.3f,%.3f,%.3f,%d,%d,%.6f,%.3f\n" % \
                                    (name,
                                     read_type_counts["F1"],
                                     read_type_counts["F2"],
                                     read_type_counts["F4"],
                                     read_type_counts["sample_variance"],
                                     vital_stats["insert_mean"],
                                     vital_stats["insert_sd"],
                                     vital_stats["read_len"],
                                     vital_stats["initial_read_len"],
                                     coverage,
                                     num_tel))

  def get_parser(self):
    parser = self.default_parser()
    parser.description = textwrap.dedent(
      """\
    %s
    %s

      The telbam2length command allows the user to genereate a telomere
      length estimate from a previously generated TELBAM file.

      Example useage:

      telomerecat telbam2length /path/to/some_telbam.bam

      This will generate a .csv file with an telomere length estimate
      for the `some_telbam.bam` file.

    %s
    """
      % (self.instance_name, self.header_line, self.header_line,)
    )

    arg_list = ['input_telbam',
                'output_csv',
                'trim',
                'nreads_for_task',
                'file_input',
                'pseudobulk',
                'error_path',
                'error_list',
                'cov_ntel']

    for arg_name in arg_list:
      add_arg[arg_name](parser)

    return parser


if __name__ == "__main__":
  print("Do not run this script directly. Type `telomerecat` for help.")
