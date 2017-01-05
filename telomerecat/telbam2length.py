import sys
import textwrap
import time
import os
import re
import parabam
import pdb

import numpy as np
import pandas as pd

from shutil import copy
from itertools import izip
from collections import namedtuple,Counter
from pprint import pprint

from telomerecat import Csv2Length
from telomerecat.core import TelomerecatInterface

######################################################################
##
##      Create a length estimate given a set of TELBAMS 
##
##      Author: jhrf
##
######################################################################

class SimpleReadFactory(object):

    def __init__(self, vital_stats=None, trim_reads=0):
        self._SimpleRead = namedtuple("SimpleRead","seq qual" +
                                      " five_prime pattern mima_loci n_loci"+
                                      " avg_qual")

        if vital_stats:
            self._read_len = vital_stats["read_len"]
            self._phred_offset = vital_stats["phred_offset"]
        else:
            self._read_len = 100
            self._phred_offset = 33

        self._trim_reads = trim_reads
        self._compliments = {"A":"T","T":"A",
                             "C":"G","G":"C",
                             "N":"N"}

        self.mima_logic = MismatchingLociLogic()

    def get_simple_read(self,read):
        seq,qual = self.__flip_and_compliment__(read)
        if self._trim_reads > 0:
            seq,qual = (seq[:self._trim_reads],
                        qual[:self._trim_reads])

        (mima_loci,frameshift_loci),pattern = \
                self.mima_logic.get_telo_mismatch(seq)

        avg_qual,n_loci = self.__get_phred_score__(qual, 
                                                   mima_loci, 
                                                   frameshift_loci)

        simple_read = self._SimpleRead(
            seq,
            qual,
            self.__get_five_prime__(pattern),
            pattern,
            mima_loci,
            n_loci,
            avg_qual)

        return simple_read

    def __get_phred_score__(self, qual, mima_loci, frameshift_loci):
        if len(mima_loci) + len(frameshift_loci) == 0:
            return 0, 0

        remove_mimas = []
        phreds = []

        for start, end in frameshift_loci:
            fuse_phreds = []

            for i in xrange(start,end):
                if i in mima_loci:
                    remove_mimas.append(i)
                fuse_phreds.append(ord(qual[i])-self._phred_offset)

            phreds.append(min(fuse_phreds))

        for loci in mima_loci:
            if loci not in remove_mimas:
                phreds.append(ord(qual[loci])-self._phred_offset)

        return np.mean(phreds), len(phreds)

    def __get_average_qual__(self,qual,mima_loci):
        if len(mima_loci) == 0:
            return 0
        phreds = np.array([ord(qual[i])-self._phred_offset for i in mima_loci])
        return np.mean(phreds)

    def __trim_seq__(self,seq,qual):
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

        return seq[:cutoff],qual[:cutoff]

    def __get_five_prime__(self,pattern):
        if pattern is None:
            return None
        else:
            return pattern == "CCCTAA"

    def __get_pattern__(self,seq):
        cta,tag = "CCCTAA","TTAGGG"
        pattern = None
        if cta in seq or tag in seq:   
            if seq.count(cta) > seq.count(tag):
                pattern = cta
            else:
                pattern = tag
        return pattern

    def __flip_and_compliment__(self,read):
        if read.is_reverse:
            compliments = self._compliments
            seq_compliment = map(lambda base: compliments[base],read.seq)
            seq_compliment = "".join(seq_compliment)
            return(seq_compliment[::-1],read.qual[::-1])
        else:
            return (read.seq,read.qual)

class MismatchingLociLogic(object):

    def get_telo_mismatch(self, seq):

        c_rich_count = seq.count("CCCTAA")
        g_rich_count = seq.count("TTAGGG")

        if c_rich_count > g_rich_count:
            return self.get_mismatching_loci(seq, "CCCTAA"), "CCCTAA"
        elif g_rich_count > c_rich_count:
            return self.get_mismatching_loci(seq, "TTAGGG"), "TTAGGG"
        else:
            c_mima, c_fuse = self.get_mismatching_loci(seq, "CCCTAA")
            g_mima, g_fuse = self.get_mismatching_loci(seq, "TTAGGG")

            c_score = len(c_mima) + len(c_fuse)
            g_score = len(g_mima) + len(g_fuse)

            if c_score < g_score:
                return (c_mima, c_fuse), "CCCTAA"
            else:
                return (g_mima, g_fuse), "TTAGGG"

    def get_mismatching_loci(self, seq, pattern):
        segments = re.split("(%s)"%(pattern,), seq)
        segments = self.filter_segments(segments)
        segments = self.collapse_segments(segments, pattern)
        segments = self.merge_segments(segments, pattern)
        mima_loci, fuse_loci = self.find_mima_loci(segments, pattern)

        return mima_loci, fuse_loci

    def find_mima_loci(self, segments, pattern):
        mima_loci = []
        fuse_loci = []

        offset = 0

        for i,seg in enumerate(segments):
            neighbours = self.get_neighbour_segments(segments,i)
            
            score, segment_loci = \
                self.get_best_offset_score(seg, pattern, neighbours)
            
            mima_loci.extend(np.array(segment_loci)+offset)

            if pattern in neighbours[0] or pattern in seg:
                if offset-1 not in mima_loci and \
                    offset+1 not in mima_loci and \
                     offset not in mima_loci:

                    if neighbours[0] != "":
                        fuse_loci.append( (offset-1,offset+2,))
            offset += len(seg)

        return mima_loci,fuse_loci

    def get_neighbour_segments(self, segments, i):
        prev_segment = ""
        next_segment = ""

        if i > 0:
            prev_segment = segments[i-1]
        if i < (len(segments)-1):
            next_segment = segments[i+1]

        return prev_segment, next_segment

    def get_best_offset_score(self, seq, 
                                    pattern, 
                                    neighbours, 
                                    override = False):

        best_score = float("Inf")
        best_mima = []

        if pattern in seq and not override:
            best_score, best_mima = 0,[]
        elif len(seq) == 1:
            best_score, best_mima =  1,[0]
        elif len(seq) < len(pattern)*.66:
            for position,neighbour in enumerate(neighbours):
                if len(neighbour) > 0:
                    if position == 0:
                        temp_seq = neighbour+seq
                    else:
                        temp_seq = seq+neighbour
                    score, mima = self.get_best_offset_score(temp_seq, 
                                                             pattern,
                                                             neighbours,
                                                             override=True)
                    if score < best_score:
                        best_score = score
                        if position == 0:
                            best_mima = \
                                (np.array(mima) - len(neighbour)).tolist()
                        else:
                            best_mima = mima

        else:
            for i in xrange(len(pattern)):
                comparison_seq = self.telo_sequence_generator(pattern, i)
                mima_loci = []
                score = 0
                for s, (l, c) in izip(seq, comparison_seq):
                    if s != c:
                        mima_loci.append(l)
                        score += 1
                        if score > best_score:
                            break

                if score < best_score:
                    best_score = score
                    best_mima = mima_loci

        return best_score, best_mima
            
    def telo_sequence_generator(self, pattern, offset):

        offset_pattern = (pattern[offset:] + pattern)[:len(pattern)]
        i = 0
        while True:
            yield i, offset_pattern[i % len(pattern)]
            i += 1

    def filter_segments(self, segments):
        return [seg for seg in segments if seg != '']

    def merge_segments(self, segments, pattern):
        merge_segments = list(segments)

        for i,seg in enumerate(segments):
            if pattern not in seg:
                top_trim = 0
                bottom_trim = len(seg)
                
                #See if valid sequence continues into the corrupted segment
                if i > 0:
                    #Check preccedding segment
                    top_trim = self.compare_to_pattern(seg, pattern)
                if i < len(segments)-1:
                    #Check succeeding segment
                    bottom_trim = self.compare_to_pattern(seg,
                                                          pattern,
                                                          reverse = True)

                #ensure that loci are not merged into both
                # i-1 and i+1 simealteanously
                if bottom_trim < top_trim:
                    bottom_trim = top_trim

                #do merge operation
                if top_trim > 0:
                    merge_segments[i-1] += seg[:top_trim]
                if bottom_trim < len(seg):
                    merge_segments[i+1] = seg[bottom_trim:] +\
                                             merge_segments[i+1]
            
                #correct current segment
                merge_segments[i] = seg[top_trim:bottom_trim]

        return self.filter_segments(merge_segments)

    def compare_to_pattern(self, seq, pattern, reverse = False):
        i = 0
        generator = self.get_sequence_generator(seq, pattern, reverse)

        for c,s in generator:
            if c == s:
                i += 1
            else:
                break

        if reverse:
            return_index =  len(seq) - i
        else:
            return_index = i

        return return_index

    def get_sequence_generator(self, seq, pattern, reverse): 
        def forward_gen(seq, pattern):
            for c,s in zip(seq,pattern):
                yield c, s
            return

        def reverse_gen(seq, pattern):
            for c,s in izip(seq[::-1], pattern[::-1]):
                yield c, s
            return

        if reverse:
            generator = reverse_gen(seq, pattern) 
        else: 
            generator = forward_gen(seq, pattern) 

        return generator

    def collapse_segments(self, segments, pattern):
        collapsed_segments = []
        current_segment = ''
        for seg in segments:
            if seg == pattern:
                current_segment += seg
            else:
                if current_segment != '':
                    collapsed_segments.append(current_segment)
                collapsed_segments.append(seg)
                current_segment = ''

        if current_segment != '':
            collapsed_segments.append(current_segment)

        return collapsed_segments

    def print_mima(self, seq, qual, pat):
        print "-"
        loci_status, mima_loci, fuse_loci =\
            self.get_loci_status(seq, pat)

        print "Mima:",mima_loci
        print "Fuse:",fuse_loci
        print seq
        print loci_status
        print qual

    def get_loci_status(self, seq, pat):
        mima_loci,fuse_loci = self.get_mismatching_loci(seq,pat)
        loci_status = []

        for i in xrange(len(seq)):
            if i in mima_loci:
                loci_status.append("X")
            else:
                loci_status.append("_")

        for start,end in fuse_loci:
            loci_status[start:end] = ["F"] * (end - start)

        return "".join(loci_status), mima_loci, fuse_loci

class VitalStatsFinder(object):

    def __init__(self,temp_dir,total_procs,task_size):
        self.temp_dir = temp_dir
        self._total_procs = total_procs
        self._task_size = task_size

    def __csv_to_dict__(self,stats_path):
        insert_dat = np.genfromtxt(stats_path,delimiter=",",
                                names=True,dtype=("S256",float,float,
                                                         float,float,
                                                         float,float,
                                                         float,float,
                                                         float))

        ins_N = int(insert_dat['N'])
        if ins_N == 0:
            insert_mean = -1
            insert_sd = -1
        else:
            ins_sum = int(insert_dat['sum'])
            ins_power_2 = int(insert_dat['power_2'])

            insert_mean,insert_sd = \
                        self.__get_mean_and_sd__(ins_sum, ins_power_2, ins_N)
    
        min_qual = int(insert_dat['min_qual'])
        qual_mean,qual_sd = self.__get_mean_and_sd__(insert_dat["qual_sum"],
                                                    insert_dat["qual_power_2"],
                                                    insert_dat["qual_N"])

        return {"insert_mean":insert_mean, 
                "insert_sd": insert_sd,
                "min_qual":min_qual,
                "max_qual":int(insert_dat['max_qual']),
                "read_len":int(insert_dat['read_len']),
                "qual_mean":qual_mean,
                "qual_sd":qual_sd}

    def __get_mean_and_sd__(self,x_sum,x_power_2,x_N):
        x_mean = x_sum / x_N
        x_sd = np.sqrt( (x_N * x_power_2) - x_sum**2) / x_N

        return x_mean,x_sd

    def get_vital_stats(self,sample_path):

        vital_stats_csv = self.__run_vital_rule__(sample_path)
        vital_stats = self.__csv_to_dict__(vital_stats_csv)
        vital_stats["phred_offset"] = vital_stats["min_qual"]
 
        return vital_stats

    def __run_vital_rule__(self,sample_path,keep_in_temp=True):
        def rule(read,constants,master):
            stats = {}

            hash_count = read.qual.count("#")
            if read.is_read1 and read.is_proper_pair and read.mapq > 38:
                insert_size = abs(read.template_length)
                stats["sum"] = insert_size
                stats["power_2"] = insert_size**2
                stats["N"] = 1
            
            stats["read_len"] =  len(read.seq)
            byte_vals = map(ord,read.qual)
            min_qual = min(byte_vals)
            max_qual = max(byte_vals)

            qual_mean = np.mean(byte_vals)
            stats["qual_sum"] = qual_mean
            stats["qual_power_2"] = qual_mean**2
            stats["qual_N"] = 1

            stats["min_qual"] = min_qual
            stats["max_qual"] = max_qual

            return stats

        structures = {}

        structures["sum"] = {"data":0,"store_method":"cumu"}
        structures["power_2"] = {"data":0,"store_method":"cumu"}
        structures["N"] = {"data":0,"store_method":"cumu"}
        structures["read_len"] = {"data":0,"store_method":"max"}

        structures["min_qual"] = {"data":999,"store_method":"min"}
        structures["max_qual"] = {"data":0,"store_method":"max"}


        structures["qual_sum"] = {"data":0,"store_method":"cumu"}
        structures["qual_power_2"] = {"data":0,"store_method":"cumu"}    
        structures["qual_N"] = {"data":0,"store_method":"cumu"}  

        stat_interface = parabam.Stat(temp_dir=self.temp_dir,
                                      total_procs = self._total_procs,
                                      task_size = 10000,
                                      keep_in_temp=keep_in_temp)

        out_paths = stat_interface.run(input_paths= [sample_path],
                                       constants = {},
                                       rule = rule,
                                       struc_blueprint = structures)

        return out_paths["global"]["stats"]

class ReadStatsFactory(object):

    def __init__(self,temp_dir,
                      total_procs=4,
                      task_size = 5000,
                      trim_reads = 0,
                      debug_print=False):

        self.temp_dir = temp_dir
        self._total_procs = total_procs
        self._task_size = task_size

        self._debug_print = debug_print
        self._trim = trim_reads

    def get_read_counts(self,path,vital_stats):
        read_stat_paths = self.run_read_stat_rule(path, vital_stats)

        read_array = self.__path_to_read_array__(read_stat_paths["read_array"])


        random_counts = pd.read_csv(read_stat_paths["random_counts"],
                                    header=None).values 
        mima_counts = pd.read_csv(read_stat_paths["mima_counts"],
                                  header=None).values 

        read_counts = self.read_array_to_counts(read_array,
                                                mima_counts,
                                                random_counts,
                                                0)

        self.__delete_analysis_paths__(read_stat_paths)
        
        return read_counts

    def __delete_analysis_paths__(self, read_stat_paths):
        for analysis, path in read_stat_paths.items():
            os.remove(path)

    def __paths_to_error_profile__(self,read_stat_paths):
        random_counts = pd.read_csv(read_stat_paths["random_counts"],
                                    header=None).values 
        read_counts = pd.read_csv(read_stat_paths["mima_counts"],
                                  header=None).values 

        sample_variance = 3
        
        error_profile = self.__array_to_profile__(read_counts, 
                                                  random_counts, 
                                                  0)
        return error_profile, sample_variance

    def __get_sample_variance__(self, read_counts):
        read_counts[0,:] = 0
        mask = read_counts > 0
        mask[40:,:] = False

        read_variance = (read_counts[mask].std() / read_counts[mask].mean())
        return read_variance

    def __array_to_profile__(self, read_counts, random_counts, thresh = None):

        dif_counts = read_counts - random_counts
        dif_counts = dif_counts * (1 * (dif_counts > 0))

        ten_percent = int(dif_counts.shape[0] * .1)

        mask = np.ones(dif_counts.shape)
        mask[:,98:] = 0
        mask[75:,:] = 0
        
        dif_counts = dif_counts * mask

        if thresh is None:
            threshold = self.__get_threshold__(dif_counts)
        else:
            threshold = thresh

        error_profile = (dif_counts > 0) * 1
        error_profile[:ten_percent,:] = 1

        error_profile = self.__prune_error_profile__(error_profile)

        return error_profile

    def __get_threshold__(self, dif_counts, percent):
        relevant = dif_counts[dif_counts > 0]
        threshold = np.percentile(relevant, percent)
        print threshold, percent
        return threshold

    def __remove_noise__(self, error_profile):
        row_max, col_max = error_profile.shape
        error_profile[int(row_max*.11):,int(col_max*.7):] = 0
        error_profile[int(row_max*.35):,0] = 0
        error_profile[int(row_max*.55):,:] = 0 

        return error_profile

    def __get_exclusion_mask__(self, read_counts):
        mask = np.zeros(read_counts.shape)
        x_start = int(read_counts.shape[0] * .2)
        y_start = int(read_counts.shape[1] * .5)
        mask[x_start:,y_start:] = 1
        return mask

    def __prune_error_profile__(self, error_profile):
        first_locis = self.__get_first_loci__(error_profile)
        prune_mask = np.ones(error_profile.shape)
        for row_i in xrange(1,error_profile.shape[0]):
            if first_locis[row_i] == -1:
                continue 
            else:
                for col_i in xrange(error_profile.shape[1]):
                    if (not error_profile[row_i,col_i]) or \
                        col_i == 0:
                            continue
                    elif self.__prune_decision__(row_i,col_i,error_profile):
                        prune_mask[row_i,col_i] = 0
        return error_profile * prune_mask

    def __prune_decision__(self, row_i, col_i, error_profile):
        try:
            return self.__get_neighbor_sum__(row_i,col_i,error_profile) < 4
        except IndexError:
            return False

    def __get_neighbor_sum__(self,row_i, col_i, error_profile):

        neighbours = [(row_i-1,col_i+1),
                      (row_i-1,col_i),
                      (row_i-1,col_i-1),
                      (row_i,col_i+1),
                      (row_i,col_i-1),
                      (row_i+1,col_i+1),
                      (row_i+1,col_i),
                      (row_i+1,col_i-1),]

        neighbours_sum = sum([ error_profile[r,c] for (r,c) in neighbours])
        return neighbours_sum

    def __get_first_loci__(self,error_profile):
        first_loci = []
        for row_i in xrange(error_profile.shape[0]):
            if any(error_profile[row_i,:]):
                for col_i in xrange(error_profile.shape[1]):
                    if error_profile[row_i,col_i]:
                        first_loci.append(col_i)
                        break
            else:
                first_loci.append(-1)
        return first_loci

    def __rationalise_error_profile__(self, error_profile):
        start_row = np.where(error_profile)[0].max()
        global_loci = 0
        for i in reversed(xrange(0,start_row+1)):
            error_bins_in_row = np.where(error_profile[i,:])[0]
            if len(error_bins_in_row) > 0:
                cur_loci = error_bins_in_row.max()
            else:
                cur_loci = 0

            #if cur_loci > global_loci:
            global_loci = cur_loci
            error_profile[i,:global_loci+1] = True
        
        return error_profile

    def __array_to_file__(self,array,unique):
        df = pd.DataFrame(array)
        out_path = "./%s-tmctout.csv" % (unique)
        df.to_csv(out_path,index=False,header=False)
        return out_path

    def __path_to_read_array__(self,read_array_path):
        return pd.read_csv(read_array_path,header=None).values

    def read_array_to_counts(self, read_array, 
                                   mima_counts, 
                                   random_counts, 
                                   sample_variance):

        error_profile = self.__array_to_profile__(mima_counts, 
                                                  random_counts, 
                                                  0)
        counts = self.__counts_for_error_profile__(read_array, 
                                                   error_profile)
        print counts

        counts["sample_variance"] = 3
        return counts

    def __counts_for_error_profile__(self, read_array, error_profile):
        complete_reads,boundary_reads = \
                        self.__get_complete_status__(read_array,error_profile)

        f2_count,f4_count = self.__get_boundary_counts__(boundary_reads)
        f1_count = self.__get_f1_count__(complete_reads)

        return_dat = { "F2":int(f2_count),
                       "F1":int(f1_count),
                       "F4":f4_count,}
        return return_dat

    def __get_f1_count__(self,complete_reads):
        return float(complete_reads.shape[0]) / 2

    def __get_boundary_counts__(self,boundary_reads):
        f2_count,f4_count,total_reads = \
                 self.__get_read_counts__(boundary_reads)

        observed_f2a_ratio = (f2_count-f4_count) / float(total_reads)
        observed_f2a_weight = f2_count
            
        return f2_count, f4_count

    def __get_read_counts__(self,boundary_reads):
        f2_count = sum(boundary_reads[:,3] == 1)
        f4_count = sum(boundary_reads[:,3] == 0)
        total_reads = boundary_reads.shape[0]
        return f2_count,f4_count,total_reads
        
    def __get_complete_status__(self,read_array,error_profile):
        boundary_indicies = []
        complete_indicies = []

        for i in xrange(int(read_array.shape[0])):
            read_info = map(int,read_array[i,[0,-2]])
            pair_info = map(int,read_array[i,[2,-1]])

            read = error_profile[read_info[0],read_info[1]] 
            pair = error_profile[pair_info[0],pair_info[1]]

            if read and pair:
                complete_indicies.append(i)
            elif (not read) and pair:
                boundary_indicies.append(i)

        return read_array[complete_indicies,:],\
                read_array[boundary_indicies,:]

    def run_read_stat_rule(self,path,
                                vital_stats,
                                keep_in_temp=True):

        simple_read_factory = SimpleReadFactory(vital_stats,
                                                trim_reads=self._trim)
        phred_offset = vital_stats["phred_offset"]

        maxtrix_max = (vital_stats["max_qual"] - phred_offset)+1
        matrix_shape = (vital_stats["read_len"]+1,101)

        def get_return_stats(reads,coords):

            return_stats = [len(reads[0].mima_loci),
                            int(reads[0].five_prime),
                            len(reads[1].mima_loci),
                            int(reads[1].five_prime),
                            reads[0].avg_qual,
                            reads[1].avg_qual,
                            coords[0],
                            coords[1]]

            return return_stats

        def qual_coord(loci,quals,offset):
            if len(loci) == 0:
                coord = 0
            elif len(loci) == len(quals):
                coord = 100
            elif len(loci) > 0:
                loci_quals = []
                all_quals = []

                for i in xrange(len(quals)):
                    qual_byte = ord(quals[i]) - offset
                    all_quals.append(qual_byte)
                    if i in loci:
                        loci_quals.append(qual_byte)

                mean_qual = np.mean(all_quals)
                mean_loci = np.mean(loci_quals)

                qual_ratio = mean_loci / mean_qual

                if qual_ratio > 1:
                    coord = 100
                else:
                    #print mean_loci, mean_qual, mean_loci / mean_qual, len(loci)
                    coord = int(round((mean_loci / mean_qual) * 100,0))

            return coord


        def rule(reads, constants, master):
            simple_reads = [simple_read_factory.get_simple_read(read) \
                                                         for read in reads]
            
            random_counts = np.zeros(matrix_shape)
            mima_counts = np.zeros(matrix_shape)

            coords = []

            for read in simple_reads:
                mima_coord = qual_coord(read.mima_loci, read.qual, phred_offset)
                coords.append(mima_coord)
                mima_counts[read.n_loci,mima_coord] += 1

                #sample_size = int(np.random.uniform(1,80))
                sample_size = len(read.mima_loci)
                if sample_size > 0:
                    rand_quals  = \
                        np.random.choice(xrange(len(read.qual)),sample_size)
                    rand_coord = \
                            qual_coord(rand_quals, read.qual, phred_offset)

                    random_counts[int(sample_size),rand_coord] += 1


            return_dat = np.zeros((2,8))
            return_dat[0,:] = get_return_stats(simple_reads,coords)
            return_dat[1,:] = get_return_stats(simple_reads[::-1],coords[::-1])

            results = {"read_array":np.array(return_dat),
                      "random_counts":random_counts,
                      "mima_counts":mima_counts}

            return results

        structures = {"read_array":{"data":np.zeros((2,8)),
                                    "store_method":"vstack"},
                     "mima_counts":{"data":np.zeros(matrix_shape),
                                    "store_method":"cumu"},
                     "random_counts":{"data":np.zeros(matrix_shape),
                                    "store_method":"cumu"},}

        stat_interface = parabam.Stat(temp_dir=self.temp_dir,
                                      pair_process=True,
                                      total_procs = self._total_procs,
                                      task_size = self._task_size,
                                      keep_in_temp = keep_in_temp,
                                      verbose=0)

        out_paths = stat_interface.run(
            input_paths = [path],
            constants = {},
            rule = rule,
            struc_blueprint = structures)

        return out_paths[path]

class Telbam2Length(TelomerecatInterface):

    def __init__(self,
                 temp_dir=None,
                 task_size=10000,
                 total_procs=8,
                 reader_n=2,
                 verbose=False,
                 announce=True,
                 cmd_run=False):

        super(Telbam2Length,self).__init__(instance_name = \
                                            "telomerecat telbam2length", 
                                        temp_dir=temp_dir,
                                        task_size=task_size,
                                        total_procs=total_procs,
                                        reader_n=reader_n,
                                        verbose=verbose,
                                        announce = announce,
                                        cmd_run=cmd_run)

    def run_cmd(self):
        self.run(input_paths = self.cmd_args.input,
                 trim = self.cmd_args.trim,
                 inserts_path = self.cmd_args.insert,
                 correct_f2a = not self.cmd_args.disable_correction,
                 output_path = self.cmd_args.output)

    def run(self,input_paths,
                 trim = 0,
                 output_path = None,
                 correct_f2a = True,
                 inserts_path = None):
        
        """The main function for invoking the part of the 
           program which creates a telbam from a bam

        Arguments:
            inputs_paths (list): The TELBAMs that we wish to estimate TL estimates for
            output_path (string): Specify a path to output results to (optional)
            inserts_path (string): A path to a file containing insert length estimates
                                   for each TELBAM. Formatted as follows:
                                        example_telbam.bam, insert_mean, insert_sd
        """

        self.__introduce__()

        names = map(lambda b: os.path.basename(b),input_paths)
        names = map(lambda nm: nm.replace("_telbam",""),names)

        output_csv_path = self.__get_output_path__(output_path)
        temp_csv_path = self.__get_temp_path__()

        insert_length_generator = self.__get_insert_generator__(inserts_path)
        
        self.__output__(" Collecting meta-data for all samples | %s\n" \
                            % (self.__get_date_time__(),),1)

        vital_stats_finder = VitalStatsFinder(self.temp_dir, 
                                        self.total_procs,
                                        self.task_size)

        for sample_path,sample_name, in izip(input_paths,names):
            sample_intro = "\t- %s | %s\n" % (sample_name, 
                                              self.__get_date_time__())

            self.__output__(sample_intro,2)

            vital_stats = vital_stats_finder.get_vital_stats(sample_path)
            if trim > 0:
                vital_stats["read_len"] = trim

            self.__check_vital_stats_insert_size__(inserts_path,
                                                    insert_length_generator,
                                                    vital_stats)

            read_type_counts = self.__get_read_types__(sample_path,
                                                       vital_stats,
                                                       self.total_procs,
                                                       trim)

            self.__write_to_csv__(read_type_counts,
                                    vital_stats,
                                    temp_csv_path,
                                    sample_name)
        
        self.__output__("\n",1)
        length_interface = Csv2Length(temp_dir=self.temp_dir,
                                       total_procs=self.total_procs,
                                       verbose=self.verbose,
                                       announce=False,
                                       cmd_run=False)

        length_interface.run(input_paths=[temp_csv_path], 
                             output_paths=[output_csv_path],
                             correct_f2a=correct_f2a)

        self.__print_output_information__(output_csv_path)
        self.__goodbye__()

        return output_csv_path

    def __print_output_information__(self, output_csv_path):
        self.__output__((" Length estimation results "
                         "written to the following file:\n"),1)
        self.__output__("\t./%s\n\n" % (os.path.basename(output_csv_path,)))

    def __get_insert_generator__(self,inserts_path):
        if inserts_path:
            with open(inserts_path,"r") as inserts_file:
                for line in inserts_file:
                    yield map(float,line.split(","))

    def __check_vital_stats_insert_size__(self,inserts_path,
                                        insert_length_generator,vital_stats):
        if inserts_path:
            insert_mean,insert_sd = insert_length_generator.next()
            vital_stats["insert_mean"] = insert_mean
            vital_stats["insert_sd"] = insert_sd
            self.__output__("\t\t+ Using user defined insert size: %d,%d\n" \
                                                    % (insert_mean,insert_sd),2)
        elif vital_stats["insert_mean"] == -1:
            default_mean,default_sd = 350,25
            vital_stats["insert_mean"] = 350
            vital_stats["insert_sd"] = 25
            self.__output__("\t\t+ Failed to estimate insert size. Using default: %d,%d\n"\
                                                % (default_mean,default_sd),2)

    def __get_read_types__(self,sample_path,
                                vital_stats,
                                total_procs,
                                trim,
                                read_stats_factory=None):

        if read_stats_factory is None:
            read_stats_factory = ReadStatsFactory(temp_dir=self.temp_dir,
                                                  total_procs=total_procs,
                                                  trim_reads=trim,
                                                  debug_print=False)
            
        read_type_counts = read_stats_factory.get_read_counts(sample_path,
                                                              vital_stats)
        return read_type_counts 

    def __get_temp_path__(self):
        temp_path = os.path.join(self.temp_dir,"telomerecat_temp_%d.csv" \
                                                            % (time.time()))
        self.__create_output_file__(temp_path)
        return temp_path

    def __get_output_path__(self, user_output_path):
        tmct_output_path = None
        if user_output_path is None:
            tmct_output_path = os.path.join("./telomerecat_length_%d.csv" \
                                                            % (time.time(),))
        else:
            tmct_output_path = user_output_path

        return tmct_output_path

    def __create_output_file__(self,output_csv_path):
        with open(output_csv_path,"w") as total:
            header = "Sample,F1,F2,F4,Psi,Insert_mean,Insert_sd,Read_length\n"
            total.write(header)
        return output_csv_path

    def __write_to_csv__(self,
                         read_type_counts,
                         vital_stats,
                         output_csv_path,
                         name):
        with open(output_csv_path,"a") as counts:
            counts.write("%s,%d,%d,%d,%.3f,%.3f,%.3f,%d\n" %\
                (name,
                read_type_counts["F1"],
                read_type_counts["F2"],
                read_type_counts["F4"],
                read_type_counts["sample_variance"],
                vital_stats["insert_mean"],
                vital_stats["insert_sd"],
                vital_stats["read_len"]))

    def get_parser(self):
        parser = self.default_parser()
        parser.description = textwrap.dedent(
        '''\
        %s
        %s

            The telbam2length command allows the user to genereate a telomere
            length estimate from a previously generated TELBAM file.

            Example useage:

            telomerecat telbam2length /path/to/some_telbam.bam

            This will generate a .csv file with an telomere length estimate
            for the `some_telbam.bam` file.

        %s
        ''' % (self.instance_name,
               self.header_line,
               self.header_line,))

        parser.add_argument('input', metavar='TELBAM(S)', nargs='+',
            help="The TELBAM(s) that we wish to analyse")
        parser.add_argument('--output',
            metavar='CSV',type=str, nargs='?', default=None,
            help=('Specify output path for length estimation CSV.\n'
                  'Automatically generated if left blank [Default: None]'))
        parser.add_argument('-s',type=int, nargs='?', default=10000
            ,help="The amount of reads considered by each\n"
                    "distributed task. [Default: 10000]")
        parser.add_argument('-t',"--trim", type=int, nargs='?', default=0
            ,help="Use only the amount of sequence specified by this  \n"
                  "option (i.e if the value 90 is supplied\n"
                  "then only the first 90 bases are\n"
                  "considered) [Default: Whole read]")

        return parser

if __name__ == "__main__":
    print "Do not run this script directly. Type `telomerecat` for help."
