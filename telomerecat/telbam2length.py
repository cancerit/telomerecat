import sys
import textwrap
import time
import os
import re
import parabam
import pdb
import pysam

import numpy as np
import pandas as pd

from shutil import copy
from itertools import izip
from collections import namedtuple,Counter
from pprint import pprint

from telomerecat import Csv2Length
from telomerecat.core import TelomerecatInterface
from telomerecat import readmodel

######################################################################
##
##      Create a length estimate given a set of TELBAMS 
##
##      Author: jhrf
##
######################################################################

class TeloReadFactory(object):

    def __init__(self, allow = 0, trim_reads=0):
        self._TeloRead = namedtuple("SimpleRead","seq qual" +
                                      " five_prime pattern mima_loci n_loci"+
                                      " avg_qual")

        self._allow = allow
        self._trim_reads = trim_reads

        self._compliments = {"A":"T","T":"A",
                             "C":"G","G":"C",
                             "N":"N"}

        self.mima_logic = MismatchingLociLogic()

    def get_telo_read(self,read):
        seq,qual = self.__flip_and_compliment__(read)
        if self._trim_reads > 0:
            seq,qual = (seq[:self._trim_reads],
                        qual[:self._trim_reads])
        return self.sequence_to_simpleread(seq,qual)

    def sequence_to_simpleread(self, seq, qual):
        (mima_loci,frameshift_loci),pattern = \
                self.mima_logic.get_telo_mismatch(seq)


        all_loci, avg_qual,n_loci = self.__get_phred_score__(qual, 
                                                   mima_loci, 
                                                   frameshift_loci)

        all_loci, avg_qual, n_loci = self.adjust_mima_for_thresh(all_loci, 
                                                                qual, 
                                                                n_loci)

        simple_read = self._TeloRead(
            seq,
            qual,
            self.__get_five_prime__(pattern),
            pattern,
            all_loci,
            n_loci,
            avg_qual)

        return simple_read

    def adjust_mima_for_thresh(self, all_loci, qual, n_loci):

        new_loci = []
        new_phreds = []
        new_n_loci = 0

        qual_bytes = [ord(b) for b in qual]
        read_qual_mean = np.mean(qual_bytes)
        read_qual_std = np.std(qual_bytes)

        for i in all_loci:

            if qual_bytes[i] > (read_qual_mean - read_qual_std):
                new_n_loci+=1
                new_loci.append(i)
                new_phreds.append(qual_bytes[i])

        if new_n_loci == 0:
            new_phreds_mean = 0
        else:
            new_phreds_mean = np.mean(new_phreds)

        return new_loci, int(new_phreds_mean), new_n_loci

    def __get_phred_score__(self, qual, mima_loci, frameshift_loci):
        if len(mima_loci) + len(frameshift_loci) == 0:
            return [],0, 0

        remove_mimas = []
        phreds = []
        all_loci = []

        for start, end in frameshift_loci:
            fuse_phreds = []

            for i in xrange(start,end):
                if i in mima_loci:
                    remove_mimas.append(i)
                fuse_phreds.append( (i,ord(qual[i])) )

            min_loci, min_phred = min(fuse_phreds,key=lambda x : x[1])

            phreds.append(min_phred)
            all_loci.append(min_loci)

        for locus in mima_loci:
            if locus not in remove_mimas:
                phreds.append(ord(qual[locus]))
                all_loci.append(locus)

        all_loci.sort()

        return all_loci, np.mean(phreds), len(phreds)

    def __get_average_qual__(self,qual,mima_loci):
        if len(mima_loci) == 0:
            return 0
        phreds = np.array([ord(qual[i]) for i in mima_loci])
        return np.mean(phreds)

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

class SampleStatsFinder(object):

    def __init__(self,temp_dir,total_procs,task_size,trim_length=0):
        self.temp_dir = temp_dir
        self._total_procs = total_procs
        self._task_size = task_size
        self._trim_length = trim_length

    def get_sample_stats(self,sample_path):

        read_stats_csv_path = self.__run_read_stat_rule__(sample_path)
        sample_stats = self.csv_to_sample_stats(read_stats_csv_path)

        return sample_stats

    def csv_to_sample_stats(self, read_stats_csv_path):
        read_stats = pd.read_csv(read_stats_csv_path).values

        sample_stats = {}
        self.__load_insert_dat__(sample_stats, read_stats)
        sample_stats["read_len"] = round(read_stats[:,8].mean(),1)
        sample_stats["initial_read_len"] = round(read_stats[:,9].mean(),1)
        sample_stats["min_phred"] = int(read_stats[:,6].min())
        sample_stats["max_phred"] = int(read_stats[:,7].min())

        sample_stats["read_stats_path"] = read_stats_csv_path

        return sample_stats

    def __load_insert_dat__(self, sample_stats, read_stats):
        relevant_mask = np.arange(read_stats.shape[0]) % 2 == 0
        insert_stats = read_stats[relevant_mask,:]

        proper_pair_mask = read_stats[:,3] == 1
        insert_stats = read_stats[proper_pair_mask,:]
        
        well_map_thresh = np.percentile(insert_stats[:,4],66)
        well_map_mask = insert_stats[:,4] >= well_map_thresh

        well_mapped_reads = insert_stats[well_map_mask,:]

        sample_stats["insert_mean"] = round(well_mapped_reads[:,5].mean(),3)
        sample_stats["insert_std"] = round(well_mapped_reads[:,5].std(),3)

    def __run_read_stat_rule__(self,sample_path,keep_in_temp=True):

        telo_read_factory = TeloReadFactory()#allow=self._allow,
                                            #trim_reads=self._trim)


        def telomere_heuristic(telo_read):
            prev = None
            distances = []
            for i in telo_read.mima_loci:
                if prev is not None:
                    distances.append(i - prev)
                prev = i
            if len(distances) == 0:
                distances.append(100)
            return np.mean(distances)

        def prepare_stats(telo_read, read):
            qual_bytes = [ord(q) for q in telo_read.qual]

            stats = (telo_read.n_loci,
                     int(telo_read.five_prime),
                     telo_read.avg_qual,
                     int(read.is_proper_pair),
                     read.mapq,
                     abs(read.template_length),
                     min(qual_bytes),
                     max(qual_bytes),
                     len(telo_read.seq),
                     len(read.seq),
                     telomere_heuristic(telo_read))

            return stats

        def get_read_entries(reads):
            telo_reads = \
                [telo_read_factory.get_telo_read(read) for read in reads]

            print_check = [t.n_loci == 10 for t in telo_reads]
            if any(print_check):
                for t in telo_reads:
                    telo_read_factory.mima_logic.print_mima(t.seq,
                                                            t.qual,
                                                            t.pattern)
                print "--"
            read_1_stats = prepare_stats(telo_reads[0],reads[0])
            read_2_stats = prepare_stats(telo_reads[1],reads[1])

            read_1_entry = list(read_1_stats)
            read_1_entry.extend(read_2_stats)

            read_2_entry = list(read_2_stats)
            read_2_entry.extend(read_1_stats)

            return_dat = np.zeros((2,len(read_1_entry)))

            return_dat[0,:] = read_1_entry
            return_dat[1,:] = read_2_entry

            return return_dat


        def rule(reads, constants, master):
            stats = {"read_stats":get_read_entries(reads)}
            return stats

        structures = {"read_stats":{"data":np.zeros((2,22)),
                                    "store_method":"vstack"}}

        stat_interface = parabam.Stat(temp_dir=self.temp_dir,
                                      total_procs = self._total_procs,
                                      task_size = 7000,
                                      pair_process = True,
                                      keep_in_temp=keep_in_temp)

        out_paths = stat_interface.run(input_paths= [sample_path],
                                       constants = {},
                                       rule = rule,
                                       struc_blueprint = structures)

        return out_paths[sample_path]["read_stats"]

class ReadStatsFactory(object):

    def __init__(self,temp_dir,
                      total_procs=4,
                      task_size=5000,
                      trim_reads=0,
                      allow=0,
                      debug_print=False):

        self.temp_dir = temp_dir
        self._total_procs = total_procs
        self._task_size = task_size

        self._debug_print = debug_print
        self._trim = trim_reads
        self._allow = allow

    def get_read_counts(self, path, sample_stats):
        read_stats = self.__path_to_read_stats__(
            sample_stats["read_stats_path"])

        # read_model = readmodel.TelomereReadModel(sample_stats,
        #                                         read_stats=read_stats).run()
        read_model = readmodel.get_read_model(sample_stats,
                                              read_stats,
                                              self._total_procs)

        read_counts = self.read_stats_to_counts(read_stats, read_model)

        pdb.set_trace()
        # self.__delete_analysis_paths__(read_stats_path)
        return read_counts

    def dist_matrix(self, read_model, exclude_dists, count):

        results = []
        for _ in xrange(count):
            results.append(read_model.get_dist(exclude_dists=exclude_dists))

        return np.array(results).reshape(count, -1)

    def __delete_analysis_paths__(self, read_stat_paths):
        for analysis, path in read_stat_paths.items():
            os.remove(path)

    def __path_to_read_stats__(self, read_stats_path):
        return pd.read_csv(read_stats_path, header=None).values

    def read_stats_to_counts(self, read_stats, read_model):

        complete_reads, boundary_reads = \
            self.__get_complete_status__(read_stats, read_model)

        f2_count,f4_count = self.__get_boundary_counts__(boundary_reads)
        f1_count = self.__get_f1_count__(complete_reads)

        return_dat = { "F2": int(f2_count),
                       "F1": int(f1_count),
                       "F4": f4_count,
                       "sample_variance": 0}

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
        
    def __get_complete_status__(self,read_stats, read_model):
        boundary_indicies = []
        complete_indicies = []

        telo_est = self.dist_matrix(read_model, ["subtelo", "nontelo"], 100)
        subtelo_est = self.dist_matrix(read_model, ["atelo", "nontelo"], 100)

        for i in xrange(int(read_stats.shape[0])):
            read_info = map(int,read_stats[i,[0,-2]])
            pair_info = map(int,read_stats[i,[2,-1]])

            read = error_profile[read_info[0],read_info[1]] 
            pair = error_profile[pair_info[0],pair_info[1]]

            if read and pair:
                complete_indicies.append(i)
            elif (not read) and pair:
                boundary_indicies.append(i)

        return read_stats[complete_indicies,:],\
                read_stats[boundary_indicies,:]

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
        self.run(input_paths=self.cmd_args.input,
                 trim=self.cmd_args.trim,
                 allow=self.cmd_args.allow,
                 inserts_path=self.cmd_args.insert,
                 correct_f2a=not self.cmd_args.disable_correction,
                 output_path=self.cmd_args.output)

    def run(self, input_paths,
                  trim=0,
                  allow=0,
                  output_path=None,
                  correct_f2a=True,
                  inserts_path=None):
        
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

        sample_stats_finder = SampleStatsFinder(self.temp_dir,
                                                self.total_procs,
                                                self.task_size,
                                                trim)

        for sample_path,sample_name, in izip(input_paths,names):
            sample_intro = "\t- %s | %s\n" % (sample_name, 
                                              self.__get_date_time__())

            self.__output__(sample_intro,2)

            sample_stats = sample_stats_finder.get_sample_stats(sample_path)

            read_type_counts = self.__get_read_types__(sample_path,
                                                       sample_stats,
                                                       self.total_procs,
                                                       trim,
                                                       allow)

            f2a = read_type_counts["F2"] - read_type_counts["F4"]
            length = (float(read_type_counts["F1"]) / f2a) * vital_stats["insert_mean"]
            print "F1:%d,F2a:%d,L:%d" % (read_type_counts["F1"],f2a,length,)

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

    def __get_read_types__(self,sample_path,
                                sample_stats,
                                total_procs,
                                trim,
                                allow):

        read_stats_factory = ReadStatsFactory(temp_dir=self.temp_dir,
                                              total_procs=total_procs,
                                              trim_reads=trim,
                                              allow=allow,
                                              debug_print=False)

        read_type_counts = read_stats_factory.get_read_counts(sample_path,
                                                              sample_stats)
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
            header = ("Sample,F1,F2,F4,Psi,Insert_mean,Insert_sd,"
                      "Read_length,Initial_read_length\n")
            total.write(header)
        return output_csv_path

    def __write_to_csv__(self,
                         read_type_counts,
                         vital_stats,
                         output_csv_path,
                         name):
        with open(output_csv_path,"a") as counts:
            counts.write("%s,%d,%d,%d,%.3f,%.3f,%.3f,%d,%d\n" %\
                (name,
                read_type_counts["F1"],
                read_type_counts["F2"],
                read_type_counts["F4"],
                read_type_counts["sample_variance"],
                vital_stats["insert_mean"],
                vital_stats["insert_sd"],
                vital_stats["read_len"],
                vital_stats["initial_read_len"]))

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
        parser.add_argument('-a','--allow',metavar="PERCENT",
                            type=int, nargs='?', default=5
            ,help="A threshold on the `genuine` mismatches to allow\n"
                  "in each seqeuncing read. Expressed as a percentage of\n"
                  "read total [Default: 0]")
        parser.add_argument('-t',"--trim", type=int, metavar="INT",nargs='?', default=0
            ,help="Use only the amount of sequence specified by this  \n"
                  "option (i.e if the value 90 is supplied\n"
                  "then only the first 90 bases are\n"
                  "considered) [Default: Whole read]")

        return parser

if __name__ == "__main__":
    print "Do not run this script directly. Type `telomerecat` for help."
