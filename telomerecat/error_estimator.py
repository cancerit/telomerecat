#! /usr/bin/env python

import sys
import os
import pdb
import pysam
import time
import random
from collections import Counter
import numpy as np

class ReadRecord(object):
    def __init__(self, read_id, base_id, base):
        self.read_id = read_id
        self.base_id = base_id
        self.base = base

class BaseRecord(object):
    def __init__(self, nucleotide, count, freq):
        self.nucleotide = nucleotide
        self.count = count
        self.freq = freq

class Base(object):

    def __init__(self):
        self.observations = 0
        self.read_records = []

        self.nucleotide_counts = Counter()
        self.most_common_freq = 0

    def add(self, nucleotide, read_record):
        self.nucleotide_counts[nucleotide] += 1
        self.observations += 1

        self.read_records.append(read_record)

        top_base,top_count, freq = self.get_most_common()
        self.most_common_freq = freq

    def get_most_common(self):
        if self.observations > 0:
            top_base,top_count = self.nucleotide_counts.most_common()[0]
            freq = (float(top_count) / self.observations)
            return top_base, top_count, freq
        else:
            return "N",0,0

    def get_sorted_base_records(self):
        records = []
        for nucleotide,count in self.nucleotide_counts.items():
            freq = float(count) / self.observations
            records.append(BaseRecord(nucleotide, count, freq))

        records.sort(key= lambda x: x.count)
        return records

class Sequence(object):
    def __init__(self, reads):
        self._compliments = {"A":"T",
                     "C":"G",
                     "G":"C",
                     "T":"A",
                     "N":"N"}
        
        self.reads = reads
        self.sequence = self.sequence_setup(reads)

    def sequence_setup(self, reads):
        sequence = []
        offset = reads[0].pos

        for read_id, read in enumerate(reads):
            seq_gen = self.__sequence_generator__(read)

            for mapped_index, true_index, base in seq_gen:
                sequence_index = self.__get_offset_locus__(read,
                                                           mapped_index,
                                                           offset)
                self.__prepare_sequence__(sequence, sequence_index)

                if sequence_index >= 0:
                    read_record = ReadRecord(read_id, true_index, base)
                    sequence[sequence_index].add(base, read_record)

        return sequence

    def __sequence_generator__(self, read):

        if read.is_unmapped:
            return

        cigar_list = self.__get_cigar_list__(read)

        mapped_index = 0
        mapped_offset=0

        if read.cigar[0][0] != 0:
            mapped_offset -= read.cigar[0][1]

        for true_index, base, cigar_entry in zip(xrange(len(read.seq)),
                                     read.seq,
                                     cigar_list):

                yield mapped_index+mapped_offset, true_index, base
                mapped_index += 1
        return

    def get_consensus_sequence(self, thresh = .8):
        consensus = []
        for base in self.sequence:
            base_char, count, freq = base.get_most_common()

            record_char = base_char.upper()
            if freq < thresh:
                record_char = base_char.lower()
            consensus.append(record_char)

        return consensus

    def print_consensus_sequence(self):
        consensus_sequence = self.get_consensus_sequence()
        print "".join(consensus_sequence)

    def __get_cigar_list__(self, read):
        cigar_list = []
        offset = 0

        cigar_iter = iter(read.cigar)

        for operation,length in cigar_iter:
            for _ in xrange(offset, offset+length):
                cigar_list.append(operation == 0)
            offset = offset + length
        return cigar_list

    def print_bases(self):
        for i,base in enumerate(self.sequence):
            top_base,top_count,freq = base.get_most_common()
            print "%d: %s %d %.1f%%" % (i, top_base, top_count, freq*100), \
                  dict(base.nucleotide_counts)

    def __prepare_sequence__(self, sequence, sequence_index):
        while len(sequence) <= sequence_index:
            sequence.append(Base())

    def __get_offset_locus__(self, read, base_index, offset):
        return (read.pos - offset) + base_index

    def print_base_segments(self, base_id, surround = 3):
        for record in self.sequence[base_id].read_records:
            read = self.reads[record.read_id]
            base_id = record.base_id
            left_flank = read.seq[base_id-surround:base_id]
            right_flank = read.seq[base_id+1:base_id+1+surround]

            left_flank_qual = read.qual[base_id-surround:base_id]
            right_flank_qual = read.qual[base_id+1:base_id+1+surround]

            print record.read_id,left_flank, read.seq[base_id], right_flank
            print "  ",left_flank_qual, read.qual[base_id], right_flank_qual
            print ""

    def get_error_breakdown(self):
        consensus_sequence = self.get_consensus_sequence()

        test_bases = self.__test_bases__

        nucleotide_changes = Counter()

        prev_base = self.sequence[0]
        base = self.sequence[1]

        for i, next_base in enumerate(self.sequence[2:]):
            if test_bases(prev_base, base, next_base):

                base_records = base.get_sorted_base_records()
                most_common = base_records[-1]

                for record in base_records:
                        if record.nucleotide != most_common.nucleotide:
                            error_key = \
                                "%s%s>%s%s" % (consensus_sequence[i-1],
                                               most_common.nucleotide,
                                               record.nucleotide,
                                               consensus_sequence[i+1])
                            nucleotide_changes[error_key] += 1

                        for read_record in base.read_records:
                            if read_record.base == record.nucleotide:
                                pass
        return nucleotide_changes

    def __test_bases__(self,prev_base, base, next_base, thresh = .8):
        return base.observations >= 5 and \
                base.most_common_freq != 1  and \
                 base.most_common_freq >= thresh and \
                  prev_base.most_common_freq >= thresh and \
                    prev_base.most_common_freq >= thresh

    def populate_error_matricies(self, correct, error):

        for i, base in enumerate(self.sequence):
            if base.observations >=5 and \
                base.most_common_freq >= .8 and \
                 base.most_common_freq != 1:

                base_records = base.get_sorted_base_records()
                most_common = base_records[-1]

                correct_loci = []
                error_count = 0
                for read in base.read_records:
                    if read.base != most_common.nucleotide:
                        error[0,read.base_id] += 1
                        error_count += 1
                    else:
                        correct_loci.append(read.base_id)

                for i in random.sample(correct_loci, error_count):
                    correct[0,i] += 1

    def error_ratio(self):
        error_bases = 0
        for base in self.sequence:
            if base.observations >= 5:
                if base.most_common_freq >= .8 and base.most_common_freq != 1:
                    correct_reads = base.observations * base.most_common_freq
                    error_bases += base.observations - correct_reads 

        return error_bases

class ScanBam(object):

    def scan(self, path):
        bam_file = pysam.AlignmentFile(path,"rb")
        references = zip(bam_file.references,bam_file.lengths)

        sample_count = 0

        bad_base_ratios = []
        read_counts = []

        while sample_count < 300:
            chrm_name, chrm_len = random.sample(references, 1)[0]
            if chrm_len > (10 ** 7 * 5):
                random_location = random.randint(0, chrm_len-5000000)

                region = "%s:%d-%d" % (chrm_name,
                                       random_location,
                                       random_location + 5000)

                reads = self.get_reads(bam_file, region)
                if len(reads) > 100:

                    read_len = len(reads[0].seq)

                    read_counts.append(len(reads))
                    sample_count += 1
                    sequence = Sequence(reads)

                    unmapped = 0
                    for read in sequence.reads:
                        unmapped += int(read.is_unmapped)

                    error_bases = sequence.error_ratio()

                    bad_base_ratio = (len(reads)*read_len) / float(error_bases+(unmapped*read_len))
                    bad_base_ratios.append(bad_base_ratio)

        bam_file.close()
        return read_counts, bad_base_ratios

    def get_reads(self, bam_file, region):
        all_reads = []
        for read in bam_file.fetch(region = region):
            all_reads.append(read)
        return all_reads

def modify_telbam_header(telbam_path, 
                         read_counts, 
                         bad_base_ratios, 
                         temp_dir_path):

    telbam_file = pysam.AlignmentFile(telbam_path, "rb")
    telbam_header = telbam_file.header
    if "CO" not in telbam_header.keys():
        telbam_header["CO"] = []

    error_est_string = "telomerecat_error_ratio:%.3f,%.3f,%.3f" % (np.mean(bad_base_ratios),
                                                                   np.median(bad_base_ratios),
                                                                   np.std(bad_base_ratios))
    error_count_string = "telomerecat_error_counts:%.3f,%.3f,%.3f" % (np.mean(read_counts),
                                                                   np.median(read_counts),
                                                                   np.std(read_counts))
    telbam_header["CO"].extend((error_est_string, error_count_string))

    temp_telbam_path = os.path.join(temp_dir_path,"temp_dir_path_%d.bam" % (time.time(),))
    temp_telbam_file = pysam.AlignmentFile(temp_telbam_path,"wb",header=telbam_header)

    for read in telbam_file.fetch(until_eof=True):
        temp_telbam_file.write(read)

    temp_telbam_file.close()
    telbam_file.close()

    os.rename(temp_telbam_path, telbam_path)

def add_error_to_telbam_header(bam_path, telbam_path, temp_dir_path):
    read_counts,bad_base_ratios = ScanBam().scan(bam_path)

    modify_telbam_header(telbam_path, 
                         read_counts, 
                         bad_base_ratios, 
                         temp_dir_path)

    #sequence.print_consensus_sequence()
    #sequence.print_base_segments(874)

    #sequence.get_error_breakdown()

if __name__ == "__main__":

    if len(sys.argv) >= 2:
        for path in sys.argv[1:]:
            add_error_to_telbam_header(path, None)
    else:
        print "Useage:\n\terror_model.py <bam_paths...>"
