#!/usr/bin/env python
#Once upon a time...

#Create a bam file with all the putatative telomere reads

import sys,os
import argparse
import textwrap
import pysam
import time
import gc
import pdb

import parabam
import error_estimator

import Queue as Queue2

from multiprocessing import Queue,Process
from itertools import izip

######################################################################
##
##      Create telbams for list of given BAM files.
##
##      Author: jhrf
##
######################################################################

#Here we define what constitutes a telomereic read.
#If both reads in a pair match the critea we store
#them to the telbam
def rule(reads, constants, master):
    tel_pats = constants["tel_pats"]
    thresh = constants["thresh"]
    telomere_status = False

    for read in iter(reads):
        for pattern in tel_pats:
            if read.seq.count(pattern) > thresh:
                telomere_status = True
                break

    results = []
    if telomere_status:
        results = [("telbam",reads[0]),("telbam",reads[1])]
    return results

#This whole package is essentially just a wrapper for a call to parabam subset
class Bam2Telbam(parabam.core.Interface):
    '''Interface to interact with telomerecat programatically. This interface
    is also called by the `telomerecat` script'''
    def __init__(self,
                 temp_dir=None,
                 task_size=250000,
                 total_procs=8,
                 reader_n=2,
                 verbose=False,
                 announce=True,
                 cmd_run=False):

        super(Bam2Telbam,self).__init__(instance_name = \
                                            "telomerecat bam2telbam", 
                                        temp_dir=temp_dir,
                                        task_size=task_size,
                                        total_procs=total_procs,
                                        reader_n=reader_n,
                                        verbose=verbose,
                                        announce = announce,
                                        cmd_run=cmd_run)

    def run_cmd(self):
        """Called from the master `telomerecat` script which handels the
        cmd line interface. Requires an argparse parser. Users should call
        the `run` function."""
        self.run(input_paths = self.cmd_args.input)

    def run(self,input_paths, keep_in_temp = False):
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

        subset_types=["telbam"]
        tel_pats = ["TTAGGG","CCCTAA"]

        #need to define my constants and engine here:
        telbam_constants = {"thresh":2,
                            "tel_pats":tel_pats}

        final_output_paths = {}

        for input_path in input_paths:

            if self.verbose:
                sys.stdout.write(" Generating TELBAM from: %s\n" % (input_path,))
                sys.stdout.write("\t- TELBAM generation started %s\n" %\
                                    (self.__get_date_time__(),))

            subset_interface = parabam.Subset(temp_dir=self.temp_dir,
                                              total_procs=self.total_procs,
                                              task_size=self.task_size,
                                              reader_n=self.reader_n,
                                              verbose=self.verbose,
                                              pair_process=True,
                                              include_duplicates=True,
                                              keep_in_temp=keep_in_temp)
            #call to parabam subset
            telbam_paths = subset_interface.run(input_paths=[input_path],
                                                subsets=subset_types,
                                                constants = telbam_constants,
                                                rule = rule)

            if self.verbose:
                sys.stdout.write("\t- Adding error estimation to TELBAM header\n")

            error_estimator.add_error_to_telbam_header(input_path, 
                                                       telbam_paths[input_path]['telbam'], 
                                                       self.temp_dir)

            if self.verbose:
                sys.stdout.write("\t- TELBAM generation finished %s\n\n"\
                                     % (self.__get_date_time__(),))

            gc.collect()
            final_output_paths.update(telbam_paths)

        self.__goodbye__()
        return final_output_paths

    def get_parser(self):
        parser = self.default_parser()
        parser.description = textwrap.dedent('''\
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
         ''' % (self.instance_name, 
               self.header_line, 
               self.header_line,))

        parser.add_argument('input',metavar='BAM(S)', nargs='+',
            help='BAM file(s) for which we wish to generate TELBAMS')

        return parser


if __name__ == "__main__":
    print "Type telomerecat -h for help!"

#....happily ever after.
