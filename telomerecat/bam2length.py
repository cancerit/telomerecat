#!/usr/bin/env python
#Once upon a time...

import textwrap
import parabam

from telomerecat.core import TelomerecatInterface

######################################################################
##
##      Genereate length estimates for given BAM files
##
##      This simply strings together bam2telbam 
##      and telbam2length packages
##        
##      Author: jhrf
##
######################################################################

class Bam2Length(TelomerecatInterface):
    def __init__(self,temp_dir=None,
                 task_size=250000,
                 total_procs=8,
                 reader_n=1,
                 verbose=False,
                 announce=True,
                 cmd_run=False):

        super(Bam2Length,self).__init__(instance_name="telomerecat bam2length",
                                         temp_dir=temp_dir,
                                         task_size=task_size,
                                         total_procs=total_procs,
                                         reader_n=reader_n,
                                         verbose=verbose,
                                         announce=announce,
                                         cmd_run=cmd_run)

    def run_cmd(self):
        self.run(input_paths = self.cmd_args.input,
                 output_path = self.cmd_args.output,
                 inserts_path = self.cmd_args.insert,
                 discard_telbams = self.cmd_args.discard_telbams)

    def run(self, input_paths, 
                  output_path=None, 
                  inserts_path=None,
                  discard_telbams=False):
        
        #Import here to avoid infinite loop on import
        from telomerecat import Bam2Telbam
        from telomerecat import Telbam2Length

        self.__introduce__()

        telbam_interface = Bam2Telbam(temp_dir=self.temp_dir,
                                      total_procs=self.total_procs,
                                      task_size=self.task_size,
                                      reader_n=self.reader_n,
                                      verbose=self.verbose,
                                      announce=False)

        out_files = telbam_interface.run(input_paths=input_paths, 
                                         keep_in_temp = discard_telbams)

        length_paths = self.collapse_out_files(out_files)
        length_interface = Telbam2Length(temp_dir=self.temp_dir,
                                         total_procs=self.total_procs,
                                         task_size=self.task_size,
                                         reader_n = self.reader_n,
                                         announce=False,
                                         verbose=self.verbose)

        length_interface.run(input_paths=length_paths,
                             output_path=output_path,
                             inserts_path=inserts_path)

        self.__goodbye__()
        self.interface_exit()

    def collapse_out_files(self,out_files):
        length_paths = []
        for path in out_files.keys():
            length_paths.append(out_files[path]["telbam"])
        return length_paths

    def get_parser(self):
        parser = self.default_parser()
        parser.description = textwrap.dedent(
        '''\
        %s
        %s

            The bam2length command allows the user to genereate a telomere
            length estimate from a BAM file.

            The majority of the time taken running the bam2length script is
            spent collecting telomeric reads from the BAM file. If you wish
            to run multiple analyses on this file be sure to keep the TELBAMS
            that this run creates by using the `-k` option.

            If you wish to generate TELBAMS seperately from length estimation
            you should use the bam2telbam command.
            
            Type `telomerecat bam2telbam` to find out more.

        %s
        ''' % (self.instance_name, 
               self.header_line, 
               self.header_line,))        

        parser.add_argument('input',metavar='BAM(S)', nargs='+',
            help=('BAM file(s) for which we wish to\n'
                  'generate telomere length estimates'))
        parser.add_argument('-x','--discard_telbams',
            action="store_true",default=True,
            help=('The program will NOT save any TELBAMs\n'
                   'generated as part of the analysis'))

        return parser

if __name__ == "__main__":
    print "Please do not run this script directly."\
     " Type telomerecat -h for more information."


#....happily ever after.