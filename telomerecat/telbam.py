from typing import List
import os
import pysam
import logging
import tempfile
from easytimer import tick, tock

from telomerecat.constants import TEL_PATS, HTS_EXT_TO_AF_MODE
from pysam import AlignmentFile

def _log_setup(loglevel):
    logging.basicConfig(level=getattr(logging, loglevel.upper()), format='%(levelname)s: %(message)s')

def get_align_file(path:str, mode='r', template=None, expectIndex=True, threads=1):
    hts_ext = os.path.splitext(path)[-1]
    xam_type = HTS_EXT_TO_AF_MODE[hts_ext]
    if expectIndex is False:
        save = pysam.set_verbosity(0)
    af = pysam.AlignmentFile(path, f'{mode}{xam_type}', template=template, threads=threads)
    if expectIndex is False:
        pysam.set_verbosity(save)
    return af

def get_processes(processes):
    if processes == -1:
        processes = len(os.sched_getaffinity(0))
    return processes

def get_hts_processes(processes):
    if processes > 2:
        processes = 2
    return processes

def telbam_path(dest_dir:str, xam_file:str):
    (base, ext) = os.path.splitext(os.path.basename(xam_file))
    return os.path.join(dest_dir, f'{base}_telbam{ext}')

def collate_pairs(tmpdir: str, xam_file: str, processes: int):
    pairs = os.path.join(tmpdir, 'pairs.bam')
    collate_opts = ['--no-PG', '-f', '-n', '1', '-r', '5000000']  # ~8GB RAM
    if processes > 1:
        collate_opts.extend(['-@', str(processes)])
    collate_opts.extend(['-o', pairs, xam_file, os.path.join(tmpdir, 'collate')])
    pysam.collate(*collate_opts) # don't forget to splat, this isn't perl
    return pairs

def pairs_to_telbam(af_pairs:AlignmentFile, af_telbam:AlignmentFile):
    read_iter = af_pairs.fetch(until_eof=True)
    while True:
        read_a = next(read_iter, None)
        if read_a is None:
            break
        read_b = next(read_iter)
        qseq = read_a.query_sequence
        if TEL_PATS[0] in qseq or TEL_PATS[1] in qseq:
            af_telbam.write(read_a)
            af_telbam.write(read_b)
        else:
            qseq = read_b.query_sequence
            if TEL_PATS[0] in qseq or TEL_PATS[1] in qseq:
                af_telbam.write(read_a)
                af_telbam.write(read_b)
    return


def process_alignments(outbam_dir:str, processes:int, tmpdir:str, files:List[str], quick=False, verbose=0):
    if verbose > 0:
        # for compatibility with original telomerecat
        _log_setup('DEBUG')
    telbam_paths = {}
    # general max CPUs (works correctly for cgroups)
    processes = get_processes(processes)
    # max threads for HTS actions
    hts_processes = get_hts_processes(processes)

    for xam_file in files:
        with tempfile.TemporaryDirectory(prefix='telomerecat_', dir=tmpdir) as tmpdir:
            # could parallel on contig if wanted but more difficult to pull reads together

            if quick is True:
                pairs = xam_file
            else:
                logging.debug('Collating read-pairs')
                tick('Collating read-pairs')
                pairs = collate_pairs(tmpdir, xam_file, hts_processes)

            logging.debug('Collecting telomere candidates')
            tick('Collecting telomere candidates')

            telbam = telbam_path(outbam_dir, xam_file)
            telbam_paths[xam_file] = {'telbam': telbam}

            af_pairs = get_align_file(pairs, mode='r', expectIndex=False, threads=hts_processes)
            af_telbam = get_align_file(telbam, mode='w', template=af_pairs)
            pairs_to_telbam(af_pairs, af_telbam)

            tock()
    return telbam_paths
