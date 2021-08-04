from typing import List
import os
import sys
import pysam
import logging
import tempfile
import multiprocessing as mp
import subprocess
from functools import partial

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
    (base, _) = os.path.splitext(os.path.basename(xam_file))
    return os.path.join(dest_dir, f'{base}_telbam.bam')


def collate_pairs(xam_file: str, tmpdir:str, processes=1, reference=None):
    pairs = os.path.join(tmpdir, 'pairs.bam')
    collate_opts = ['--no-PG', '-f', '-n', '1', '-l', '1', '-r', '5000000']  # ~8GB RAM
    if processes > 1:
        collate_opts.extend(['-@', str(processes)])
    if reference:
        collate_opts.extend(['--reference', reference])
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


def to_telbam(xam_file:str, outbam_dir:str, tmpdir:str, hts_processes=1, reference=None):
    telbam = telbam_path(outbam_dir, xam_file)
    telbam_paths = {xam_file: {'telbam': telbam}}

    with tempfile.TemporaryDirectory(prefix='telomerecat_', dir=tmpdir) as this_tmp:
        collate_wrap = ['pysam_collate', '-t', str(this_tmp), '-@', str(hts_processes)]
        if reference:
            collate_wrap.extend(['--reference', reference])
        collate_wrap.append(xam_file)
        named_fifo = os.path.join(this_tmp, 'collatepipe.bam')
        os.mkfifo(named_fifo)
        collate_wrap.extend(['--output', named_fifo])
        collate_proc = subprocess.Popen(collate_wrap)
        af_pairs = get_align_file(named_fifo, mode='r', expectIndex=False, threads=hts_processes)
        af_telbam = get_align_file(telbam, mode='w', template=af_pairs, threads=hts_processes)
        pairs_to_telbam(af_pairs, af_telbam)

        # we should be finished anyway
        (out_proc, err_proc) = collate_proc.communicate(timeout=15)
        collate_exit_code = collate_proc.returncode
        if collate_exit_code != 0:
            print(f'ERROR: samtools collate process exited with a non-zero value ({collate_exit_code})')
            sys.exit(collate_exit_code)

    return telbam_paths


def process_alignments(outbam_dir:str, processes:int, tmpdir:str, files:List[str], reference=None, verbose=0):
    if verbose > 0:
        # for compatibility with original telomerecat
        _log_setup('DEBUG')
    telbam_paths = {}
    # general max CPUs (works correctly for cgroups)
    processes = get_processes(processes)
    # max threads for HTS actions
    hts_processes = get_hts_processes(processes)

    results = []
    ctx = mp.get_context('forkserver')
    logging.debug(f'Starting telbam generation for {len(files)} inputs')

    parallel_processes = int(processes/hts_processes)
    # allow non-multiprocess use when minimal files, helpful for debugging
    if parallel_processes <= 1 or len(files) == 1:
        logging.debug(f'Inline execution of telbam generation, {hts_processes} threads per process')
        # allow max threads for helpers
        hts_processes = processes
        for xam_file in files:
            telbam_path = to_telbam(xam_file, outbam_dir=outbam_dir, tmpdir=tmpdir, hts_processes=hts_processes, reference=reference)
            telbam_paths.update(telbam_path)
    else:
        logging.debug(f'Multiprocessing execution of telbam generation, {hts_processes} threads per process, {parallel_processes} in parallel')
        with ctx.Pool(processes=parallel_processes) as pool:
            results = pool.map(partial(to_telbam, outbam_dir=outbam_dir, tmpdir=tmpdir, hts_processes=hts_processes, reference=reference), files)
        for telbam_path in results:
            telbam_paths.update(telbam_path)

    return telbam_paths
