#!/usr/bin/env python3
import os
import pkg_resources
import pysam
import click

@click.command()
@click.version_option(pkg_resources.require("telomerecat")[0].version)
@click.option('--output', required=True, type=str, help='File or named pipe')
@click.option('-@', '--threads', required=False, type=int, default=1, show_default=True, help='Helper threads')
@click.option('--reference', required=False, type=str, default=None, show_default=True, help='Reference FASTA, required for CRAM')
@click.option('-t', '--tmpdir', required=True, type=click.Path(exists=True, resolve_path=True), help='Large preexisting tempdir')
@click.argument('in_bam_cram', nargs=1, required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True))
def thin_wrap(output, threads, reference, tmpdir, in_bam_cram):
    """Simple wrapper around pysam.collate to avoid requiring system installed samtools"""

    prefix = os.path.join(tmpdir, 'collate')
    # although this is intended to be used via a named pipe this will reduce the chance of it being filled
    collate_opts = ['--no-PG', '-r', '2500000', '-f', '--output-fmt', 'BAM', '-l', '1']
    collate_opts.extend(['-o', output])
    collate_opts.extend(['-@', str(threads)])
    if reference:
        collate_opts.extend(['--reference', reference])
    collate_opts.append(in_bam_cram)
    collate_opts.append(prefix)

    pysam.collate(*collate_opts)


