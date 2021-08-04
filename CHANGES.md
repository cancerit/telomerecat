# CHANGES

## 4.0.1

* Add minimal compression to `pysam_collate` stability issue, suspect filling pipe buffer
* Add threads to hts writer of telbam
* Remove requirements file as controlled by `setup.py`
* Migrate to circleci (automate pypi deploy)

## 4.0.0

* `telomerecat` no longer depends on parabam for BAM processing in bam2telbam (due to first point)
  * Multiprocessing `fork` is slowly being dropped:
    * [see here](https://docs.python.org/3/library/multiprocessing.html#multiprocessing.get_context)
    * At present all subcommands are working under py 3.8
  * `samtools collate` can rapidly bring pairs together, rule migrated to dedicated code
  * Slightly faster
* CRAM inputs now supported, output alignments remain BAM

## 3.4.0

* Fixed float being passed to randint function.
* Now installs cython as CI was breaking.

## 3.3.0

* Code are migrated to python3, and does not support python2 anymore.
* added command line option `--temp_dir` to allow user to specify a folder for intermediate files.
* added command line option `--outbam_dir` to allow user to specify a folder to save telBAMs to.
* added Dockerfile
* added hidden command line option `--seed_randomness` to produce stable output, but it runs very slow and in risk of voilating telomerecat's original design. ***Not*** recommand to use other than validating outputs between code changes.

2017/12/8 - Telomerecat v3.2 -jhrf

- Fixed a bug in the command line interface where the specificed number
  of runs of the simulation was not passed to the simulator. The default value
  would always be used instead.
- It is now default to NOT use F2a correction. The user must now specify the use
  of F2a correction with the `-e` parameter. The `-d` parameter (disable-correction)
  has been deprecated
- General code format corrections across all files.


2017/5/26 - Telomerecat v3.1.2 - jhrf

- Added a check in bam2telbam rule that ensures the read.seq is loaded into
  the read and does not return None. This was causing a hang in PCAWG samples
- Tidied up the code surrounding mismatching loci following code review.
  Process should be easier to follow now
- Minor bug fixes and restructing in telbam2length
- Changed the requirements such that a fixed version of pysam is required
  this version is proven stable and more recent versions have caused crashes

2017/1/24 - Telomerecat v3.1.1 - jhrf

- Removed debug print message from telbam2length

2017/1/17 - Telomerecat v3.1 - jhrf

- Added the ability to specify a read trim via the command line
- Made the telbam2length estimation protocol more responsive to differing
  read lengths (i.e removed magic numbers)
- Fixed command line parameterisation bugs in csv2length and telbam2length
