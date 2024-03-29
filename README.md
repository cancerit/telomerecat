# Telomerecat (Telomere Computational Analysis Tool)

[![cancerit](https://circleci.com/gh/cancerit/telomerecat.svg?style=svg)](https://circleci.com/gh/cancerit/telomerecat)

Telomerecat is a tool for estimating the average telomere length (TL) for a paired end, whole genome sequencing (WGS) sample.

Telomerecat is adaptable, accurate and fast. The algorithm accounts for sequencing amplification artifacts, anneouploidy (common in cancer samples) and noise generated by WGS. For a high coverage WGS BAM file of around 100GB telomerecat can produce an estimate in ~1 hour.

## Docker container

Telomerecat is available as a Docker container on [Quay.io][quay-tags].

No "latest" image is defined, you need to specify the version you require, e.g.:

```bash
export VERSION_TEL=3.4.1 # update as appropriate
docker pull quay.io/wtsicgp/telomerecat:${VERSION_TEL}
```

## Singularity

The docker container is known to work with singularity, save the image locally via:

```bash
export VERSION_TEL=3.4.1 # update as appropriate
singularity pull docker://quay.io/wtsicgp/telomerecat:${VERSION_TEL}
```

## INSTALL

Installation is via `pip`.  Simply execute with the URL to a package release, e.g.:

```bash
pip3 install telomerecat
```

## Basic usage

Please see the command line help:

```bash
telomerecat --help
```

### Processes

When selecting the number of processes/threads the following should be considered:

* Single sample/input - 1, 2 or 4 recommended
* Multi sample/input - even values
  * parallel bam2telbam processes will be started with 2 cpus each (assuming >2 processes)

### Package Dependancies

`pip` will install the relevant dependancies, listed here for convenience:

* [parabam](https://github.com/cancerit/parabam)
* [numpy](https://numpy.org/)
* [pysam](https://www.scipy.org/)
* [pandas](https://pandas.pydata.org/)

## Development Dependencies

You will need virtualenv available on your system.

### Create a virtual python environement

```bash
cd $PROJECTROOT
python3 -m venv env
source env/bin/activate
python setup.py develop # so bin scripts can find module
```

## Cutting a release

1. Check version in `setup.py` has been updated
2. Check parabam version in `setup.py`/`Dockerfile`
3. Follow standard Hubflow release process (within cancerit)

CircleCI will handle docker image push to quay.io and package deployment to pypi.

<!-- Quay.io -->
[quay-tags]: https://quay.io/repository/wtsicgp/telomerecat?tab=tags
