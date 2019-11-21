from setuptools import setup, find_packages
import os
import sys
import shutil

from telomerecat.version import __version__

setup(
  name="telomerecat",
  description="Telomere Computational Analysis Tool",
  version=__version__,
  author="JHR Farmery",
  license="GPL",
  author_email="cgpit@sanger.ac.uk",
  packages=["telomerecat"],
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=2.2", "argparse", "numpy", "pysam", "pandas"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  zip_safe=False,
)
