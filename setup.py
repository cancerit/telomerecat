#!/usr/bin/env python3

from setuptools import setup

setup(
  name="telomerecat",
  description="Telomere Computational Analysis Tool",
  version='3.4.1',
  author="JHR Farmery",
  license="GPL",
  python_requires='>= 3.6',
  author_email="cgphelp@sanger.ac.uk",
  packages=["telomerecat"],
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=2.3.0", "numpy", "pysam", "pandas"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  zip_safe=False,

)
