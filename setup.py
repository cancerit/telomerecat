#!/usr/bin/env python3

from setuptools import setup

setup(
  name="telomerecat",
  version='3.4.1',
  description="Telomere Computational Analysis Tool",
  url='https://github.com/cancerit/telomerecat',
  author="JHR Farmery",
  license="GPL",
  python_requires='>= 3.6',
  author_email="cgphelp@sanger.ac.uk",
  packages=["telomerecat"],
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=2.3.2", "numpy", "pysam", "pandas", "easytimer"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  zip_safe=False,

)
