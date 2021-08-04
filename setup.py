#!/usr/bin/env python3

from setuptools import setup

setup(
  name="telomerecat",
  version='4.0.1',
  description="Telomere Computational Analysis Tool",
  long_description=open('README.md').read(),
  long_description_content_type='text/markdown',
  url='https://github.com/cancerit/telomerecat',
  author="JHR Farmery",
  license="GPL",
  python_requires='>= 3.7',
  author_email="cgphelp@sanger.ac.uk",
  packages=["telomerecat"],
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=3.0.1", "numpy", "pysam", "pandas", "click"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  entry_points={'console_scripts': ['pysam_collate=telomerecat.pysam_collate:thin_wrap'],},
  zip_safe=False,
)
