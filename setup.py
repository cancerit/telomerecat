from setuptools import setup, find_packages

from telomerecat.version import __version__

setup(
  name="telomerecat",
  description="Telomere Computational Analysis Tool",
  version=__version__,
  author="JHR Farmery",
  license="GPL",
  python_requires': '>= 3.6',
  author_email="cgpit@sanger.ac.uk",
  packages=["telomerecat"],
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=2.2.5", "argparse", "numpy", "pysam", "pandas"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  zip_safe=False,
)
