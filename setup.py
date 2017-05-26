from setuptools import setup,find_packages
import os 
import sys
import shutil

#Instalise the __version__ variable
execfile('./telomerecat/_version.py')

setup(name='telomerecat',
	description='Telomere Computational Analysis Tool',
      version=__version__,
      author="JHR Farmery",
      license='GPL',
      author_email = 'jhrf2@cam.ac.uk',
      packages = ['telomerecat'],
      package_dir = {'telomerecat': 'telomerecat'},
      install_requires = ['parabam>=2.2',
                          'argparse',
                          'numpy',
                          'pysam==0.10.0',
                          'pandas'],
      include_package_data = True,
      scripts = ['./telomerecat/bin/telomerecat'],
      zip_safe = False
    )
