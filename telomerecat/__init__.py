#__init__.py
import core

from telomerecat.csv2length import Csv2Length
from telomerecat.bam2telbam import Bam2Telbam
from telomerecat.telbam2length import Telbam2Length
from telomerecat.bam2length import Bam2Length

import bam2telbam
import bam2length
import telbam2length
import csv2length
import error_estimator
import readmodel

from _version import __version__