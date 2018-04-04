#!/usr/bin/env python

import sys

sys.dont_write_bytecode = True

import datetime
import urllib
from astropy.io import fits
import numpy as np
import json
import pickle
import time
import urllib
import glob
import pandas as pd
# import feather
import os.path
from sunpy.net import jsoc, attrs, vso
import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import shutil

HMI_tmp_path = '/d1/flares/tmp'

for count in range(0,6000):

    try:
        os.rmdir(HMI_tmp_path + '/%05d' % count)
        print('Removing empty folder: ' + HMI_tmp_path + '/%05d' % count)
    except: 
        empty = 1

