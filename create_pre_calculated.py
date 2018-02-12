#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:44:25 2018

@author: eejvt
"""

import numpy as np
from random import random
from random import uniform
import random
import numpy.random as nprnd
#import Jesuslib as jl
import matplotlib.pyplot as plt
from glob import glob
import poisson_error_function as pef



header='Fraction frozen, Temperatures, k, k_upper, k_lower'
for i in range(101,1000):
    data=pef.poisson_errors(i)
    np.savetxt('/nfs/see-fs-01_users/eejvt/PYTHON_CODE/INP_poisson_errors/pre_calculated/%i.csv'%i,data, delimiter=',',header=header)
