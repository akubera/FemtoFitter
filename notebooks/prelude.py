#
# notebooks/prelude.py
#

import os
from pathlib import Path

os.chdir(Path(__file__).parent.parent)

import sys

from scipy.interpolate import interp1d, CubicSpline
from statistics import mean

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import yaml
import feather

import stumpy
import lmfit
from femtofitter.pyfitter import Data3D as PyDat3D, FitterGauss as PyFitterGauss

import ROOT
from ROOT import gROOT, cppyy, TFile, TCanvas, TH1