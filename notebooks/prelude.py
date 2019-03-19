#
# notebooks/prelude.py
#

import os
from pathlib import Path

os.chdir(Path(__file__).parent.parent)

import sys
from itertools import chain, repeat, cycle, islice
from functools import partial, reduce
from copy import copy

from scipy.interpolate import interp1d, CubicSpline
from statistics import mean

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import yaml
import feather

import stumpy
from stumpy.utils import iter_tobject, walk_matching
from stumpy.rhist import normalize_hist, average_hist_contents

import lmfit
from femtofitter import FitResults, PathQuery
from femtofitter.plotting import QuadPlot
from femtofitter.pyfitter import (
    PyData3D,
    FitterGauss as PyFitterGauss
)

import ROOT
from ROOT import gROOT, cppyy, TFile
from ROOT import TCanvas, TLegend, TLine, TText
from ROOT import TH1, TH3, TF1
