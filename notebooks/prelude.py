#
# notebooks/prelude.py
#

import os
from pathlib import Path

os.chdir(Path(__file__).parent.parent)

import sys

# fix alibuild environment
pathlist = [
    '/home/akubera/alice/root-venv/lib/python3.6/site-packages',
    '/home/akubera/development/physics/stumpy',
]
for p in map(Path, pathlist):
    if p.exists() and str(p) not in sys.path:
        sys.path.insert(1, str(p))
del p, pathlist

from itertools import chain, repeat, cycle, islice, product
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
from femtofitter.plotting import QuadPlot, normalize_subcanvases
from femtofitter.plotting import plot_outside, plot_outside_tdir
from femtofitter.pyfitter import (
    PyData3D,
    FitterGauss as PyFitterGauss
)

import ROOT
from ROOT import gROOT, gStyle, cppyy, TFile
from ROOT import TMinuit
from ROOT import TCanvas, TLegend, TLine, TText
from ROOT import TGraph, TGraphErrors
from ROOT import TH1, TH3, TF1
from ROOT import TH1D, TH1F, TH2D, TH2F, TH3D, TH3F
from ROOT import TProfile, TProfile2D
TH1.AddDirectory(False)

from ROOT import AliFemtoConfigObject
from ROOT import Data3D, Data1D
from ROOT import Fitter3DGaussLcms, Fitter3DLevy
from ROOT import Fitter1DGauss, Fitter1DGaussLin, Fitter1DGaussPolyBg
from ROOT import Fitter1DLevy, Fitter1DLevyPolyBg
from ROOT import FsiKFile, FsiGamov, FsiStatic
from ROOT import Mrc1DRatio, Mrc1DMatrix
from ROOT import Mrc3DRatio, Mrc3DHypercube
