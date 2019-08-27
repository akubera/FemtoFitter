#!/usr/bin/env python


import sys
import pandas as pd
import numpy as np
from stumpy import Histogram
from stumpy.utils import walk_matching
from ROOT import gSystem, TFile, TMinuit

gSystem.Load('build/libFemtoFitter.so')

from ROOT import Fitter3DGaussLcms


def series_from_result(result, **kwd):

    m = result.as_map()
    s = pd.Series(dict(m, **kwd))
    return s


try:
    fname = sys.argv[1]
except IndexError:
    fname = 'data-5882.root'

tfile = TFile.Open(fname)
if not tfile:
    exit(1)


cfg = 'cfg47AFEC7914D0B886'
pair = 'pim'
cent = '40_50'
kt = '0.8_1.0'
mfield = '--'

tdir = tfile.Get(f"Q3D/{cfg}/{pair}/{cent}/{kt}/{mfield}")
assert tdir

import asyncio
from concurrent.futures import ProcessPoolExecutor


def run_fit(path):
    from ROOT import TFile
    tfile = TFile.Open(fname)
    tdir = tfile.Get(path)

    _, cfg, pair, cent, kt, mfield = path.split('/')
    num, den = map(Histogram.BuildFromRootHist, map(tdir.Get, ("num", "den")))
    lim = .18
    fitter = Fitter3DGaussLcms.From(tdir, lim)

    results_pml = fitter.fit_chi2()
    ser = series_from_result(results_pml, cfg=cfg, pair=pair, cent=cent, kt=kt, field=mfield)

    return ser


executor = ProcessPoolExecutor()
loop = asyncio.get_event_loop()

path = f"Q3D/{cfg}/{pair}/*/*/{mfield}"
paths = (path for path, _ in walk_matching(tfile, path))
jobs = [loop.run_in_executor(executor, run_fit, path) for path in paths]

finished, _pending = loop.run_until_complete(asyncio.wait(jobs))
results = [r.result() for r in finished]

rdf = pd.DataFrame(results).sort_values(['cent', 'kt']).reset_index(drop=True)
print(rdf)
print(rdf[['cent', 'kt', 'Ro','Rs','Rl', 'lam']])

exit(0)


results = []
path = f"Q3D/{cfg}/{pair}/*/*/{mfield}"
for path, tdir in walk_matching(tfile, path):
    _, cfg, pair, cent, kt, mfield = path.split('/')
    print('  ', tdir)
    num, den = map(Histogram.BuildFromRootHist, map(tdir.Get, ("num", "den")))
    lim = .18
    fitter = Fitter3DGaussLcms.From(tdir, lim)

    results_pml = fitter.fit_pml()
    ser = series_from_result(results_pml, cfg=cfg, pair=pair, cent=cent, kt=kt, field=mfield)

    print(ser)
    results.append(ser)
    continue
    s = num.get_slice(*((-lim, lim), )*3)
    nps = tuple(slice(sl.start, sl.stop+1) for sl in s)




rdf = pd.DataFrame(results)
print(rdf)
print(rdf[['Ro','Rs','Rl', 'lam']])

exit(0)


# def build_fit_hist(histname):
#     h = tdir.Get(histname)
#     h


num, den = map(Histogram.BuildFromRootHist, map(tdir.Get, ("num", "den")))
lim = .18

# slices and numpy-slices
s = num.get_slice(*((-lim, lim), )*3)
nps = tuple(slice(sl.start, sl.stop+1) for sl in s)

n = np.array(num[nps].T.flat, dtype=int)
d = den[nps].T.flat

fitter = FitterGbaussOSL.From(tdir, lim)
numer = np.array(fitter.num_as_vec())
print(">>", numer.shape)

# cfg = 'cfg47AFEC7914D0B886'
# pair = 'pim'
# cent = '40_50'
# kt = '0.8_1.0'
# mfield = '--'

# exit(0)
results_pml = fitter.fit_pml()
ser = series_from_result(results_pml, cfg=cfg, pair=pair, cent=cent, kt=kt, field=mfield)
print(ser)
input()
print("#")
print("# BEGIN CHI2")
print("#")
results_chi2 = fitter.fit_chi2()


print("PML:")
results_pml.print()
print("Chi2:")
results_chi2.print()

# print(fitter)

qspace = num.meshgrid_centerbin

qout = qspace[0][nps].T.flat
qside = qspace[1][nps].T.flat
qlong = qspace[2][nps].T.flat

# print(qout.shape)
# print(qspace[0].shape)
print('------')
#/
#Fit-Result: Ro=5.01354 Rs=5.31392 Rl=5.03968
#Fit-Result: Ro=5.04053 Rs=5.29246 Rl=5.01207
#
#Fit-Result: Ro=5.01354 Rs=5.31392 Rl=5.03968
#Fit-Result: Ro=5.04027 Rs=5.29291 Rl=5.01106
#------

def calc_chi2(n, d, c):
    r = np.divide(n, d, where=d!=0.0, out=np.zeros_like(c))
    v = np.divide((1.0 + r) * r * r, n, out=np.zeros_like(c))

    return (r - c) ** 2 / v


# print(np.any(d == n))

# for i in range(34, 49):
#   print(i, fitter.data.num[i], n[i] / d[i],
#     n[i] > d[i]
#   )
#           # fitter.data.qspace[0][i], qout[i],
#           # fitter.data.qspace[1][i], qside[i],
#           # fitter.data.qspace[2][i], qlong[i])

# print()
# exit(0)

# # minuit = TMinuit()
