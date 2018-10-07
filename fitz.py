#!/usr/bin/env python
#
# fitz.py
#

import re
import sys
import asyncio
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

import lmfit
from lmfit import Minimizer, Parameters
import numpy as np
import pandas as pd

from stumpy import Histogram
from stumpy.utils import walk_matching

from fitterz import FitterGauss, FitterGauss4, FitterLevy

FITTER = FitterGauss
FITTER = FitterLevy


def arg_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('filename',
                        default="~/Physics/data/1234813315/AnalysisResults.root",
                        help='input rootfile')
    parser.add_argument("-o", "--output",
                        default=datetime.now().strftime(r"FitResult-%Y%m%d%H%M%S.json"))
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = arg_parser().parse_args(argv)

    loop = asyncio.get_event_loop()
    results = loop.run_until_complete(start(args.filename))

    df = pd.DataFrame(results)
    df.to_json(args.output)

    print(args.output)

    return 0


async def start(filename,
                count=-1,
                fitters=('Fitter', ),
                *,
                fit_range=0.22,
                kt='*',
                centrality='*',
                pair_type='*',
                cfg='cfg*'):

    from ROOT import gROOT, TFile
    gROOT.SetBatch()

    tfile = TFile.Open(filename)
    if not tfile:
        raise FileNotFoundError(filename)

    loop = asyncio.get_event_loop()
    executor = ProcessPoolExecutor()

    jobs = []

    path = f"Q3D/{cfg}/{pair_type}/{centrality}/{kt}"

    if isinstance(fit_range, float):
        fit_range = (fit_range, )

    for path, container in walk_matching(tfile, path):
        for fitrange in fit_range:
            for field in ('--', '++'):
                job = loop.run_in_executor(executor, run_fit, filename, path + "/" + field, fitrange)
                jobs.append(job)
                count -= 1

                if not count:
                    break

    job_count = len(jobs)
    results = []
    while jobs:
        finished, jobs = await asyncio.wait(jobs, timeout=20)
        results.extend(finished)
        print('%d / %d' % (len(results), job_count))

    fit_results = []
    for r in (x.result() for x in results):
        if not r:
            print("bad r")
        elif isinstance(r, list):
            fit_results.extend(r)
        else:
            fit_results.append(r)

    return fit_results


def run_fit(filename, path, fit_range, fitter_class=FITTER):
    from ROOT import TFile
    tfile = TFile.Open(filename)
    data = tfile.Get(path)
    if not data:
        print("Error reading container %r from %r" % (filename, path))
        return

    r = re.compile(r'Q3D/(?P<CFG>cfg\w+)/(?P<pair_type>\w+)/'
                   r'(?P<centrality>[^/]+)/(?P<kt>[^/]+)/(?P<field>[^/]+)')
    try:
        path_info = r.match(path).groupdict()
    except AttributeError:
        return

    Hist = Histogram.BuildFromRootHist

    num, den, qinv = map(lambda name: Hist(data.Get(name)), ('num', 'den', 'qinv'))

    fitter = fitter_class(num, den, qinv, fit_range)

    q3d_params = fitter.default_parameters()
    # mini = lmfit.Minimizer(fitter.resid_loglike, q3d_params, reduce_fcn=np.sum, nan_policy='omit')
    # mini = lmfit.Minimizer(fitter.resid, q3d_params, reduce_fcn=fitter.chi2, nan_policy='omit')

    # prev_result = mini.minimize(method='nelder')
    # prev_result = mini.minimize(method='leastsq')

    # result = mini.minimize(method='leastsq', params=prev_result.params)
    
    perlim_result = lmfit.minimize(fitter.resid_loglike, q3d_params, reduce_fcn=np.sum, method='nelder')
    result = lmfit.minimize(fitter.resid, perlim_result.params, reduce_fcn=fitter.chi2, method='leastsq')
    lam_ratio = perlim_result.params['lam'].value / result.params['lam'].value
    if lam_ratio < .9:
        print('lam %f %s -> %s' % (lam_ratio, perlim_result.params['lam'], result.params['lam']))
        next_result = lmfit.minimize(fitter.resid,
                                     result.params,
                                     reduce_fcn=fitter.chi2,
                                     method='leastsq')
        if np.isclose(next_result.params['lam'].value, result.params['lam'].value, 1e-3):
            print('not much change')
            for k in result.params:
                result.params[k].value = perlim_result.params[k].value 
        else:
            print("next... ", next_result.params)
            result = next_result

    iterations = 1
    # minimizer
    # result = mini.leastsq(prev_result.params)

    # chi2 = fitter.resid_chi2(result.params).sum()
    # rchi2 = chi2 / ndof
    #print(path, chi2, rchi2)

    # while abs(result.params['Rl'] - prev_result.params['Rl']) > 1e-4:
    # rchi2 = 1000
    # while rchi2 > 5:
    #     iterations += 1

    #     prev_result = result
    #     result = mini.minimize(method=method, params=prev_result.params)
    #     rchi2 = fitter.reduced_chi2(result.params)

    #     keys = ("Ro", "Rs", "Rl")
    #     a = np.array([[res[key].value for key in keys] for res in (prev_result.params, result.params)])
    #     # convergence! all parameters are "close" (i.e. within 1e-8)
    #     if np.allclose(*a):
    #         break

    #     if iterations > 50:
    #         break

    results = path_info

    for k in result.params:
        p = result.params[k]
        results[k] = p.value
        results[k + "_err"] = p.stderr

    results['fit_range'] = fit_range
    results['iterations'] = iterations

    results['chi2'] = fitter.resid_chi2(result.params).sum()
    results['ndof'] = fitter.ndof
    results['chi2_per_ndof'] = results['chi2'] / results['ndof']
    results['pml_per_ndof'] = fitter.reduced_pml(result.params)

    return results


if __name__ == "__main__":
    main()
