#!/usr/bin/env python
#
# fit-quick.py
#


import re
import sys
from datetime import datetime
import pandas as pd

MRC_FILE = "mrc.root"


def arg_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('filename',
                        default="~/Physics/data/1234813315/AnalysisResults.root",
                        help='input rootfile')
    parser.add_argument("-b", "--bounds",
                        type=float,
                        default=0.21)
    parser.add_argument("-o", "--output",
                        default=datetime.now().strftime(r"FitResults/fit-quick-%Y%m%d%H%M%S.json"))
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = arg_parser().parse_args(argv)

    from ROOT import gSystem, TFile
    assert 0 == gSystem.Load('build/libFemtoFitter.so')

    cfg = 'cfgB64295A7EE5BD10E'
    pair = 'pim'
    cent = '00_10'
    kt = '0.200_0.300'
    mfield = '--'

    path = f"{cfg}/{pair}/{cent}/{kt}/{mfield}"
    path = 'cfgB64295A7EE5BD10E/pim/0_10/0.200_0.300/--'
    print(path)

    results = run_fit(args.filename, path, args.bounds)

    df = pd.DataFrame(results)
    df.to_json(args.output)

    print(args.output)

    return 0


def run_fit(filename: str, path: str, fit_range: float, mrc_path: str=None):
    from ROOT import TFile
    tfile = TFile.Open(filename)

    tdir = tfile.Get(path)
    assert tdir

    from ROOT import apply_momentum_resolution_correction
    from ROOT import FitterGaussOSL, FitterGaussFull
    from ROOT import FitterLevy

    FITTER = FitterGauss

    fitter = FITTER.From(tfile, path, fit_range)

    fit_results = fitter.fit()
    results = dict(fit_results.as_map())

    results['fit_range'] = fit_range

    results['chi2'] = fitter.resid_chi2(fit_results)
    results['ndof'] = fitter.ndof
    # results['chi2_per_ndof'] = results['chi2'] / results['ndof']
    # results['pml_per_ndof'] = fitter.reduced_pml(result.params)
    return [results]


if __name__ == "__main__":
    main()
