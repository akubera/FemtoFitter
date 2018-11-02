#!/usr/bin/env python
#
# fit-quick.py
#


import re
import sys
from datetime import datetime
import pandas as pd

from fit import run_fit


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
    with open(args.output, 'w') as f:
        df.to_json(df, f, indent=True)

    print(args.output)

    return 0


if __name__ == "__main__":
    main()
