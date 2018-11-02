#!/usr/bin/env python
#
# fit.py
#
"""
Fitting routines
"""

from multiprocessing import Pool, Process
from statistics import mean
from collections import namedtuple
from typing import Dict, NamedTuple
import json

import pandas as pd

# PathQuery = namedtuple("PathQuery", ('cfg', 'pair', 'cent', 'kt', 'mfield'))

class PathQuery(NamedTuple):

    cfg: str
    pair: str
    cent: str
    kt: str
    magfield: str

    def as_path(self):
        return '/'.join(self)

    def as_dict(self):
        return {key: getattr(self, key)
                for key in ('cfg', 'pair', 'cent', 'kt', 'magfield')}

    @classmethod
    def from_path(cls, path):
        if isinstance(path, str):
            path = path.split("/")
        return cls(*path)

# def run_fit1(filename: str, cfg: str, pair: str, cent: str, kt: str, mfield: str, fit_range: float, mrc_path: str=None):
def find_and_fit(filename: str, query: PathQuery, fit_range: float):

    path = query.as_path()

    results = run_fit(filename, path, fit_range, mrc_path)
    results.update(query.as_dict())
    results['kt'] = mean(map(float, query['kt'].split("_")))

    return results


# def run_fit_gauss(*args, **kwargs):
#     from ROOT import FitterGaussOSL
#     return run_fit(FitterGaussOSL, *args, **kwargs)

def run_fit_gauss(filename, path, fit_range):
    from ROOT import FitterGaussOSL
    return run_fit(FitterGaussOSL, filename, path, fit_range)

def run_fit_levy(*args, **kwargs):
    from ROOT import FitterLevy
    return run_fit(FitterLevy, *args, **kwargs)

def run_fit_gauss_full(*args, **kwargs):
    from ROOT import FitterGaussFull
    return run_fit(FitterGaussFull, *args, **kwargs)


def run_fit(fitter_class,
            filename: str,
            query: PathQuery,
            fit_range: float,
            mrc_path: str=None):

    import ROOT
    # TFile, FITTER = bfrom ROOT import TFile
    tfile = TFile.Open(filename)
    assert tfile

    path = query.as_path()
    tdir = tfile.Get(path)
    assert tdir

    from ROOT import apply_momentum_resolution_correction

    fitter = fitter_class.From(tfile, path, fit_range)

    fit_results = fitter.fit()
    results = dict(fit_results.as_map())
    results.update(query.as_dict())

    results['fit_range'] = fit_range

    results['chi2'] = fitter.resid_chi2(fit_results)
    results['ndof'] = fitter.size()
    results['rchi2'] = results['chi2'] / results['ndof']

    return results


def parallel_fit_all(tfile, ofilename=None):
    """
    """

    from stumpy.utils import walk_matching
    from datetime import datetime
    from pathlib import Path

    filename = Path(str(tfile.GetName()))

    cfg = 'cfg*'
    pair = cent = kt = mfield = '*'
    search = f"{cfg}/{pair}/{cent}/{kt}/{mfield}"

    paths = []
    for path, _ in walk_matching(tfile, search):
        query = PathQuery.from_path(path)
        assert path == query.as_path()
        paths.append(query)

    configuration_information = get_configuration_json(tfile, paths)

    fitrange = 0.21
    pool = Pool()
    results = pool.starmap(run_fit_gauss, ((str(filename), p, fitrange) for p in paths))
    df = pd.DataFrame(results)
    output_data = {
        'filename': str(filename.absolute()),
        'timestamp': datetime.now().isoformat(timespec='milliseconds'),
        'df': df.to_dict(orient='records'),
        'config': configuration_information,
    }

    if ofilename:
        with Path(ofilename).open('w') as outfile:
            json.dump(output_data, outfile, indent=True)
    else:
        from pprint import pprint
        pprint(output_data)


def get_configuration_json(tfile, queries):
    from ROOT import AliFemtoConfigObject

    result = {}
    for c in {q.cfg for q in queries}:
        path = f"{c}/config"
        config = tfile.Get(path)
        # if not config:
        if not isinstance(config, AliFemtoConfigObject):
            print("Missing AliFemtoConfigObject 'config' in %r" % path)
            continue

        result[c] = json.loads(config.as_JSON_string())

    return result


if __name__ == "__main__":
    import sys
    try:
        filename = sys.argv[1]
    except IndexError:
        filename = 'data.root'

    from ROOT import gSystem, TFile
    assert 0 == gSystem.Load('build/libFemtoFitter.so')

    data = TFile.Open(filename)
    parallel_fit_all(data)
