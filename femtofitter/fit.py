#!/usr/bin/env python
#
# fit.py
#
"""
Fitting routines
"""

import sys
from itertools import chain
from multiprocessing import Pool, Process
from statistics import mean
import json

import pandas as pd

from femtofitter import PathQuery


def find_and_fit(filename: str, query: PathQuery, fit_range: float, mrc_path=None):
    query = PathQuery.From(query)

    path = query.as_path()

    results = run_fit(filename, path, fit_range, mrc_path)
    results.update(query.as_dict())
    results['kt'] = mean(map(float, query['kt'].split("_")))

    return results


# def run_fit_gauss(*args, **kwargs):
#     from ROOT import FitterGaussOSL
#     return run_fit(FitterGaussOSL, *args, **kwargs)

def run_fit_gauss(*args, **kwargs):
    from ROOT import Data3D
    from ROOT import FitterGaussOSL
    return run_fit(FitterGaussOSL, *args, **kwargs)


def run_fit_levy(*args, **kwargs):
    from ROOT import Data3D
    from ROOT import FitterLevy
    return run_fit(FitterLevy, *args, **kwargs)


def run_fit_gauss_full(*args, **kwargs):
    from ROOT import Data3D
    from ROOT import FitterGaussFull
    return run_fit(FitterGaussFull, *args, **kwargs)

def run_fitter(fitter, *args, **kwargs):
    import ROOT
    from ROOT import Data3D
    Fitter = getattr(ROOT, fitter)

    return run_fit(Fitter, *args, **kwargs)

def run_fit(fitter_class,
            filename: str,
            query: PathQuery,
            fit_range: float,
            fit_chi2: bool = False,
            mrc_path: str = None):

    from ROOT import TFile
    from ROOT import Data3D
    tfile = TFile.Open(filename)
    assert tfile

    path = query.as_path()
    tdir = tfile.Get(path)
    assert tdir

    fitter = fitter_class.From(tfile, path, fit_range)

    if mrc_path is not None:
        from ROOT import apply_momentum_resolution_correction
        mrc_file = TFile.Open("Data-MRC-1987.root")
        mrc = mrc_file.Get(mrc_path)
        mrc_file.Close()
        # print("MRC:",mrc)
        if mrc == None:
            print("missing mrc path %r" % mrc_path)
            return {}
        apply_momentum_resolution_correction(fitter.data, mrc)

    fit_results = fitter.fit_chi2() if fit_chi2 else fitter.fit_pml()
    if not fit_results:
        print(f"Could not fit: {query.as_path()}")
        return {}

    results = dict(fit_results.as_map())
    results.update(query.as_dict())

    results['fit_range'] = fit_range

    results['kT'] = mean(map(float, query.kt.split("_")))
    results['chi2'] = fitter.resid_chi2(fit_results)
    results['ndof'] = fitter.size()
    results['rchi2'] = results['chi2'] / results['ndof']
    results['mrc'] = mrc_path

    return results


def parallel_fit_all(tfile,
                     ofilename=None,
                     fitter_t='FitterGausOSL',
                     mrc=False,
                     fitrange=0.21,
                     chi2=False):
    """
    """

    from stumpy.utils import walk_matching
    from datetime import datetime
    from pathlib import Path

    import ROOT
    try:
        getattr(ROOT, fitter_t)
    except AttributeError:
        print(f"Could not load fitter {fitter_t!r}", file=sys.stderr)
        return

    filename = Path(str(tfile.GetName()))

    cfg = 'cfg*'
    pair = cent = kt = mfield = '*'
    search = f"AnalysisQ3D/{cfg}/{pair}/{cent}/{kt}/{mfield}"
    if mrc:
        if isinstance(mrc, PathQuery):
            mrc_path = mrc.as_path()
        else:
            mrc_cfg = mrc if isinstance(mrc, str) else "cfg5348379EA4DD77C6"
            mrc_path = "AnalysisTrueQ3D/%s/{pair}/00_90/{kt}/{magfield}/mrc" % mrc_cfg
    else:
        mrc_path = None

    paths = []
    mrc_paths = []
    for path, _ in walk_matching(tfile, search):
        query = PathQuery.from_path(path)
        assert path == query.as_path()
        paths.append(query)
        if mrc_path:
            mrc_paths.append((query, mrc_path.format(**query.as_dict())))

    configuration_information = get_configuration_json(tfile, paths)

    filename = str(filename.absolute())
    pool = Pool()
    # results = pool.starmap(run_fit_gauss, ((filename, p, fitrange) for p in paths[:4]))

    work = chain(
        ((fitter_t, filename, p, fitrange, chi2) for p in paths),
        ((fitter_t, filename, p, fitrange, chi2, m) for p, m in mrc_paths),
    )

    results = pool.starmap(run_fitfitter, work)
    #results += pool.starmap(run_fit_gauss, ])
    #results += pool.starmap(run_fit_levy, [(filename, p, fitrange) for p in paths[:1]])
    df = pd.DataFrame(results)
    output_data = {
        'filename': filename,
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
    for c in {"%s/%s" % (q.analysis, q.cfg) for q in queries}:
        path = f"{c}/config"
        config = tfile.Get(path)
        # if not config:
        if not isinstance(config, AliFemtoConfigObject):
            print("Missing AliFemtoConfigObject 'config' in %r" % path)
            continue

        result[c] = json.loads(config.as_JSON_string())

    return result


if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except IndexError:
        filename = 'data.root'

    from ROOT import gSystem, TFile
    assert 0 == gSystem.Load('build/libFemtoFitter.so')

    data = TFile.Open(filename)
    parallel_fit_all(data, "fitres-%s.json" % filename.rpartition('.')[0])
