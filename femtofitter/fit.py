#!/usr/bin/env python
#
# fit.py
#
"""
Fitting routines
"""

import sys
from copy import copy
from itertools import chain
from dataclasses import dataclass
from multiprocessing import Pool, Process
from statistics import mean
import json

import pandas as pd

from femtofitter import PathQuery

from typing import Any, Tuple, Optional


def run_fitter(fitter, *args, **kwargs):
    import ROOT
    from ROOT import Data3D
    Fitter = getattr(ROOT, fitter)

    return run_fit(Fitter, *args, **kwargs)


def run_fit(fitter_classname: str,
            filename: str,
            fsi_class: str,
            fsi_args: Tuple[Any],
            query: PathQuery,
            fit_range: float,
            minimum: float = 0.0,
            fit_chi2: bool = False,
            mrc_path: str = None,
            subset: str = None,
            ) -> dict:
    """
    Finds data at (filename, query), and fits using the
    remaining parameters
    """

    import ROOT
    from ROOT import TFile
    from ROOT import Data3D

    # load fitter-classname
    if isinstance(fitter_classname, str):
        fitter_class = getattr(ROOT, fitter_classname)
    elif isinstance(fitter_classname, type):
        fitter_class = fitter_classname
    else:
        raise TypeError(f"Unexpected type {type(fitter_classname)}")

    if isinstance(fsi_class, str):
        fsi_class = getattr(ROOT, fsi_class)
    elif fsi_class is None:
        fsi_class = getattr(ROOT, 'FsiGamov')
    elif not isinstance(fsi_class, type):
        raise TypeError(f"Unexpected type of FSI calculator: {type(fsi_class)}")

    tfile = TFile.Open(filename)
    if not tfile:
        raise FileNotFoundError()

    path = query.as_path()
    tdir = tfile.Get(path)
    assert tdir

    if mrc_path is not None:
        #from ROOT import apply_momentum_resolution_correction
        if mrc_path is True:
            mrc_filename = filename
        elif mrc_path is ...:
            mrc_filename = "Data-MRC-1987.root"
        else:
            #mrc_filename = 'Data-sbin.root'
            mrc_filename, _, cfg = mrc_path.partition(':')
            mrc_query = query(cent='00_90')
            mrc_analyses = [
                'AnalysisTrueQ3D',
                'TrueQ3D',
            ]
            if cfg:
                mrc_query.cfg = cfg

            if mrc_filename and mrc_filename != filename:
                mrc_tfile = TFile.Open(mrc_filename)
            else:
                mrc_tfile = tfile

            for analysis in mrc_analyses:
                mrc_query.analysis = analysis
                mrc_rootpath = mrc_query.as_path()
                mrc_path = mrc_rootpath + "/mrc"

                mrc = mrc_tfile.Get(mrc_path)
                if mrc:
                    break
            else:
                print(f"Could not find MRC at {mrc_filename}:{mrc_rootpath}")
                return {}

            if mrc_tfile is not tfile:
                mrc_tfile.Close()

            # print(f"Loaded MRC from file {mrc_filename} {mrc_rootpath}", )

        data = Data3D.FromDirectory(tdir, mrc, fit_range, minimum)
        # apply_momentum_resolution_correction(fitter.data, mrc)
    else:
        data = Data3D.FromDirectory(tdir, fit_range, minimum)

    if subset == 'sailor':
        data = data.cowboy_subset()
    elif subset == 'cowboy':
        data = data.sailor_subset()

    fitter = fitter_class(data)
    fitter.fsi = fsi_class.new_shared_ptr(*fsi_args)
    fitter.SetParamHintsFromDir(tdir)

    results = {
        'fsi': str(fitter.fsi.ClassName()),
        'mrc': mrc_path,
    }
    fit_results = fitter.fit_chi2() if fit_chi2 else fitter.fit_pml()
    if not fit_results:
        print(f"Could not fit: {query.as_path()}")
        return {}

    results.update(fit_results.as_map())
    results.update(query.as_dict())

    results['fit_range'] = fit_range
    results['subset'] = subset

    results['cent'] = query.cent
    results['kT'] = mean(map(float, query.kt.split("_")))
    results['chi2'] = fitter.resid_chi2(fit_results)
    results['ndof'] = fitter.degrees_of_freedom()
    results['rchi2'] = results['chi2'] / results['ndof']
    #results['mrc'] = mrc_path
    results['gamma'] = fitter.data.gamma
    #results['fsi'] = f'{fsi_class.__name__}{fsi_args}'

    return results


@dataclass
class ParallelFitArgs:
    tfile: str
    fsi: str
    fsi_args: Tuple[any] = ()
    output_path: str = None
    fitter_t: str ='FitterGausOSL'
    mrc: bool = False
    fitrange: float = 0.11
    ratio_min: float = 0.0
    chi2: bool = False
    limit: Optional[int] = None
    threads: Optional[int] = None

    @classmethod
    def FromCli(cls, tfile, cli_args):
        """
        Build from command line arguments
        """
        fsi_args = cli_args.fsi_args
        if isinstance(fsi_args, str):
            # try interpreting as python literal
            import ast
            try:
                fsi_args = ast.literal_eval(fsi_args)
            except :
                pass

        if not isinstance(fsi_args, tuple):
            fsi_args = (fsi_args, )

        return cls(tfile,
                   cli_args.fsi,
                   fsi_args,
                   cli_args.output,
                   cli_args.fitter,
                   cli_args.mrc_path,
                   cli_args.fitrange,
                   cli_args.ratio_minimum,
                   cli_args.chi2,
                   cli_args.limit,
                   cli_args.threads)


def pfit_all(args):
    return parallel_fit_all(args.tfile,
                            (args.fsi, args.fsi_args),
                            args.output_path,
                            args.fitter_t,
                            args.mrc,
                            args.fitrange,
                            args.ratio_min,
                            args.chi2,
                            args.limit,
                            args.threads)


def parallel_fit_all(tfile,
                     fsi: Tuple[str, Tuple[Any]],
                     output_path=None,
                     fitter_t='FitterGausOSL',
                     mrc=False,
                     fitrange=0.21,
                     ratio_min=0.0,
                     chi2=False,
                     limit=None,
                     threads=None):
    """
    """
    print("MINIMUM", ratio_min)

    from stumpy.utils import walk_matching
    from datetime import datetime
    from pathlib import Path

    import ROOT
    try:
        fitter_t = getattr(ROOT, fitter_t)
    except (AttributeError, TypeError):
        print(f"Could not load fitter {fitter_t!r}", file=sys.stderr)
        return

    filename = Path(str(tfile.GetName()))

    if mrc:
        if isinstance(mrc, PathQuery):
            mrc_path = mrc.as_path()
        else:
            #mrc_cfg = mrc if isinstance(mrc, str) else "cfgBDC0F09B1F286D46"
            #mrc_path = "AnalysisTrueQ3D/%s/{pair}/00_90/{kt}/{magfield}/mrc" % mrc_cfg
            mrc_path = mrc
    else:
        mrc_path = None

    cfg = 'cfg*'
    pair = cent = kt = mfield = '*'

    paths = []
    mrc_paths = []
    limit_reached = False

    valid_analysis_keys = {'AnalysisQ3D', 'Q3DLCMS'}
    for analysis in valid_analysis_keys:
        search = f"{analysis}/{cfg}/{pair}/{cent}/{kt}/{mfield}"
        for path, _ in walk_matching(tfile, search):
            query = PathQuery.from_path(path)
            assert path == query.as_path()
            paths.append(query)
            if mrc_path:
                mrc_paths.append((query, mrc_path))
            if limit is not None:
                if len(paths) >= limit:
                    limit_reached = True
                    break
        if limit_reached:
            break

    configuration_information = get_configuration_json(tfile, paths)

    filename = str(filename.absolute())
    pool = Pool(threads)
    # results = pool.starmap(run_fit_gauss, ((filename, p, fitrange) for p in paths[:4]))

    work = chain(
        ((fitter_t, filename, *fsi, p, fitrange, ratio_min, chi2) for p in paths),
        ((fitter_t, filename, *fsi, p, fitrange, ratio_min, chi2, m) for p, m in mrc_paths),
    )

    results = pool.starmap(run_fit, work)

    df = pd.DataFrame(results)
    output_data = {
        'filename': filename,
        'timestamp': datetime.now().isoformat(timespec='milliseconds'),
        'command': sys.argv,
        'df': df.to_dict(orient='records'),
        'config': configuration_information,
    }

    if output_path:
        with Path(output_path).open('w') as outfile:
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
