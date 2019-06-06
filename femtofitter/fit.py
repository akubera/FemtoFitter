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


@dataclass
class FitArgs:
    filename: str
    query: PathQuery
    fitter_classname: str
    fsi_classname: str
    fsi_args: Tuple[Any]
    fit_range: float
    minimum: float = 0.0
    fit_chi2: bool = False
    mrc_path: str = None
    subset: str = None

    @classmethod
    def BuildFrom(cls, datalist):
        for data in datalist:
            yield cls(data)

    def keys(self):
        from dataclasses import fields
        yield from (f.name for f in fields(self))

    def values(self):
        from dataclasses import astuple
        yield from astuple(self)

    def items(self):
        yield from vars(self).items()

    def __iter__(self):
        yield from self.keys()

    def __len__(self):
        from dataclasses import fields
        return len(fields(self))

    def __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError as e:
            raise KeyError() from e

    def run_fit(self):
        return run_fit(**self)


def run_fit(fitter_classname: str,
            filename: str,
            fsi_classname: str,
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

    if isinstance(fsi_classname, str):
        fsi_class = getattr(ROOT, fsi_classname)
    elif fsi_class is None:
        fsi_class = getattr(ROOT, 'FsiGamov')
    elif isinstance(fsi_classname, type):
        fsi_class = fsi_classname
    else:
        raise TypeError(f"Unexpected name of FSI calculator: {type(fsi_class)}")

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
                mrc_path = mrc_query.as_path()
            #     mrc_path = mrc_rootpath + "/mrc"

                mrc = mrc_tfile.Get(mrc_path)
                if mrc:
                    break
            else:
                print(f"Could not find MRC at {mrc_filename}:{mrc_rootpath}")
                return {}

            # if mrc_tfile is not tfile:
            #     mrc_tfile.Close()

            # print(f"Loaded MRC from file {mrc_filename} {mrc_rootpath}", )

        data = Data3D.FromDirectory(tdir, mrc, fit_range)
        # apply_momentum_resolution_correction(fitter.data, mrc)
    else:
        data = Data3D.FromDirectory(tdir, fit_range, minimum)

    if not subset:
        subset = None
    elif subset.startswith("sail"):
        data = data.sailor_subset()
        subset = 'sailor'
    elif subset.startswith("cow"):
        data = data.sailor_subset()
        subset = 'cowboy'
    elif subset.startswith("cone"):
        mf = query.magfield
        pr = query.pair

        same = ((mf == '--' and pr == 'pim') or
                (mf == '++' and pr == 'pip'))
        data = data.cone_subset(same)
        subset = 'cone'
    else:
        raise ValueError("Unexpected subset value: %r" % subset)

    fitter = fitter_class(data)
    fitter.fsi = fsi_class.new_shared_ptr(*fsi_args)
    # fitter.SetParamHintsFromDir(tdir)

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
    results['kT'] = float('%g' % mean(map(float, query.kt.split("_"))))
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
    output_path: str = None
    fsi_classname: str = 'FsiKFile'
    fsi_args: Tuple[any] = ('KFile2.root', )
    fitter_t: str = 'FitterGausOSL'
    mrc: bool = False
    mrc_only: bool = False
    fitrange: float = 0.11
    subset: Optional[str] = None
    ratio_min: float = 0.0
    chi2: bool = False
    limit: Optional[int] = None
    skip: Optional[int] = None
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
            except Exception:
                pass

        if not isinstance(fsi_args, tuple):
            fsi_args = (fsi_args, )

        return cls(tfile,
                   cli_args.output,
                   cli_args.fsi,
                   fsi_args,
                   cli_args.fitter,
                   cli_args.mrc_path,
                   cli_args.mrc_only,
                   cli_args.fitrange,
                   cli_args.subset,
                   cli_args.ratio_minimum,
                   cli_args.chi2,
                   cli_args.limit,
                   cli_args.skip,
                   cli_args.threads)


def pfit_all(args):
    return parallel_fit_all(args.tfile,
                            args.output_path,
                            (args.fsi_classname, args.fsi_args),
                            args.fitter_t,
                            args.mrc,
                            args.mrc_only,
                            args.fitrange,
                            args.subset,
                            args.ratio_min,
                            args.chi2,
                            args.limit,
                            args.skip,
                            args.threads)


def parallel_fit_all(tfile,
                     output_path=None,
                     fsi: Tuple[str, Tuple[Any]]=('FsiKFile', 'KFile2.root'),
                     fitter_t='FitterGausOSL',
                     mrc=False,
                     mrc_only=False,
                     fitrange=0.11,
                     subset=None,
                     ratio_min=0.0,
                     chi2=False,
                     limit=None,
                     skip=None,
                     threads=None):
    """
    """
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

    skip = skip or 0

    valid_analysis_keys = {'AnalysisQ3D', 'Q3DLCMS', 'Q3DPosQuad'}
    for analysis in valid_analysis_keys:
        search = f"{analysis}/{cfg}/{pair}/{cent}/{kt}/{mfield}"
        for path, _ in walk_matching(tfile, search):
            query = PathQuery.from_path(path)
            assert path == query.as_path()
            if skip > 0:
                skip -= 1
                continue
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
    # results = pool.starmap(run_fit_gauss, ((filename, p, fitrange) for p in paths[:4]))

    # work = chain(
    #     ((fitter_t, filename, *fsi, p, fitrange, ratio_min, chi2) for p in paths),
    #     ((fitter_t, filename, *fsi, p, fitrange, ratio_min, chi2, m) for p, m in mrc_paths),
    # )

    from dataclasses import replace
    fit_args = FitArgs(filename=filename,
                       query=None,
                       fitter_classname=fitter_t,
                       fsi_classname=fsi[0],
                       fsi_args=fsi[1],
                       fit_range=fitrange,
                       minimum=ratio_min,
                       fit_chi2=chi2,
                       mrc_path=None,
                       subset=subset)

    nomrc_fits = ()
    if not mrc_only:
        nomrc_fits = (replace(fit_args, query=p) for p in paths)
    mrc_fits = (replace(fit_args, query=p, mrc_path=m) for p, m in mrc_paths)

    work = chain(nomrc_fits, mrc_fits)

    with Pool(threads) as pool:
        results = pool.map(FitArgs.run_fit, work)

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
