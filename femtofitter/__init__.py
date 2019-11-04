#
# femtofitter/__init__.py
#

import json
# from typing import NamedTuple
from dataclasses import dataclass
import pandas as pd


def __getattr__(key):
    import ROOT
    return getattr(ROOT, key)


def load_femtofitter_lib(exit_on_err=1):
    """
    Attempt to load libFemtoFitter.so
    """
    from sys import stderr
    from pathlib import Path
    from ROOT import gInterpreter

    fitlibrary = Path(__file__).parent.parent / 'build' / 'libFemtoFitter.so'
    fitlibrary = str(fitlibrary.absolute())
    retcode = gInterpreter.Load(fitlibrary)
    if retcode < 0:
        print("Could not load libFemtoFitter.so", fitlibrary, file=stderr)
        if exit_on_err:
            raise SystemExit(exit_on_err)

    return retcode


@dataclass
class PathQuery:
    analysis: str
    cfg: str
    pair: str
    cent: str
    kt: str
    magfield: str

    def mean_kt(self):
        from statistics import mean
        return mean(map(float, self.kt.split("_")))

    def estimate_gamma(self):
        """
        Return the estimated gamma in the transverse (out) direction

        $\\gamma_T^2
            = \\frac{p_T^2 + m^2}{m^2}
            = \\frac{4k_T^2 + m^2}{m^2}$
        """
        kt = self.mean_kt()
        return pd.np.sqrt(1.0 + (2 * kt / 0.13957) ** 2)

    def as_path(self):
        return '/'.join(self)

    @classmethod
    def from_path(cls, path):
        if isinstance(path, str):
            path = path.split('/')
        return cls(*path)

    def as_dict(self):
        keys = ('analysis', 'cfg', 'pair', 'cent', 'kt', 'magfield')
        return {key: getattr(self, key) for key in keys}

    @classmethod
    def from_dict(cls, data):
        from dataclasses import fields
        keys = {f.name for f in fields(cls)}

        return cls(**{k: v
                      for k, v in data.items()
                      if k in keys})

    def find(self, df, *keys):
        """
        Return dataframe matching selected keys from this query
        """
        from functools import reduce
        from operator import and_

        if not keys:
            keys = ('analysis', 'cfg', 'pair', 'cent', 'kt', 'magfield')

        query = reduce(and_, (df[key] == getattr(self, key) for key in keys))

        return df[query]

    @classmethod
    def From(cls, obj):
        from pathlib import Path

        if isinstance(obj, cls):
            return obj

        if isinstance(obj, dict):
            return cls.from_dict(obj)

        if isinstance(obj, str):
            return cls.from_path(obj)

        if isinstance(obj, Path):
            return cls.from_path(str(obj))

        if isinstance(obj, pd.Series):
            from dataclasses import fields
            return cls(*obj[[f.name for f in fields(cls)]])

        if isinstance(obj, pd.DataFrame):
            return (cls.From(r[1]) for r in obj.iterrows())

        from ROOT import TDirectory
        if isinstance(obj, TDirectory):
            path = obj.GetName()
            return cls.from_path(path)

        raise TypeError(f"Cannot build pathquery from {obj.__class__}")

    def __call__(self, **kwargs):
        from copy import copy
        res = copy(self)
        for k, v in kwargs.items():
            setattr(res, k, v)
        return res

    def __iter__(self):
        yield from self.values()

    def keys(self):
        from dataclasses import fields
        yield from (f.name for f in fields(self))

    def values(self):
        yield from (getattr(self, f) for f in self.keys())

    def items(self):
        yield from zip(self.keys(), self.values())


def get_momentum_resolution_correction_map(path='mrcdata.yaml'):
    import yaml

    with open(path) as f:
        data = yaml.load(f)

    # invert the dict
    result = {}
    for cfg, datas in data.items():
        # implicit identity
        result[cfg] = cfg

        for d in datas:
            if d in result:
                print(f"WARNING: `{d}` specified multiple times")
            result[d] = cfg

    return result


def flatten_config(cfg: dict, delim: str = '.') -> dict:
    result = {}
    for key, val in cfg.items():
        if isinstance(val, dict):
            for subkey, subval in flatten_config(val).items():
                result[f'{key}{delim}{subkey}'] = subval
        else:
            result[key] = val
    return result


def fitresult_data_df(json_data):
    return FitResults.get_data_df(json_data)


def fitresult_config_df(json_data):
    return FitResults.get_config_df(json_data)


class FitResults:

    # tfile cache
    _opened_files = {}

    def __init__(self, path):

        with open(path) as f:
            data = json.load(f)

        self.df = self.extract_data_df(data)
        self.config = self.extract_config_df(data)
        self._data = data

    def __getattr__(self, key):
        if key in self._data:
            return self._data[key]

        raise AttributeError

    def get_merged_df(self):
        return pd.merge(self.df, self.config, on='cfg')

    def print_row_summary(self, row):
        dat = self.df.loc[row]

        fitter = 'FitterGauss6'
        keys = {
            'FitterGauss6': ("Ro", "Rs", "Rl", "lam", "Ros"),
        }[fitter]
        print(f" | kT: {dat.kt} cent: {dat.cent}  pair: {dat.pair} ")
        print(f" | subset: {dat.subset or None}")
        print(" |-------- ")
        for key in keys:
            print(f" | {key}: {dat[key]:.4g} ± {dat[key + '_err']:.4g}")

    def data_from_row(self, row):
        query = PathQuery.From(self.df.loc[row])

        tfile = self.get_file(self._data['filename'])
        tdir = tfile.Get(query.as_path())
        return tdir

    @staticmethod
    def merge_table(tbl, *keys, xkey='kT', as_dataframe=True):
        result = {xkey: pd.np.array(tbl.index)}

        for key in keys:
            v = tbl[key]
            e = tbl[key + '_err']

            weights = 1.0 / e
            vals = (v * weights).sum(axis=1) / weights.sum(axis=1)
            err = pd.np.sqrt((e**2).sum(axis=1))
            result[key] = pd.np.array(vals)
            result[key + '_err'] = pd.np.array(err)

        if as_dataframe:
            return pd.DataFrame(result)

        return result

    def combined_data(self,
                      df=None,
                      group='cent',
                      xkey='kT',
                      keys=("Ro", "Rs", "Rl", "lam"),
                      cols=('pair', 'magfield')):
        """
        Merge data correctly with statistical errors
        """

        keys = set(keys)

        pivot_keys = keys | {key + "_err" for key in keys}

        df = df or self.df

        tbls = []
        for subkey, sdf in df.groupby(group):
            pivot = sdf.pivot_table(index=xkey,
                                    columns=cols,
                                    values=pivot_keys)
            tbl = self.merge_table(pivot, *keys)
            if isinstance(xkey, str):
                tbl.insert(0, group, subkey)
            tbls.append(tbl)

        return pd.concat(tbls).reset_index(drop=True)

    @classmethod
    def get_file(cls, filename):
        from ROOT import TFile

        try:
            return cls._opened_files[filename]
        except KeyError:
            pass

        tfile = TFile.Open(filename)
        if not tfile:
            raise FileNotFoundError(f'Requested file {filename}')
        cls._opened_files[filename] = tfile
        return tfile

    @staticmethod
    def extract_data_df(json_data):
        data = json_data['df'] if 'df' in json_data else json_data
        return pd.DataFrame(data)

    @staticmethod
    def extract_config_df(json_data):
        config = json_data['config'] if 'config' in json_data else json_data
        result = []

        if isinstance(config, list):
            return pd.DataFrame(config)

        for v, c in config.items():
            flat = flatten_config(c)
            flat['cfg'] = v.partition("/")[2]
            result.append(flat)

        return pd.DataFrame(result).set_index('cfg')


def unique_histnames():
    i = 1
    while True:
        yield "hist_%d" % i
        i += 1


histname = unique_histnames()


def yield_projections(*hists, lim=0.2, scale=False):
    from ROOT import TH3
    projections = TH3.ProjectionX, TH3.ProjectionY, TH3.ProjectionZ
    axis_getters = TH3.GetXaxis, TH3.GetYaxis, TH3.GetZaxis

    for project, get_axis in zip(projections, axis_getters):
        results = ()
        for h in hists:
            ax = get_axis(h)
            lobin, hibin = map(ax.FindBin, (-lim, lim))
            ph = project(h, next(histname), *(lobin, hibin)*2)
            ph.SetStats(False)
            if scale:
                ph.Scale(1.0 / (hibin - lobin + 1) ** 2)
            results += (ph, )
        yield results


def normalize_yaxis(canvas, offset_percent=2):
    from ROOT import TH1

    lofact = 1 - offset_percent / 100
    hifact = 1 + offset_percent / 100

    pads = canvas.GetListOfPrimitives()
    hists = [pad.GetListOfPrimitives().At(0) for pad in pads]
    hists = [h for h in hists if isinstance(h, TH1)]
    range_min, range_max = 100, 0

    range_min = min(range_min, *(h.GetMinimum() for h in hists))
    range_max = max(range_max, *(h.GetMaximum() for h in hists))

    for h in hists:
        h.GetYaxis().SetRangeUser(range_min * lofact, range_max * hifact)


def interpolate_mrc_hist(h, mrc, warn=False):
    xaxis = h.GetXaxis()
    yaxis = h.GetYaxis()
    zaxis = h.GetZaxis()

    for k in range(zaxis.GetFirst(), zaxis.GetLast() + 1):
        z = zaxis.GetBinCenter(k)
        for j in range(yaxis.GetFirst(), yaxis.GetLast() + 1):
            y = yaxis.GetBinCenter(j)
            for i in range(xaxis.GetFirst(), xaxis.GetLast() + 1):
                x = xaxis.GetBinCenter(i)

                mbin = mrc.FindBin(x, y, z)
                if mrc.IsBinUnderflow(mbin) or mrc.IsBinOverflow(mbin):
                    if warn:
                        print("bin", mbin, "is outside mrc range")
                    continue

                f = mrc.Interpolate(x, y, z)
                # print("factor", f)
                h.SetBinContent(i, j, k, f * h.GetBinContent(i, j, k))
