#
# femtofitter/__init__.py
#

from typing import NamedTuple


class PathQuery(NamedTuple):
    analysis: str
    cfg: str
    pair: str
    cent: str
    kt: str
    magfield: str

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
        return cls(**data)

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

        from ROOT import TDirectory
        if isinstance(obj, TDirectory):
            path = obj.GetName()
            return cls.from_path(path)

        raise TypeError(f"Cannot build pathquery from {obj.__class__}")


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


def flatten_config(cfg: dict, delim: str='.') -> dict:
    result = {}
    for key, val in cfg.items():
        if isinstance(val, dict):
            for subkey, subval in flatten_config(val).items():
                result[f'{key}{delim}{subkey}'] = subval
        else:
            result[key] = val
    return result


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
                ph.Scale(1.0 / (hibin - lowbin + 1) **2)
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
