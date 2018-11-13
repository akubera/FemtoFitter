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
        if isinstance(obj, cls):
            return obj

        if isinstance(obj, dict):
            return cls.from_dict(obj)

        return cls.from_path(obj)


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
