#
# femtofitter/fitspec.py
#


import yaml
from pathlib import Path


class FitSpec:
    """
    """

    @classmethod
    def from_path(cls, path: Path):
        path = Path(path)

        with path.open('r') as f:
            return cls.from_file(f)

    @classmethod
    def from_file(cls, fd):
        data = yaml.load(fd)

        self = cls()
        self.data = data
        return self


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('filename', help='filename')
    args = parser.parse_args()

    # s = FitSpec.FromPath(args.filename)
    s = FitSpec.from_path(args.filename)

    print(vars(s))
