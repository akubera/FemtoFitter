#
# femtofitter/fit1d.py
#
"""
Fitting routines for 1D correlation functions
"""


from .fit import parallel_fit_all, run_fit


def run_fit_gauss1d(*args, **kwargs):
    from ROOT import Data3D, FitterGauss1D
    return run_fit(FitterGauss1D, *args, **kwargs)


def _main(argv=None):
    if argv is None:
        import sys
        argv = sys.argv[1:]

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input",
                        help='Input rootfiles')
    parser.add_argument("output",
                        nargs='?',
                        default=None,
                        help='Input rootfiles')
    args = parser.parse_args(argv)

    from ROOT import TFile
    tfile = TFile.Open(args.input)
    if not tfile:
        return 1

    if args.output:
        dest = args.output
    else:
        dest = 'fitres-%s.json' % args.output

    parallel_fit_all(tfile, dest)


if __name__ == "__main__":
    exit(_main())
