#
# femtofitter/plotting.py
#

import numpy as np
import pandas as pd


def projections3D():
    from ROOT import TH3

    yield TH3.ProjectionX
    yield TH3.ProjectionY
    yield TH3.ProjectionZ


def plot_projections(fit_results, data, mrc, limit=0.1, ylim=(0.98, 1.15), c=None):
    import ROOT
    from ROOT import TFile, TCanvas, TH3
    from random import randint

    if c is None:
        c = TCanvas()
    c.Divide(3, 1, 0, 0)

    num = data.get_numerator()
    den = data.get_denominator()

    axis = num.GetXaxis()
    binlo, binhi = map(axis.FindBin, (-limit, limit))
    limits = (binlo, binhi)*2

    try:
        norm = fit_results['norm']
    except TypeError:
        norm = fit_results.norm

    if not isinstance(norm, float):
        norm = norm.value

    scale_factor = 1.0 / norm

    basename = 'h%03d' % randint(0, 10000)
    names = ('_x', '_y', '_z')
    titles = ('q_{out}', 'q_{side}', 'q_{long}')
    for i, (name, projection) in enumerate(zip(names, projections3D()), 1):
        n = projection(num, basename + 'n' + name, *limits)
        d = projection(den, basename + 'd' + name, *limits)
        r = n.Clone(basename + "r" + name)
        if r.GetSumw2N() == 0:
            r.Sumw2()
        r.SetStats(False)
        r.Divide(d)
        r.Scale(scale_factor)
        r.GetYaxis().SetRangeUser(*ylim)
        r.SetTitle(titles[i-1])
        r.SetTitleSize(29)

        c.cd(i)
        r.DrawCopy("HE")

    RED = 2
    fithist = den.Clone(basename + "fithist")
    fithist.Divide(mrc)
    fithist.SetLineColorAlpha(RED, 0.7)

    from ROOT import FitterGaussOSL
    params = FitterGaussOSL.FitParams(fit_results)
    params.norm = 1.0
#     params.gamma = data.data.gamma

    qinv = data.get_qinv()
    params.apply_to(fithist, qinv)

    for i, (name, projection) in enumerate(zip(names, projections3D()), 1):
        f = projection(fithist, basename + 'f' + name, *limits)
        d = projection(den, basename + 'd' + name, *limits)

        if f.GetSumw2N() == 0:
            f.Sumw2()
        f.Divide(d)
        f.SetLineWidth(2)
        f.SetLineColorAlpha(RED, 0.5)

        c.cd(i)
        f.DrawCopy("HIST C SAME")

    c.Draw()
    return c


def stat_mean(vals, errs):
    """
    Correct statistical error combining
    """
    weights = errs ** -2
    value = (vals * weights).sum() / weights.sum()
    error = np.sqrt(1.0 / weights.sum())
    return value, error


def extract_values(df, ykey, ekey=None, xkey='kT', ):
    # ensure x
    odf = df.cent_data.sort_values(xkey)
    if ekey is None:
        ekey = ykey + "_err"

    res = np.array([(x, *stat_mean(d[ykey], d[ekey])) for x, d in odf.groupby(xkey)])
    return res.T


class QuadPlot:

    def __init__(self, fr):
        from . import FitResults
        if isinstance(fr, FitResults):
            self.df = fr.df
        else:
            self.df = fr

    def show(self, *groups, limit_alpha=True, cfg=None, sysdf=None, xshift=None):
        import matplotlib.pyplot as plt
        from itertools import repeat

        if cfg is not None:
            df = self.df[self.df.cfg == cfg]
        else:
            df = self.df

        default_plot_ops = {
            "linestyle": "--",
            'marker': 's',
            'fmt': '',
            'capsize': 4,
#             'solid_capstyle': 'butt',
            'xlim': (0.2, 1.0),
        }

        plt.close()

        keys = ["Ro", "Rs", "Rl", 'lam']
        if 'alpha' in df.columns:
            keys.append('alpha')

        legend_key = 'lam'

        if xshift is None:
            xshift = repeat(0)

        def _do_makeplot(data):
            fig, axs = plt.subplots(2, 2, figsize=(12, 12))
            xshift_it = iter(xshift)
            for cent, cent_data in data.groupby("cent"):
                shift = next(xshift_it)
                cent = cent.replace('_', '-') + "%"
                for ax, key in zip(axs.flat, keys):
                    plot_ops = default_plot_ops.copy()
                    if limit_alpha and key == 'alpha':
                        plot_ops['ylim'] = (1.0, 2.2) if limit_alpha is True else limit_alpha

                    if key == 'lam':
                        title = r"$\lambda$"
                    elif key == 'Ro':
                        title = "$R_{out}$"
                    elif key == 'Rs':
                        title = "$R_{side}$"
                    elif key == 'Rl':
                        title = "$R_{long}$"
                    else:
                        title = key

                    X, Y, E = extract_values(cent_data, key)

                    if sysdf is not None:
                        subsysdf = sysdf[sysdf.cent == cent]
                        E + subsysdf[key + "_sys"]

                    xlim = plot_ops.pop("xlim")
                    plot_ops['label'] = cent
                    plot_ops['lw'] = 2

                    ax.errorbar(X, Y, E, label=cent, **plot_ops)
                    X += shift
                    eb = ax.errorbar(X, Y, yerr=E, label=cent, lw=2, **plot_ops)
                    eb[-1][0].set_linestyle('--')

                    ax.set_title(title)
                    ax.set_xlim(*xlim)

                    if key == legend_key:
                        leg = ax.legend(numpoints=1, loc='best', fontsize=16)
                        if leg:
                            leg.set_title("Centrality", prop={"size": 16, 'weight': 'bold'})

            # for ax, key in zip(axs.flat, ("Ro", "Rs", "Rl", 'lam')):
            #     if key.startswith("R"):
            #         ax.set_ylim(2.0, 8.0)
            #     else:
            #         ax.set_ylim(0.2, 0.8)
            return fig

        if groups:
            result = [(group_val, _do_makeplot(pair_data))
                      for group_val, pair_data in map(self.df.groupby, groups)]
        else:
            result = _do_makeplot(df)

        return result
