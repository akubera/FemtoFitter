#
# femtofitter/plotting.py
#

import numpy as np
import pandas as pd

from dataclasses import dataclass

from typing import Optional


def projections3D():
    from ROOT import TH3

    yield TH3.ProjectionX
    yield TH3.ProjectionY
    yield TH3.ProjectionZ


def yield_axes(h):
    yield h.GetXaxis()
    yield h.GetYaxis()
    yield h.GetZaxis()


def plot_projected_cf(num, den,
                      lim=0.05,
                      c=None,
                      opt='HIST HE',
                      linecolor=None,
                      linestyle=None):
    from ROOT import TCanvas
    from stumpy.rhist import normalize_hist

    if c is None:
        c = TCanvas()
        c.SetCanvasSize(800, 800)
        c.Divide(2, 2)

    lb, hb = map(num.GetXaxis().FindBin, (-lim, lim))
    l = (lb, hb)*2

    projections = ((p(num, "pn", *l), p(den, "pd", *l))
                   for p in projections3D())

#     normrange = tuple(map(num.GetXaxis().FindBin, (-.11, -0.08)))

    for i, (n, d) in enumerate(projections, 1):
        c.cd(i)
        if n.GetSumw2N() == 0:
            n.Sumw2()
        n.Divide(d)
        n.SetStats(False)
        if linecolor is not None:
            n.SetLineColor(linecolor)
        if linestyle is not None:
            n.SetLineStyle(linestyle)
        normalize_hist(n, (-.11, -0.08))
        n.DrawCopy(opt)

    return c


def plot_projections(fit_results,
                     data,
                     mrc,
                     limit=0.1,
                     ylim=(0.98, 1.15),
                     c=None,
                     FitterClass=None):
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
    titles = ('q_{out};;', 'q_{side};;', 'q_{long};;')
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

        c.cd(i)
        r.DrawCopy("HE")

    RED = 2
    fithist = den.Clone(basename + "fithist")
    fithist.Divide(mrc)
    fithist.SetLineColorAlpha(RED, 0.7)

    if FitterClass is None:
        from ROOT import Fitter3DGaussLcms as FitterClass
    params = FitterClass.FitParams(fit_results)
    params.norm = 1.0
#     params.gamma = data.data.gamma

    qinv = data.get_qinv()
    params.apply_to(fithist, qinv, 6.5)

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


def normalize_canvas(c, d=0.02):
    from ROOT import TH1
    mins, maxs = [], []
    for obj in c.GetListOfPrimitives():
        if isinstance(obj, TH1):
            maxs.append(obj.GetMaximum())
            mins.append(obj.GetMinimum())

    gmin = (1 - d) * min(mins)
    gmax = (1 + d) * max(maxs)

    for obj in c.GetListOfPrimitives():
        if isinstance(obj, TH1):
            obj.GetYaxis().SetRangeUser(gmin, gmax)


def normalize_subcanvases(c, d=0.02):
    from ROOT import TH1
    mins, maxs = [], []
    for pad in c.GetListOfPrimitives():
        for obj in pad.GetListOfPrimitives():
            if isinstance(obj, TH1):
                maxs.append(obj.GetMaximum())
                mins.append(obj.GetMinimum())

    gmin = (1 - d) * min(mins)
    gmax = (1 + d) * max(maxs)

    for pad in c.GetListOfPrimitives():
        for obj in pad.GetListOfPrimitives():
            if isinstance(obj, TH1):
                obj.GetYaxis().SetRangeUser(gmin, gmax)


def set_titles(canvas, *titles):
    from ROOT import TH1
    title_iterator = iter(titles)
    for pad in canvas.GetListOfPrimitives():
        for obj in pad.GetListOfPrimitives():
            if isinstance(obj, TH1):
                obj.SetTitle(next(title_iterator))
                break


def stat_mean(vals, errs):
    """
    Correct statistical error combining
    """
    weights = errs ** -2
    value = (vals * weights).sum() / weights.sum()
    error = np.sqrt(1.0 / weights.sum())
    return value, error


def extract_values(df, ykey, ekey=None, xkey='kT', combine=True):
    # ensure x-variable is ordered
    odf = df.sort_values(xkey)
    if ekey is None:
        ekey = ykey + "_err"

    if not combine:
        return np.array(odf[[xkey, ykey, ekey]]).T

    res = np.array([(x, *stat_mean(d[ykey], d[ekey])) for x, d in odf.groupby(xkey)])
    return res.T


class QuadPlot:

    @dataclass
    class ShowParams:
        linestyle: str = '--'
        xlim = (0.2, 1.2)
        legend_axis = ''
        capstyle: str = 'butt'

    def __init__(self, fr):
        from femtofitter import FitResults
        if isinstance(fr, FitResults):
            self.df = fr.df
        else:
            self.df = fr

    def show(self,
             *groups,
             limit_alpha=True,
             cfg=None,
             sysdf=None,
             xshift=None,
             combine=True,
             params: Optional[ShowParams] = None,
             ):
        """
        Create and display quad-plot
        """

        import matplotlib.pyplot as plt
        from itertools import repeat

        if cfg is not None:
            df = self.df[self.df.cfg == cfg]
        else:
            df = self.df

        if params is None:
            params = self.ShowParams()

        default_plot_ops = {
            'linestyle': params.linestyle,
            'marker': '.',
            'fmt': '',
            'capsize': 4,
            'solid_capstyle': params.capstyle,
            'xlim': params.xlim,
        }

        plt.close()

        keys = ["Ro", "Rs", "Rl", 'lam']
        if 'alpha' in df.columns:
            keys.append('alpha')

        legend_key = params.legend_axis

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

                    X, Y, E = extract_values(cent_data, key, combine=combine)

                    if sysdf is not None:
                        subsysdf = sysdf[sysdf.cent == cent]
                        E + subsysdf[key + "_sys"]

                    xlim = plot_ops.pop("xlim")
                    plot_ops['label'] = cent
                    plot_ops['lw'] = 2

                    X += shift
                    eb = ax.errorbar(X, Y, yerr=E, **plot_ops)
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
            result = [(grp_val, _do_makeplot(grp_data))
                      for group in groups
                      for grp_val, grp_data in self.df.groupby(group)]
        else:
            result = _do_makeplot(df)

        return result


def ratio_tree_canvas(*args, sli='z', width=0, plotsize=(200, 200)):
    """
    Create tree of histogram ratios
    """

    if len(args) == 0:
        pass

    c = TCanvas()

    pairs0 = [iter(args)] * 2

    ratios0 = []

    zbin = ng.GetXaxis().FindBin(0.0)
    start = zbin + width
    stop = ng.GetNbinsX() + width

    def _plot_2d_projection(h, title='', zrange=None):
        zax = h.GetZaxis()
        zax.SetRange(zbin, zbin)
        zz = h.Project3D("yx")
        zz.SetStats(False)
        zz.SetTitle(title)
        if zrange:
            zz.GetZaxis().SetRangeUser(*zrange)
        zz.DrawCopy("COLZ")

    for num, den in zip(*pairs0):

        ratios0.append(num.Clone("ratio_%s" % num.GetName()))
        ratios0[-1].Divide(den)

    s = "abc"
    s.upper()


    ratios

    ng, nr, dg, dr = map(datadir.Get, ("NumGen", "NumRec", "DenGen",  "DenRec"))
#     for i in range(start, stop):
#         obin = zbin - (i - zbin)
#         for j in range(1, stop):
#             oj = zbin - (j - zbin)
#             for h in (ng, nr, dg, dr):
#                 ls = h.GetBinContent(i, j, zbin)
#                 rs = h.GetBinContent(obin, oj, zbin)
#                 h.SetBinContent(i, j, zbin, ls + rs)
#                 h.SetBinContent(obin, oj, zbin, ls + rs)


    cf_g = ng.Clone("cfg_%d" % i)
    cf_g.Divide(dg)

    cf_r = nr.Clone("cfr_%d" % i)
    cf_r.Divide(dr)

    c = TCanvas()
    c.Divide(3, 1)
    c.SetCanvasSize(1200, 700)
    canvases.append(c)


    llpad = c.cd(1)
    llpad.Divide(1, 4)
    for i, h in enumerate((ng, dg, nr,dr), 1):
        llpad.cd(i)
        _plot_2d_projection(h)


    leftpad = c.cd(2)
    leftpad.Divide(1, 2)
    leftpad.cd(1)
    _plot_2d_projection(cf_g)

    leftpad.cd(2)
    _plot_2d_projection(cf_r)

    normalize_subcanvases(leftpad)

    rpad = c.cd(3)
    mrc = cf_g.Clone("mrc_%d" % i)
    mrc.Divide(cf_r)
#     rpad.SetBBoxY1(150)
    rpad.SetBBoxY1(200)
    rpad.SetBBoxCenterY(250)

#     rpad.SetBBoxY2(200)
#     rpad.SetBBoxY1(int(350*0.5))
#     rpad.SetBBoxCenterY(350)
#     rpad.SetBBoxCenterY(4300 - 350//2)

    _plot_2d_projection(mrc, 'MRC (#Delta#eta: %g, #Delta#phi: %g)' % cut, (0.6, 1.2))

    canvases.append(TCanvas())

    d_ratio = dg.Clone()
    d_ratio.Divide(dr)
    _plot_2d_projection(d_ratio, 'Denominator Ratio (Gen / Rec)')

    canvases.append(TCanvas())

    n_ratio = datadir.Get("DenGenWeight")
    n_ratio.Divide(dg)
    n_ratio.Divide(dr)
    n_ratio.Multiply(datadir.Get("DenRecWeight"))

#     n_ratio = ng.Clone()
#     n_ratio.Divide(nr)
    _plot_2d_projection(n_ratio, '"Fake" MRC', (.0, 1.2))


    return c


def build_TGraphErrors(x, y, ye=None, xe=None):
    from ROOT import TGraphErrors
    if xe is None:
        xe = np.zeros_like(x)
    else:
        xe = np.array(xe)

    if ye is None:
        ye = np.zeros_like(x)
    else:
        ye = np.array(ye)

    x = np.array(x)
    y = np.array(y)
    graph = TGraphErrors(x.size, x, y, xe, ye)
    return graph


def series_to_TGraphErrors(series, ykey, ekey=None, xkey='kT'):
    ekey = ekey or ykey + '_err'

    graph = build_TGraphErrors(series[xkey], series[ykey], series[ekey])
    return graph


def plot_outside(num, den, N=0, opts='COLZ', pad=None, norm=True):
    if pad is None:
        from ROOT import TCanvas
        pad = TCanvas()

    zbin = num.GetZaxis().FindBin(0)
    zrng = zbin-N,zbin+N

    zaxn = num.GetZaxis()
    zaxn.SetRange(*zrng)
    n = num.Project3D("yx")
    zaxn.SetRange()

    zaxd = den.GetZaxis()
    zaxd.SetRange(*zrng)
    d = den.Project3D("yx")
    zaxd.SetRange()

    n.SetStats(False)
    n.SetName("ratio")

    n.Divide(n, d)
    if norm:
        n.Scale(den.Integral() / num.Integral())

    zmax = num.GetZaxis().GetBinCenter(zrng[1])
    n.SetTitle("CF(q) #scale[0.85]{|q_{long}| < %0.1f MeV}" % (zmax * 1000))
    n.SetName("outside")
    n.Draw(opts)

    return pad

def plot_outside_tdir(tdir, N=0, pad=None, opts="COLZ", norm=True):
    keylist = [
        ("Num", "Den"),
        ("num", "den"),
    ]

    for keys in keylist:
        num, den = map(tdir.Get, keys)
        if num and den:
            break
    else:
        raise ValueError("Could not find num and den in", tdir)

    return plot_outside(num, den, N=N, pad=pad, opts=opts, norm=norm)


def plot_(tdir, N=0, pad=None, opts="COLZ", norm=True):
    keylist = [
        ("Num", "Den"),
        ("num", "den"),
    ]

    for keys in keylist:
        num, den = map(tdir.Get, keys)
        if num and den:
            break
    else:
        raise ValueError("Could not find num and den in", tdir)

    return plot_outside(num, den, N=N, pad=pad, opts=opts, norm=norm)
