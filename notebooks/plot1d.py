
import numpy as np
import pandas as pd


class PlotData:

    def __init__(self, canvas=None):
        self.canvas = canvas

    def Draw(self, opts=''):
        self.canvas.Draw(opts)


def merge_points(df):
    x = []
    y = []
    ye = []

    for kt, subdf in df.sort_values('kT').groupby('kT'):
        # x.append(subdf.kT.mean())
        x.append(subdf.kT.iloc[0])
        y.append(subdf.radius.mean())
        ye.append(subdf.radius_err.mean())

    data = pd.DataFrame([x,y,ye], ['kT', 'radius', 'radius_err']).T
    return data


def make_1d_plots(df, c=None):
    if c is None:
        from ROOT import TCanvas
        c = TCanvas()
        c.Divide(1, 2)
    import seaborn as sns
    import ROOT
    from ROOT import TGraphErrors, TLegend
    plot = PlotData()
    plot.canvas = c

    graphs = plot.graphs = []
    graphs_lam = plot.graphs_lam = []

    leg = plot.legend = TLegend(0.68, 0.6, 0.88, 0.8)
    leg.SetHeader("Centrality", 'C')

#     colors = [ROOT.kRed, ROOT.kBlue, ROOT.k]
    sns_colors = sns.color_palette("colorblind")
    tcolors = plot.tcolors = [ROOT.TColor(ROOT.TColor.GetFreeColorIndex(), *c) for c in sns_colors]
    colors = plot.colors = [c.GetNumber() for c in tcolors]


    for i, (cent, cdf) in enumerate(df.groupby('cent')):
        color = colors[i]
        data = merge_points(cdf)

        g_rinv = TGraphErrors(data.shape[0])
        np.frombuffer(g_rinv.GetX())[:] = np.array(data.kT)
        np.frombuffer(g_rinv.GetY())[:] = np.array(data.radius)
        np.frombuffer(g_rinv.GetEY())[:] =  np.array(data.radius_err) * 10

        g_lam = TGraphErrors(data.shape[0])
        np.frombuffer(g_lam.GetX())[:] = np.array(data.kT)
        np.frombuffer(g_lam.GetY())[:] = np.array(data.lam)
        np.frombuffer(g_lam.GetEY())[:] =  np.array(data.lam_err) * 10


        for g in (g_rinv, g_lam):
            g.SetMarkerStyle(21)
            g.SetMarkerSize(0.7)
            g.SetMarkerColor(color)
            g.SetLineColor(color)

        graphs.append(g_rinv)
        graphs_lam.append(g_lam)

        c.cd(1)
        if i == 0:
            g_rinv.Draw("APE")
        else:
            g_rinv.Draw("P ")

        c.cd(2)
        if i == 0:
            g_lam.Draw("APE")
        else:
            g_lam.Draw("P ")

        cent_name = cent.replace('_', '-')  + "%"
        leg.AddEntry(g_rinv, cent_name, 'P')


#     for g in graphs:
#         g.Draw('AP')
#     g.Draw('APL')

#     g.SetTitle("Radius")
    leg.Draw()
    return plot


def load_joeys_style():
    from ROOT import gStyle, gROOT
    font_id = 12
    gStyle.SetOptStat(0)
    gStyle.SetTitleFont(font_id, "X") # LaTeX typeface
    gStyle.SetTitleFont(font_id, "Y")
    gStyle.SetTitleFont(font_id, "Z")
    gStyle.SetTitleFont(font_id, "T")
    gStyle.SetLabelFont(font_id, "X") # LaTeX typeface
    gStyle.SetLabelFont(font_id, "Y")
    gStyle.SetLabelFont(font_id, "Z")
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTopMargin(0.13)
    gStyle.SetPadBottomMargin(0.25)
    gStyle.SetPadLeftMargin(0.17)
    gStyle.SetPadRightMargin(0.25)
    gStyle.SetCanvasDefW(1500)
    gStyle.SetCanvasDefH(900)
    gStyle.SetTitleSize(0.06,"X")
    gStyle.SetTitleSize(0.06,"Y")
    gStyle.SetTitleSize(0.06,"Z")
    gStyle.SetTitleSize(0.06,"T") # for title
    gStyle.SetLabelSize(0.04,"X")
    gStyle.SetLabelSize(0.04,"Y")
    gStyle.SetLabelSize(0.04,"Z")
    gStyle.SetTitleOffset(1.5,"X")
    gStyle.SetTitleOffset(1.1,"Y")
    gStyle.SetTitleOffset(1.3,"Z")
    gROOT.ForceStyle()


def make_1d_correlation_function(df, key, tfile=None, size=(800, 600), c=None):
    import ROOT
    from ROOT import TCanvas, TLine
    from femtofitter import PathQuery
    if c is None:
        c = TCanvas()
        c.SetCanvasSize(*size)

    if tfile is None:
        raise NotImplementedError

    if isinstance(key, int):
        series = df.iloc[key]
    else:
        raise NotImplementedError

    q = PathQuery.From(series)
    tdir = tfile.Get(q.as_path())
    num, den = map(tdir.Get, ("num", "den"))

    norm_rng = [num.GetXaxis().FindBin(l) for l in (0.25, 0.3)]

    ratio = num.Clone("ratio")
    if ratio.GetSumw2N() == 0:
        ratio.Sumw2()
    ratio.Divide(num, den, den.Integral(*norm_rng), num.Integral(*norm_rng))
    ratio.SetStats(0)

    ratio.SetTitle("Correlation Function; q_{inv} (GeV); CF(q_{inv})")
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(0.6)

    xax, yax = ratio.GetXaxis(), ratio.GetYaxis()
    xax.SetRangeUser(0.0, 0.3)
    xax.SetTitleSize(0.040)
    yax.SetTitleOffset(1.0)
    yax.SetTitleOffset(1.0)
    yax.CenterTitle()

    plot = PlotData()
    plot.ratio = ratio
    ratio.Draw('P')
    plot.canvas = c

    plot.unity_line = TLine(0, 1.0, 0.3, 1.0)
    plot.unity_line.Draw()
    plot.unity_line.SetLineStyle(2)
    plot.unity_line.SetLineColor(ROOT.kGray+2)

    # draw fit
    if False:
        fitter_classname = {
           'Gauss1D': "Fitter1DGauss",
        }[series.fitter]
        fitter = getattr(ROOT, fitter_classname)
        print(fitter.FitResult(series))

    return plot


def plot_ratio_data(df1, df2, data, c=None):
    from ROOT import TH2D
    if c is None:
        from ROOT import TCanvas
        c = TCanvas()

    plot = PlotData()
    plot.canvas = c
    nktbins = 5
    ncentbins = 2
    h = TH2D('hist', 'Hist; kT (GeV); Centrality',
             nktbins, -0.5, nktbins+0.5,
             ncentbins, -0.5, ncentbins+0.5)
    h.SetStats(0)

    cents = []
    kts = []
    ratios = []

    for cent, cdf in df1.groupby('cent'):
        cents.append(cent)
        for kt, ktdf in cdf.groupby('kt'):
            pass

#     h.FillRandom()
#     h2 = h.Clone('h2')
#     h2.Fill

#     h.Divide(h2)
    np.frombuffer(h.GetArray(), count=h.GetNcells())[:] = data.flatten()
    h.SetBinContent(0,0)
    h.Draw('COLZ')

    plot.h = h
    return plot


def plot_projected_fits_abc(fr):
    pass


def plot_projected_fits(fr,
                        centrality,
                        N=4,
                        size=(600, 800),
                        ylim=(0.85, 1.44),
                        c=None,
                        I=1,
                        file_cache={}):
    """
    Plot big matrix of historams and fits
    """
    import re
    import ROOT
    from femtofitter import PathQuery
    from ROOT import TFile, gStyle, gROOT
    from ROOT import Data3D, Fitter3DGaussLcms as Fitter

    if c is None:
        from ROOT import TCanvas
        c = TCanvas()
        c.SetCanvasSize(*size)

    plot = PlotData(c)

    tfile = TFile.Open(fr.filename)
    assert tfile

    integration_limits = None

    def make_projections(num, den):
        nonlocal integration_limits
        xax = num.GetXaxis()
        yax = num.GetYaxis()
        zax = num.GetZaxis()

        if not integration_limits:
            xzero = xax.FindBin(0.0)
            yzero = yax.FindBin(0.0)
            zzero = zax.FindBin(0.0)

            xrange = (xzero - N, xzero + N) if xzero > 1 else (1, N)
            yrange = (yzero - N, yzero + N)
            zrange = (zzero - N, zzero + N)

            integration_limits = [
                (xax.GetBinLowEdge(xrange[0]), xax.GetBinUpEdge(xrange[1])),
                (yax.GetBinLowEdge(yrange[0]), yax.GetBinUpEdge(yrange[1])),
                (zax.GetBinLowEdge(zrange[0]), zax.GetBinUpEdge(zrange[1]))]
            # print(integration_limits)
        else:
            xrange, yrange, zrange = integration_limits
            xrange = tuple(map(xax.FindBin, xrange))
            yrange = tuple(map(yax.FindBin, yrange))
            zrange = tuple(map(zax.FindBin, zrange))

        name = 'hist%08x' % (abs(hash(num)))

        rx = num.ProjectionX(f'ratiox{name}', *yrange, *zrange)
        ry = num.ProjectionY(f'ratioy{name}', *xrange, *zrange)
        rz = num.ProjectionZ(f'ratioz{name}', *xrange, *yrange)

        for r in (rx, ry, rz):
            r.SetStats(0)
            r.Sumw2()
            r.SetTitle('')

        dx = den.ProjectionX(f'denx{name}', *yrange, *zrange)
        dy = den.ProjectionY(f'deny{name}', *xrange, *zrange)
        dz = den.ProjectionZ(f'denz{name}', *xrange, *yrange)

        rx.Divide(dx)
        ry.Divide(dy)
        rz.Divide(dz)

        return rx, ry, rz

    def load_mrc(mrc_paths):
        mrc_ptrs = {}
        for mrc in mrc_paths:
            m = re.match(r'(?P<cls>\w+)\[(?P<file>[^:]+):(?P<path>[^\]]+)\]', mrc)
            d = m.groupdict()
            try:
                mrc_file = file_cache[d['file']]
            except KeyError:
                mrc_file = file_cache[d['file']] = TFile.Open(d['file'])

            mrc_cls = getattr(ROOT, d['cls'])
            mrc_tdir = mrc_file.Get(d['path'])
            mrc_ptrs[mrc] = mrc_cls.From(mrc_tdir)
        return mrc_ptrs

    def load_fsi(fsi_paths):
        fsi_objs = {}
        for fsi in fsi_paths:
            m = re.match(r'(?P<cls>\w+)\[(?P<file>[^\]]+)\]', fsi)
            d = m.groupdict()
            try:
                fsi_file = file_cache[d['file']]
            except KeyError:
                fsi_file = file_cache[d['file']] = TFile.Open(d['file'])

            fsi_cls = getattr(ROOT, d['cls'])
            fsi_objs[fsi] = fsi_cls.From(fsi_file)

        return fsi_objs

    df = fr.df[fr.df.cent==centrality]

    mrc_dict = load_mrc(df.mrc.unique())
    fsi_dict = load_fsi(df.fsi.unique())

    hists = plot.hists = []
    fithists = plot.fithists = []
    fits = []
    kts = []

    plot.fit_tcolor = ROOT.TColor(1.0, 0.0, 0.0, 0.5)
    fit_color = plot.fit_tcolor.GetNumber()
    fsi = ROOT.FsiStatic(1.0)
    for kt, ktdf in df.groupby('kt'):
        kts.append(kt)

        series = ktdf.iloc[I]
        query = PathQuery.From(series)
        fitres = Fitter.FitResult(series)
        fits.append(fitres)

        tdir = tfile.Get(query.as_path())
        data = Data3D.From(tdir)
        fitter = Fitter(data)
        fitter.fsi = fsi_dict[series.fsi]
        fitter.mrc = mrc_dict[series.mrc]

        num, den = map(tdir.Get, ("num", "den"))

        rx, ry, rz = make_projections(num, den)
        hists.append((rx, ry, rz))
        qinv = data.src.qinv.get()
#         fsi = fitter.fsi.get()

        fit_params = fitres.as_params()
        fit_params.norm = 1.0

        fitden = fitter.mrc.GetSmearedDen()
#         fitnum = fitter.mrc.GetSmearedFit(fit_params, qinv, fsi)
#         fitnum.Multiply(fitden)

#         np.frombuffer(fitden.GetArray(), np.float32, fitden.GetNcells())[:] = 1.0
        fitnum = fitden.Clone('fitnum')
        fit_params.fill(fitnum, qinv, fsi, 1)

        fitnum.SetLineColor(fit_color)
        fitnum.SetMarkerColor(fit_color)

        fx, fy, fz = make_projections(fitnum, fitden)
        fithists.append((fx, fy, fz))
        for f in fithists[-1]:
            f.SetLineWidth(2)
            f.Scale(1.0/(f.GetBinContent(f.GetNbinsX()) or 1))

    c.Divide(3, len(hists), 0, 0)

    for j, (hlist, flist) in enumerate(zip(hists, fithists)):
        scale = 1.0 / fits[j].norm.value
        for i, (h, f) in enumerate(zip(hlist, flist), 1):
            pad = c.cd(3 * j + i)
            pad.SetTicks(1, 1)
            h.Draw("H")
            h.Scale(scale)
            h.GetXaxis().SetLabelSize(0.12)
            h.GetYaxis().SetRangeUser(*ylim)
            h.GetYaxis().SetLabelSize(0.12)
            h.GetYaxis().SetNdivisions(30205)
            h.GetXaxis().SetNdivisions(5)
            if i == 1:
                pad.SetLeftMargin(0.2)
            elif i == 3:
                pad.SetRightMargin(0.2)
                h.GetYaxis().SetTitle(str(kts[j]).replace('_', ' - '))
                h.GetYaxis().CenterTitle()
                h.GetYaxis().RotateTitle()
                title_size = 0.12
                title_offset = -0.53 / title_size
                h.GetYaxis().SetTitleSize(title_size)
                h.GetYaxis().SetTitleOffset(title_offset)

            if j > len(hists) - 3:
                h.GetXaxis().SetTitleSize(title_size)
#                 pad.SetBottomMargin(0.5)

            f.Draw("SAME HIST L")

    return plot



