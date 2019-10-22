
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
        c.Divide(r1)
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


def make_1d_correlation_function(df, key, tfile=None, c=None):
    import ROOT
    from ROOT import TCanvas, TLine
    from femtofitter import PathQuery
    if c is None:
        c = TCanvas()

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

#     h.FillRandom()
#     h2 = h.Clone('h2')
#     h2.Fill

#     h.Divide(h2)
    np.frombuffer(h.GetArray(), count=h.GetNcells())[:] = data.flatten()
    h.SetBinContent(0,0)
    h.Draw('COLZ')

    plot.h = h
    return plot
