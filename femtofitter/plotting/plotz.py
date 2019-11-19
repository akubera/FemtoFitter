#

import sys
from plot1d import PlotData


def make_multi_projection_plot(tdir, df=None, c=None, size=(500, 1000), N=5):
    from ROOT import TCanvas
    from stumpy.utils import walk_matching

    if c is None:
        c = TCanvas()

    plot = PlotData(c)
    kt_tdirs = list(walk_matching(tdir, '*_*/++'))

    c.SetCanvasSize(*size)
    c.Divide(3, len(kt_tdirs), 0, 0)

    def gen_random_names():
        from random import choices
        from string import ascii_letters
        while True:
            yield ''.join(choices(ascii_letters, k=12))

    names = gen_random_names()
    names.send(None)

    from ROOT import Fitter3DGaussLcms

    projections = (
        lambda h: h.ProjectionX(next(names), *ybins, *zbins),
        lambda h: h.ProjectionY(next(names), *xbins, *zbins),
        lambda h: h.ProjectionZ(next(names), *xbins, *ybins))

    for row, (kt, kt_tdir) in enumerate(kt_tdirs):

        n, d = map(kt_tdir.Get, ("num", 'den'))
        if n.GetSumw2N() == 0:
            n.Sumw2()

        if row == 0:
            if abs(n.GetXaxis().GetBinLowEdge(1)) < 1e-5:
                xbins = 1, N
            else:
                bin0 = n.GetXaxis().FindBin(0.0)
                xbins = bin0 - N, bin0 + N

            bin0 = n.GetYaxis().FindBin(0.0)
            ybins = bin0 - N, bin0 + N

            bin0 = n.GetZaxis().FindBin(0.0)
            zbins = bin0 - N, bin0 + N
            del bin0

        Fitter3DGaussLcms.FitResult()

        for col in range(3):
            pad = c.cd(row * 3 + (col+1))
            if row == 0:
                pad.SetTopMargin(0.2)

            np = projections[col](n)
            np.SetStats(0)
            np.SetTitle('')
            dp = projections[col](d)
            np.Divide(dp)
            h = np.DrawCopy('HE')
            h.SetTitleSize(0.06)
#             if row == 0:
#                 h.SetTitle(["Out", "Side", "Long" ][col])
#                 h.SetTitleSize(np.GetTitleSize() * 4.2)
    return plot



class Plotz:

    def __init__(self, tfile, path='PWG2FEMTO'):
        self.tfile = tfile
        self.container = tfile.Get(path).GetListOfKeys().At(0).ReadObj()
        assert self.container, f'No containers in path {path}'
        self.analysis = self.container.GetListOfKeys().At(0).ReadObj()

    def list_containers(self, path='PWG2FEMTO/*'):
        from stumpy.utils import walk_matching

        for _, tdir in walk_matching(self.tfile, path):
            yield tdir

    def set_analysis(self, cent=None, pair=None):
        *name, clo, chi, p = self.analysis.GetName().split('_')

        if (cent, pair) == (None, None):
            return

        if cent is None:
            cent = (clo, chi)
        elif isinstance(cent, str):
            cent = cent.split('_')
        if pair is None:
            pair = p

        analysis_name = '_'.join([*name, *cent, pair])

        analysis = self.container.Get(analysis_name)
        if not analysis:
            print(f"Could not find {analysis_name!r}", file=sys.stderr)
            raise ValueError

        self.analysis = analysis

    def list_analyses(self, path='PWG2FEMTO/*/*'):
        from stumpy.utils import walk_matching

        for _, tdir in walk_matching(self.tfile, path):
            yield tdir

    def get_pair_title(self):
        name = self.analysis.GetName()
        if 'pip' in name:
            return '#pi^{+}'

        if 'pim' in name:
            return '#pi^{-}'

        return '#pi^{?}'

    def build_zvertex_canvas(self, c=None):
        from ROOT import gStyle, gROOT, TCanvas

#         gStyle.SetLabelSize(0.05)
#         gROOT.ForceStyle()

        if c is None:
            c = TCanvas()

        event_dir = self.analysis.Get("Event/pass")
        plot = PlotData(c)

        vz = event_dir.Get("VertexZ")
        vz.SetStats(0)

        xax = vz.GetXaxis()
        yax = vz.GetYaxis()

        xax.SetLabelSize(0.04)
        xax.SetTitleOffset(0.9)
        xax.SetTitleSize(0.045)

        yax.SetTitleOffset(0.9)
        yax.SetTitleSize(0.045)
        yax.SetNdivisions(505)

        vz.SetTitle("Event Vertex Z-Component")
        vz.Draw()
        plot.vz = vz

        return plot

    def build_xyvertex_canvas(self, c=None):
        from ROOT import TCanvas

        if c is None:
            c = TCanvas()
            # c.SetCanvasSize(800, 500)

        plot = PlotData(c)

        event_dir = self.analysis.Get("Event/pass")
        vxy = event_dir.Get("VertexXY")
        vxy.SetStats(0)

        xax = vxy.GetXaxis()
        xax.SetRangeUser(0.055, 0.095)
        xax.SetNdivisions(908)

        yax = vxy.GetYaxis()
        yax.SetRangeUser(0.32, 0.352)
        yax.SetNdivisions(908)
        yax.SetTitleOffset(1.6)

        c.SetRightMargin(20.0)
        vxy.SetTitle("Event Vertex XY-components")
        vxy.Draw('COL')

        plot.vxy = vxy

        return plot

    def build_vertex_canvas(self, c=None, **opts):
        from ROOT import TCanvas
        if c is None:
            c = TCanvas()
            size = opts.get('size', (900, 400))
            c.SetCanvasSize(*size)

        plot = PlotData(c)
        c.Divide(2)

        pad1 = c.cd(2)
        plot.zplot = self.build_zvertex_canvas(pad1)
        pad2 = c.cd(1)
        plot.xyplot = self.build_xyvertex_canvas(pad2)

        return plot

    def pt_phi_dist(self, c=None, opts='COLZ'):
        from ROOT import TCanvas

        if c is None:
            c = TCanvas()

        plot = PlotData(c)
        ptphi = self.analysis.Get('Tracks/pass/PtPhi')
        ptphi.SetStats(0)
        ptphi.SetTitle("%s p_{T} vs Azimuthal Angle Distribution" % self.get_pair_title())
        c.SetLogz()
        ptphi.Draw(opts)
        plot.ptphi = ptphi
        yax = ptphi.GetYaxis()
        yax.SetLabelSize(0.04)
        return plot

    def pt_dist(self, c=None):
        from ROOT import TCanvas

        if c is None:
            c = TCanvas()

        plot = PlotData(c)
        ptphi = self.analysis.Get('Tracks/pass/PtPhi')
        pt = plot.pt = ptphi.ProjectionY("trackpt")
        pt.SetStats(0)
        pt.SetTitle("%s p_{T} Distribution" % self.get_pair_title())
        pt.Draw()

        return plot

    def dedx(self, c=None, **opts):
        from ROOT import TCanvas
        from ROOT import kRed

        color = opts.get('color', kRed)
        title = opts.get('title', "TPC dE/dx ")

        if c is None:
            c = TCanvas()

        plot = PlotData(c)

        track_tdir = self.analysis.Get(f'Tracks')
        dEdX = track_tdir.Get("pass/dEdX")
        dEdX_fail = track_tdir.Get("fail/dEdX")

        plot.p = dEdX
        plot.f = dEdX_fail

        # c.SetLogz()
#         dEdX.Rebin2D(2,2)
        dEdX_fail.Rebin2D(2,2)

        dEdX_fail.SetTitle(title)

        dEdX.SetStats(0)
        dEdX_fail.SetStats(0)

        dEdX_fail.Draw()

        dEdX.Draw(' same')
        dEdX.SetMarkerColor(color)

        c.Draw()

        return plot

    def tof_time(self, c=None, **opts):
        import ROOT
        from ROOT import TCanvas

        fail_color = opts.get('failcolor', ROOT.kBlack)
        color = opts.get('color', ROOT.kRed)

        if c is None:
            c = TCanvas()

        plot = PlotData(c)
        pair_title = self.get_pair_title()

        track_pass_dir =  self.analysis.Get('Tracks/pass')
        track_fail_dir =  self.analysis.Get('Tracks/fail')

        tof_vs_p_pass = plot.p = track_pass_dir.Get("TofVsP")
        tof_vs_p_fail = plot.f = track_fail_dir.Get('TofVsP')

        pion_title = self.get_pair_title()
        tof_vs_p_fail.SetTitle(f"{pion_title} Relative TOF Time")

        tof_vs_p_fail.SetYTitle("TOF Time - TOF(#pi)")
        tof_vs_p_fail.GetYaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleSize(0.038)
        tof_vs_p_fail.SetStats(0)

        tof_vs_p_fail.SetMarkerColor(fail_color)
        tof_vs_p_fail.SetLineColor(fail_color)
        tof_vs_p_pass.SetMarkerColor(color)
        tof_vs_p_pass.SetLineColor(color)

        tof_vs_p_fail.Draw()
        tof_vs_p_pass.Draw("SAME")

        return plot

    def tof_sigma(self, c=None, **opts):
        import ROOT
        from ROOT import gStyle, TCanvas

        fail_color = opts.get('failcolor', ROOT.kBlack)
        color = opts.get('color', ROOT.kRed)

        if c is None:
            c = TCanvas()
        plot = PlotData(c)

        track_dir =  self.analysis.Get('Tracks')
        track_pass_dir = track_dir.Get('pass')
        track_fail_dir = track_dir.Get('fail')

        tof_sigma_fail = track_fail_dir.Get("NsigTof")
        tof_sigma_pass = track_pass_dir.Get("NsigTof")

        tof_sigma_pass.SetLineColor(color)
        tof_sigma_pass.SetMarkerColor(color)
        tof_sigma_fail.SetLineColor(fail_color)
        tof_sigma_fail.SetMarkerColor(fail_color)

        tof_sigma_fail.SetTitle("#sigma")
        tof_sigma_fail.SetStats(0)
#         ("#sigma")

#         tof_vs_p_fail.SetYTitle("TOF Time - TOF(#pi)")
#         tof_vs_p_fail.GetYaxis().SetTitleOffset(1.2)
#         tof_vs_p_fail.GetZaxis().SetTitleOffset(1.2)
#         tof_vs_p_fail.GetZaxis().SetTitleSize(0.038)
#         tof_vs_p_fail.SetStats(0)
#         tof_vs_p_fail.Draw()

        tof_sigma_fail.Draw()
        tof_sigma_pass.Draw('SAME')

        plot.p = tof_sigma_pass
        plot.f = tof_sigma_fail
        return plot

    def pair_detadphi(self, c=None, **opts):
        from ROOT import gStyle, TCanvas

        gStyle.SetPalette(ROOT.kTemperatureMap)

        if c is None:
            c = TCanvas()

        plot = PlotData(c)

        c.Divide(2, 1, 0,0)
        c.cd(1)

        zrng = 0.5, 1.5

        r_nocut.GetZaxis().SetRangeUser(*zrng)
        r_nocut.Draw("COL")
        c.cd(2)
        r.GetZaxis().SetRangeUser(*zrng)
        r.Draw("COLZ")

        return plot

    def pair_detadphi(self, c=None, **opts):
        from ROOT import gStyle, TCanvas

        if c is None:
            c = TCanvas()

        plot = PlotData(c)
        plot
        return plot

    def centralities(self, c=None):
        import ROOT
        from ROOT import gStyle, TCanvas

        if c is None:
            c = TCanvas()

        plot = PlotData(c)
#         for
        plot
        return plot


class ResultPlot1D:

    colors = None

    def __init__(self, fr):
        self.fr = fr

    def draw(self, c=None):
        from ROOT import TGraphErrors, TCanvas

        import seaborn as sns

        if self.colors is None:
            from ROOT import TColor
            self.tcolors = [TColor(*color)
                            for color in sns.color_palette("Paired")]
            self.colors = [tcol.GetNumber() for tcol in self.tcolors]

        if c is None:
            c = TCanvas()

        plot = PlotData(c)

        return plot
