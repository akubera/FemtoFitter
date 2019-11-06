

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
        
        if (cent, pair) is (None, None):
            return
        
        if cent is None:
            cent = (clo, chi)
        elif isinstance(cent, str):
            cent = cent.split('_')
        if pair is None:
            pair = p
            
        analysis = '_'.join([*name, *cent, pair])
        self.analysis = self.container.Get(analysis)
    
    def list_analyses(self, path='PWG2FEMTO/*/*'):
        from stumpy.utils import walk_matching
        
        for _, tdir in walk_matching(self.tfile, path):
            yield tdir

    def build_zvertex_canvas(self, c=None):
        from ROOT import gStyle, gROOT, TCanvas
        
        gStyle.SetLabelSize(0.05)
        gROOT.ForceStyle()
        
        if c is None:
            c = TCanvas()

        event_dir = self.analysis.Get("Event/pass")
        plot = PlotData(c)
        
        vz = event_dir.Get("VertexZ")
        vz.SetStats(0)
        
        xax = vz.GetXaxis()
        yax = vz.GetYaxis()
        
        xax.SetTitleOffset(1.4)
        yax.SetTitleOffset(1.4)
        vz.SetTitle("Event Vertex Z-Component")
        vz.Draw()
        plot.vz = vz

        return plot

    def build_xyvertex_canvas(self, c=None):
        
        if c is None:
            c = TCanvas()
            # c.SetCanvasSize(800, 500)

        plot = PlotData(c)

        event_dir = self.analysis.Get("Event/pass")
        vxy = event_dir.Get("VertexXY")
        vxy.SetStats(0)
        
        
        tick_code = 906
        xax = vxy.GetXaxis()
        xax.SetNdivisions(tick_code)
        xax.SetRangeUser(0.034, 0.095)
        yax = vxy.GetYaxis()
        yax.SetNdivisions(tick_code)
        
        c.SetRightMargin(20.0)
        vxy.SetTitle("Event Vertex XY-components")
        vxy.Draw('COL')
        
        plot.vxy = vxy
        
        return plot
    
    def build_vertex_canvas(self, c=None, **opts):
        from ROOT import TCanvas
        if c is None:
            c = TCanvas()
            
            size = opts.get('size', (900, 500))
            c.SetCanvasSize(*size)

        plot = PlotData(c)
        c.Divide(2)
        
        pad1 = c.cd(2)
        plot.zplot = self.build_zvertex_canvas(pad1)
        pad2 = c.cd(1)
        plot.xyplot = self.build_xyvertex_canvas(pad2)
        
        return plot
    
    def pt_dist(self, c=None):
        from ROOT import TCanvas
        
        if c is None:
            c = TCanvas()
        
        plot = PlotData(c)
        ptphi = self.analysis.Get('Tracks/pass/PtPhi')
        pt = plot.pt = ptphi.ProjectionY("trackpt")
        pt.SetStats(0)
        pt.SetTitle("#pi^{+} p_{T} Distribution")
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
        from ROOT import TCanvas
        
        if c is None:
            c = TCanvas()
        
        plot = PlotData(c)
        pair_title = 'pip' in self.analysis
        track_pass_dir =  self.analysis.Get('Tracks/pass')
        track_fail_dir =  self.analysis.Get('Tracks/fail')

        tof_vs_p_pass = plot.p = track_pass_dir.Get("TofVsP")
        tof_vs_p_fail = plot.f = track_fail_dir.Get('TofVsP')
        
        tof_vs_p_fail.SetTitle("#pi^{+} Relative TOF Time")

        tof_vs_p_fail.SetYTitle("TOF Time - TOF(#pi)")
        tof_vs_p_fail.GetYaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleSize(0.038)
        tof_vs_p_fail.SetStats(0)
        tof_vs_p_fail.Draw()

        tof_vs_p_fail.SetMarkerColor(ROOT.kBlack)
        tof_vs_p_fail.SetLineColor(ROOT.kBlack)

        tof_vs_p.SetMarkerColor(ROOT.kRed)
        tof_vs_p.SetLineColor(tof_vs_p.GetMarkerColor())

        tof_vs_p.Draw("SAME")

        return plot
    
    def tof_sigma(self, c=None, **opts):
        from ROOT import gStyle, TCanvas
        
        if c is None:
            c = TCanvas()
        plot = PlotData(c)
            
        track_pass_dir =  analysis.Get('Tracks/pass')
        track_fail_dir =  analysis.Get('Tracks/fail')

        tof_sigma = track_dir.Get("NsigTof")
        tof_vs_p_fail.SetYTitle("TOF Time - TOF(#pi)")
        tof_vs_p_fail.GetYaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleOffset(1.2)
        tof_vs_p_fail.GetZaxis().SetTitleSize(0.038)
        tof_vs_p_fail.SetStats(0)
        tof_vs_p_fail.Draw()
        
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
        