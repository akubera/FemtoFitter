#
# femtofitter/plotting/systematics.py
#

from .plotdata import PlotData


def plot_systematics(df, c=None):
    pass


def plot_systematics_plot(df, c=None):
    from ROOT import TH1D, TLine, TCanvas

    if c is None:
        c = TCanvas()
    
    c.Divide(2, 2)
    c.SetCanvasSize(800, 800)

    plot = PlotData(c)

    # point_values = list(points.values())
    hist_axis = plot.hist_axis = TH1D("_hist_axis", "R_{out}; k_{T} (GeV)", 2, 0.2, 1.2)
    hist_axis.GetYaxis().SetRangeUser(2.0, 7.0)
    hist_axis.SetStats(0)
    
    hists = plot.hists = [
        TH1D("_hist_axis", "R_{out}; k_{T} (GeV)", 2, 0.2, 1.2),
        TH1D("_hist_axis", "R_{side}; k_{T} (GeV)", 2, 0.2, 1.2),
        TH1D("_hist_axis", "R_{long}; k_{T} (GeV)", 2, 0.2, 1.2),
    ]
    
    for h in hists:
        h.SetStats(0)
    #     h.SetTitleSize(0.07)
        h.GetYaxis().SetRangeUser(1.3, 8.0)
    
    lines = plot.lines = []
    # pt_itr = (iter(point_values), ) * 2
    
    leg_cent = TLegend(0.5, 0.3, 0.745, 0.5)
    leg_cent.SetHeader("Centrality", "C")
    
    leg_temp = TLegend(0.75, 0.35, 0.95, 0.5)
    # leg_temp.SetHeader("Particlization Temperatures", "C")
    leg_temp.SetHeader("T_{Freezout}", "C")
    
    lines.append(TLine())
    lines[-1].SetLineWidth(2)
    leg_temp.AddEntry(lines[-1], '165 MeV', 'L')
    lines.append(TLine())
    lines[-1].SetLineWidth(2)
    lines[-1].SetLineStyle(2)
    leg_temp.AddEntry(lines[-1], '156 MeV', 'L')
    leg_temp.SetBorderSize(0)
    
    for i, (h, R_key) in enumerate(zip(hists, ("Ro", "Rs", "Rl")), 1):
        c.cd(i)
        h.Draw()
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['00_05']['thi'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(colors[0])
        if i == 1:
            leg_cent.AddEntry(lines[-1], "0-5%", 'L')
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['20_30']['thi'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(colors[1])
        if i == 1:
            leg_cent.AddEntry(lines[-1], "20-30%", 'L')
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['40_50']['thi'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(colors[2])
        if i == 1:
            leg_cent.AddEntry(lines[-1], "40-50%", 'L')
    
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['00_05']['tlo'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineStyle(2)
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(2)
            lines[-1].SetLineColor(colors[0])
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['20_30']['tlo'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineStyle(2)
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(2)
            lines[-1].SetLineColor(colors[1])
    
        for p1, p2 in sliding_window(2, zip(*theory_data[R_key]['40_50']['tlo'].values())):
            lines.append(TLine(*p1, *p2))
            lines[-1].Draw()
            lines[-1].SetLineStyle(2)
            lines[-1].SetLineWidth(2)
            lines[-1].SetLineColor(2)
            lines[-1].SetLineColor(colors[2])
    
    for i, key in enumerate(("Ro", "Rs", "Rl"), 1):
        c.cd(i)
        for cidx, (cent, (graph, syst)) in enumerate(hist_dict[key].items()):
            color = {'00_10': ROOT.kBlue+1,
                     '20_30': ROOT.kOrange-4,
                     '40_50': ROOT.kGreen+2}[cent]
            color = colors[cidx]
            graph.SetMarkerStyle(21)
            graph.SetMarkerSize(0.8)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            
            syst.SetFillColor(syscolors[cidx])
            syst.SetFillStyle(1001)
            syst.SetLineColor(ROOT.kBlack)
    #         syst.SetLineColor(ROOT.kBlack)
            syst.Draw("SAME P E5")
            graph.Draw('SAME P')
            
    c.cd(0)
    leg_cent.Draw()
    leg_temp.Draw()
    c.Draw()
    
    return plot
