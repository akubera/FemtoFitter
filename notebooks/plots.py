#
# notebooks/plots.py
#


def make_plots(plots):
    canvases = []
    c = TCanvas()

    c.Divide()

    return c


def random_histnames():
    from random import randint
    while True:
        yield 'h%d' % randint(0, 10000)


histname = random_histnames()

def plot_ratios(n, d, axis, lim=0.01):
    from ROOT import TCanvas, TH3, gStyle

    project = {
        'y': TH3.ProjectionY,
        'x': TH3.ProjectionX,
        'z': TH3.ProjectionZ
    }[axis]

    if isinstance(lim, float):
        l, h = map(n.GetXaxis().FindBin, (-lim, lim))
        limits = (l, h) * 2
    elif len(lim) == 2:
        l, h = map(n.GetXaxis().FindBin, lim)
        limits = (l, h) * 2
    elif len(lim) == 4:
        limits = tuple(map(n.GetXaxis().FindBin, lim))

    opt_entries_only = 10
#     saved_opt_stat = gStyle.GetOptStat()
    gStyle.SetOptStat(opt_entries_only)

    c = TCanvas()
    c.SetCanvasSize(600, 600)
    c.Divide(1, 2)
    pad = c.cd(1)
    pad.UseCurrentStyle()
    pad.Divide(2, 1)
    pad.cd(1)
    np = project(n, "n" + axis, *limits)
    np.Draw("HE")

    pad.cd(2)
    dp = project(d, "d" + axis, *limits)
    dp.Draw("HE")

#     gStyle.SetOptStat(saved_opt_stat)

    c.cd(2)
    rp = np.Clone("r" + axis)
    rp.SetTitle("Left / Right")
    if rp.GetSumw2N() == 0:
        rp.Sumw2()
    rp.Divide(dp)
    rp.SetStats(False)
    rp.DrawCopy("HE")
    c.Draw()
    return c
