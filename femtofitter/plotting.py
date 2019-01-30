#
# femtofitter/plotting.py
#


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

