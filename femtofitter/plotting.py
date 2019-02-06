#
# femtofitter/plotting.py
#

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



class QuadPlot:

    def __init__(self, fr):
        from . import FitResults
        if isinstance(fr, FitResults):
            self.df = fr.df
        else:
            self.df = fr


    def show(self, *groups, limit_alpha=True, cfg=None):
        import matplotlib.pyplot as plt

        if cfg is not None:
            df = self.df[self.df.cfg == cfg]
        else:
            df = self.df

        default_plot_ops = {
            "linestyle": "-",
            'marker': 'o',
            'xlim': (0.2, 1.0),
        }

        plt.close()

        def _do_makeplot(data):
            fig, axs = plt.subplots(2, 2, figsize=(12, 12))
            for cent, cent_data in data.groupby("cent"):
                #print(cent, cent_data)
                cent = cent.replace('_', '-') + "%"
                for ax, key in zip(axs.flat, ("Ro", "Rs", "Rl", 'lam')):
                    plot_ops = default_plot_ops.copy()
                    if limit_alpha and key == 'alpha':
                        plot_ops['ylim'] = (1.0, 2.2) if limit_alpha is True else limit_alpha

                    # plot_ops['yerr'] = key + "_err"

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

                    # cent_data.sort_values("kT").plot("kT", key, ax=ax, title=title, label=cent, **plot_ops)
                    if key == 'Ro':
                        cd = [(kT, kTdata[key].min()) for kT, kTdata in cent_data.sort_values("kT").groupby('kT')]
                    else:
                        cd = [(kT, kTdata[key].mean()) for kT, kTdata in cent_data.sort_values("kT").groupby('kT')]
                    cda = pd.np.array(cd)
                            # for _, dat in groupby]
                    plot_ops.pop("xlim")
                    ax.errorbar(*cda.T, label=cent, **plot_ops)
                    ax.plot(*cda.T, label=cent, **plot_ops)
                    ax.set_title(title)

                    # if key != 'lam':
                    #     # ax.legend_.set_visible(False)
                    # else:
                    #     leg = ax.legend(numpoints=1, loc='best', fontsize=16)
                    #     if leg:
                    #         leg.set_title("Centrality", prop={"size": 16, 'weight': 'bold'})
                    # if key.startswith("R"):
                    #     ax.set_ylim(0.0, 8.0)
                    # else:
                    #     ax.set_ylim(0.2, 0.8)
            return fig

        if groups:
            result = [(group_val, _do_makeplot(pair_data))
                      for group_val, pair_data in map(self.df.groupby, groups)]
        else:
            result = _do_makeplot(df)

        return result





        plt.show()
