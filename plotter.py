#
# plotter.py
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from stumpy.utils import walk_matching
from statistics import mean

from pathlib import Path
plt.rcParams["legend.numpoints"] = 1


class FitResults:

    def __init__(self, path):
        """ Get results from a file """
        path = Path(path).resolve()
        if path.suffix == ".json":
            import json
            with path.open() as p:
                df = pd.DataFrame(json.load(p))
        else:
            raise TypeError("Uknown type of file %r" % (str(path)))

        # ensure index is integer (not string)
        df.index = pd.Index(np.array(df.index, dtype=int))
        # interpret kt_bin as 'kt'
        df['kt_bin'] = df.kt
        df['kt'] = df.kt_bin.map(lambda x: mean(map(float, x.split('_'))))
        self.df = df.sort_index()

    @classmethod
    def FromBytes(cls, data):
        from io import BytesIO
        buff = BytesIO(data)
        self = cls.__new__(cls)
        self.df = pd.read_pickle(buff)
        return self

    @staticmethod
    def make_quadplot(df, groups=('pair_type', ), limit_alpha=True, cfg=None):
        """
        Make and return the plot
        """
        default_plot_ops = {"linestyle": "-", 'marker': 'o', 'xlim': (0.2, 1.0)}
        plt.close()
        if cfg is not None:
            df = df[df.cfg_hash == cfg]

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
                    cda = np.array(cd)
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
            result = []
            for g in groups:
                for group_val, pair_data in df.groupby(g):
                    fig = _do_makeplot(pair_data)
                    result.append((group_val, fig))
        else:
            result = _do_makeplot(df)

        return result

    @classmethod
    def show_quadplot(cls, df, groups=()):
        cls.make_quadplot(df, groups=groups)
        plt.show()

    def show_root_quadplot(self):
        from ROOT import TGraphErrors

        for cent, cent_data in self.df.groupby('centrality'):
            TGraphErrors("")

    def make_resid_plot(self, key, cent='00_05', kt_bin='0.2_0.3'):
        if not key.startswith('resid'):
            key = 'resid_%s' % key

        if isinstance(cent, (list, tuple)):
            for c in cent:
                self.make_resid_plot(key, c, kt_bin)
            return

        if isinstance(kt_bin, (list, tuple)):
            for k in kt_bin:
                self.make_resid_plot(key, cent, k)
            return

        query = (self.df.kt_bin == kt_bin) & (self.df.centrality == cent)
        data = np.array(self.df[query][key][0])
        plt.plot(data[:,0], data[:,1], marker='.', label='%s/%s' % (cent, kt_bin))


class BigProjectionPlot:
    """
    """

    def __init__(self, filename, cfg='cfgNONE', centrality="00_05", limit=0.04):
        import ROOT
        from ROOT import TFile
        from ROOT import TCanvas
        self.limit = limit
        if isinstance(filename, str):
            self.file = TFile.Open(filename)
        else:
            self.file = filename

        if not self.file:
            raise ValueError("Could not use %r" % filename)

        c = TCanvas()
        c.SetCanvasSize(900, 900)
        J = 6
        c.Divide(3, J, 0, 0)
        self._histcache = []
        for j, (ktbin, container) in enumerate(walk_matching(self.file, "%s/pip/%s/*" % (cfg, centrality))):
            print(ktbin)
            num, fit, den = map(container.Get, ("neg_num", "neg_fit_num", "neg_den"))
            #num, fit, den = map(container.Get, ("avg_num", "avg_fit_num", "avg_den"))
            data = self.projections_ratio(num, den)
            fit_data = self.projections_ratio(fit, den)
            hist_row = list(zip(data, fit_data))
            for i, (p, f) in enumerate(hist_row):
                f.SetLineColor(ROOT.kRed)
                # f.GetPainter().SetDrawOption("HIST C")

                idx = j*3+i+1
                c.cd(idx)
                #hs = ROOT.THStack(self.random_name(), "");
                #hs.Add(p)
                #hs.Add(f)
                #hs.Draw("nostack")
                #hs.SetMinimum(0.95)
                #self._histcache.append(hs)
                p.DrawCopy()
                f.DrawCopy("HIST SAME C")
            if idx >= 3*J:
                break
        self.canvas = c

    @staticmethod
    def random_name(N=8, prefix=''):
        import random
        from string import ascii_lowercase
        return prefix + ''.join(random.choice(ascii_lowercase) for _ in range(N))

    def projections(self, hist):
        from ROOT import TH3
        projs = TH3.ProjectionX, TH3.ProjectionY, TH3.ProjectionZ
        limit = self.limit
        limits = list(map(hist.GetXaxis().FindBin, (-limit, limit))) * 2
        # print("limits", limits)
        for proj in projs:
            yield proj(hist, self.random_name(), *limits)

    def projections_ratio(self, num, den):
        for n, d in zip(self.projections(num), self.projections(den)):
            n.Divide(d)
            n.SetStats(False)
            n.SetTitle("")
            # n.GetYaxis().SetRangeUser(0.96, 1.28)
            yield n
