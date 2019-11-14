#
# plotdata.py
#


class PlotData:
    """
    All the plot data
    """

    def __init__(self, canvas=None):
        self.canvas = canvas

    def Draw(self):
        self.canvas.Draw()
