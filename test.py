import ROOT
import numpy as np
from array import array

# create some example data points
n = 10
x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
y = np.array([3, 4, 2, 7, 6, 9, 5, 1, 8, 0])

np.savetxt("bar_mpvs", np.array([x, y]))

x = array("d",[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
y = array("d",[3, 4, 2, 7, 6, 9, 5, 1, 8, 0])

# create a TGraph object with triangle markers
gr = ROOT.TGraph(n, x, y)
gr.SetMarkerStyle(23)
gr.SetMarkerColor(ROOT.kBlue)

# create a TCanvas and draw the graph
c = ROOT.TCanvas("c", "c", 800, 600)
gr.Draw("AP^") # A stands for axis, P stands for points

# add some axis labels
gr.GetXaxis().SetTitle("x-axis")
gr.GetYaxis().SetTitle("y-axis")


File = ROOT.TFile.Open("test.root", "RECREATE")

gr.Write()

c.Update()
c.Draw()