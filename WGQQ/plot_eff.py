from __future__ import print_function
from builtins import range
import subprocess
import argparse
import ROOT
import glob
import sys
import os
os.chdir('/grid_mnt/data__data.polcms/cms/song/EFT2Obs')
sys.path.append(os.path.join('.', 'scripts'))

import plotting as plot
from array import array
from eftscaling import EFTScaling

x_titles = {'p_pt': 'p_{T}(\\gamma) [GeV]', 'phi_sjet': '\\phi', 'phi_quark': '\\phi_{quark}', 'w_pt': 'p_{T}(W) [GeV]'}

parser = argparse.ArgumentParser()
parser.add_argument('dir', type=str, default='ptg200',
                    help='the directory with the yoda files')
parser.add_argument('var', type=str, default='sjet_phi',
                    help='the variable to plot')
parser.add_argument('--browser', '-b', action='store_true')
args = parser.parse_args()

dir = os.path.join("./WGQQ", args.dir)
var = args.var
json1 = os.path.join(dir, 'WGQQ_' + var + '.json')
json2 = os.path.join(dir, 'WGQQ_' + var + '_matched.json')

if not os.path.isfile(json1) or not os.path.isfile(json2):
    cfg = os.path.join('./cards', 'WGQQ-' + args.dir, 'config.json')
    yoda_file = os.path.join(dir, "Rivet.yoda")
    if not os.path.isfile(yoda_file):
        yoda_files = glob.glob(os.path.join(dir, "Rivet_*.yoda"))
        subprocess.check_call(['yodamerge', '-o', yoda_file] + yoda_files)
    subprocess.check_call(['python', 'scripts/get_scaling.py', '-c', cfg, '-i', yoda_file, '--hist', '/WGQQ/' +
                        var, '--save', 'json,txt,tex', '--translate-tex', 'resources/translate_tex.json', '--dir', dir])
    subprocess.check_call(['python', 'scripts/get_scaling.py', '-c', cfg, '-i', yoda_file, '--hist', '/WGQQ/' +
                        var + "_matched", '--save', 'json,txt,tex', '--translate-tex', 'resources/translate_tex.json', '--dir', dir])

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

x_title = x_titles[var]

jhist1 = EFTScaling.fromJSON(json1)
jhist2 = EFTScaling.fromJSON(json2)
canv = ROOT.TCanvas("matching_with_" + var, '')
pads = plot.TwoPadSplit(0.30, 0.01, 0.02)

h_inc = jhist1.getNominalTH1('nominal')
h_matched = jhist2.getNominalTH1('nominal')
h_axes = [h_inc.Clone() for x in pads]
for h in h_axes:
    h.Reset()
h_axes[1].GetXaxis().SetTitle(x_title)
h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()

legend = ROOT.TLegend(0.60, 0.88 - 0.05 * 2, 0.90, 0.91, '', 'NBNDC')
legend.AddEntry(h_inc, 'inclusive', 'L')
legend.AddEntry(h_matched, 'matched', 'L')

plot.Set(h_inc, LineColor=1, LineWidth=3)
plot.Set(h_matched, LineColor=2, LineWidth=3)
h_inc.Scale(1., 'width')
h_matched.Scale(1., 'width')
h_inc.Draw('HISTSAME')
h_matched.Draw('HISTSAME')
plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.30)
legend.Draw()

pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

h_eff = plot.MakeRatioHist(h_matched, h_inc, False, False)
h_eff.Draw("HISTSAME")

plot.SetupTwoPadSplitAsRatio(pads, plot.GetAxisHist(pads[0]), plot.GetAxisHist(pads[1]), 'Eff.', True, 0, 1)
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

title_left = 'W\\gamma \\rightarrow qq\\gamma'
title_right = 'SMEFTsim3'
plot.DrawTitle(pads[0], title_left, 1.5)
plot.DrawTitle(pads[0], title_right, 3)

if args.browser:
    subprocess.check_call(['mkdir', '-p', os.path.join('plots_browser/', args.dir)])
    canv.Print(os.path.join('./plots_browser', args.dir, "matching_with_" + var + ".pdf"))
else:
    canv.Print(os.path.join(dir, "matching_with_" + var + ".pdf"))