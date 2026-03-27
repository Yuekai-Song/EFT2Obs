from __future__ import print_function
from builtins import range
from array import array
import math
import json
import argparse
import yoda
import os
from eftscaling import EFT2ObsHist, EFTScaling

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', default="Rivet.yoda")
parser.add_argument('--dir', '-d', default="WGQQ")
parser.add_argument('--output', '-o', default=None)
parser.add_argument('--config', '-c', default="Rivet.yoda")
parser.add_argument('--hist', default='/HiggsTemplateCrossSectionsStage1/HTXS_stage1_pTjet30')
parser.add_argument('--exclude-rel', default=None, help="Exclude terms with magnitude below this value relative to largest")
parser.add_argument('--rebin', default=None, help="Comma separated list of new bin edges")
parser.add_argument('--save', default='json', help="Comma separated list of output formats (json, txt, latex)")
parser.add_argument('--save-raw', action='store_true', help="Save the raw histogram information as JSON, for further processing")
parser.add_argument('--legacy', action='store_true', help="Use the legacy format for json ouput (if requested)")
parser.add_argument('--translate-tex', default=None, help="json file to translate parameter names to latex")
parser.add_argument('--translate-txt', default=None, help="json file to translate parameter names in the text file")
parser.add_argument('--bin-labels', default=None, help="json file to translate bin labels")
parser.add_argument('--overflow', action='store_true', help="Include overflow binning")
parser.add_argument('--differential', action='store_true', help="Considering differential binning when rebinning")
parser.add_argument('--nlo', action='store_true', help="Set if weights came from NLO reweighting")
parser.add_argument('--filter-params', default=None, help="Specify a subset of parameters to include")
parser.add_argument('--print-style', default="perBin", choices=["perBin", "perTerm"], help="Specify the format for printing to the screen")
parser.add_argument('--color-above', default=None, type=float, help="When using --print-style perTerm, highlight relative uncertainties above this threshold")
args = parser.parse_args()

def sumW(hist, i, overflow=False):
    if 'Estimate' in repr(type(hist)):
        return hist.bins(overflow)[i].val()
    else:
        return hist.bins(overflow)[i].sumW()
def sqr(a):
    return a * a

def sumW2(hist, i, overflow=False):
    if 'Estimate' in repr(type(hist)):
        return sqr(hist.bins(overflow)[i].err('stats')[0])
    else:
        return hist.bins(overflow)[i].sumW2()

# It is super stupid to get number of entries from the raw histo, but no choice so far
def numEntries(hist, i, overflow=False, yoda_dict=None):
    yoda_dict = aos if yoda_dict is None else yoda_dict
    if 'Estimate' in repr(type(hist)):
        hist = yoda_dict[hist.path().replace('/', '/RAW/', 1)]
    return hist.bins(overflow)[i].numEntries()

def scale(obj, sf):
    if hasattr(obj, 'scaleW'):
        obj.scaleW(sf)
    else:
        obj.scale(sf)

def rebinTo(hist, newEdges, differential=True, axis='X'):
    if differential:
        for b in hist.bins():
            scale(b, b.dVol())
    getattr(hist, f'rebin{axis.upper()}To')(newEdges)
    # there's a bug here for Yoda version before 2.1.2, now only local fix
    if differential:
        for b in hist.bins():
            scale(b, 1.0 / b.dVol())

with open(args.config) as jsonfile:
    cfg = json.load(jsonfile)
pars = cfg['parameters']
defs = cfg['parameter_defaults']

if args.output is None:
    auto_name = args.hist
    if auto_name.startswith('/'):
        auto_name = auto_name[1:]
    args.output = auto_name.replace('/', '_')
args.output = os.path.join(args.dir, args.output)
save_formats = args.save.split(',')

translate_tex = {}
if args.translate_tex is not None:
    with open(args.translate_tex) as jsonfile:
        translate_tex = json.load(jsonfile)

translate_txt = {}
if args.translate_txt is not None:
    with open(args.translate_txt) as jsonfile:
        translate_txt = json.load(jsonfile)

bin_labels = list()
if args.bin_labels is not None:
    with open(args.bin_labels) as jsonfile:
        bin_labels = json.load(jsonfile)[args.hist]

hname = args.hist

aos = yoda.read(args.input, asdict=True)

n_pars = len(pars)
n_hists = int(1 + n_pars * 2 + (n_pars * n_pars - n_pars) / 2)

filter = list()
if args.filter_params is not None:
    filter = args.filter_params.split(',')

hists = []
for i in range(n_hists):
    if args.nlo:
        hists.append(aos['%s[rw%.4i_nlo]' % (hname, i)])
    else:
        hists.append(aos['%s[rw%.4i]' % (hname, i)])

# print hists
is2D = isinstance(hists[0], yoda.Histo2D)

if args.rebin is not None and not is2D:
    rebin = [float(X) for X in args.rebin.split(',')]
    for h in hists:
        print(rebin)
        rebinTo(h, rebin, differential=args.differential, axis='X')

# so far, the overflow bin contents of Estimate from rivet are always NaN
overflow = args.overflow and not is2D and not 'Estimate' in repr(type(hists[0]))
nbins = hists[0].numBins(includeOverflows=overflow)
bin_range = range(1, nbins) if overflow else range(nbins)

if is2D:
    edges = [[[hists[0].xMins(overflow)[ib], hists[0].xMaxs(overflow)[ib]], [hists[0].yMins(overflow)[ib], hists[0].yMaxs(overflow)[ib]]] for ib in bin_range]
    # areas = list(hists[0].volumes())
else:
    edges = [[hists[0].xMins(overflow)[ib], hists[0].xMaxs(overflow)[ib]] for ib in bin_range]
areas = [sumW(hists[0], ib, overflow) for ib in bin_range]
    # print (areas,  [hists[0].bins[ib].sumW for ib in range(nbins)])

for p in pars:
    for k in defs:
        if k not in p:
            p[k] = defs[k]

n_divider = 65


def PrintEntry(label, val, err):
    print('%-20s | %12.4f | %12.4f | %12.4f' % (label, val, err, abs(err / val)))


# Generate a list of constants that need to be divided out of each entry
eftconstants = [1.] # for the SM
for ip in range(len(pars)):
    eftconstants.append(pars[ip]['val'])
    eftconstants.append(pars[ip]['val'] * pars[ip]['val'])
for ix in range(0, len(pars)):
    for iy in range(ix + 1, len(pars)):
        eftconstants.append(pars[ix]['val'] * pars[iy]['val'])
assert(len(eftconstants) == len(hists))

for ip, hist in enumerate(hists):
    scale(hist, 1. / eftconstants[ip])


def initTerms(params):
    points = list()
    points.append(list('1'))
    for i in range(len(params)):
        points.append([params[i]])
        points.append([params[i], params[i]])
    for ix in range(0, len(params)):
        for iy in range(ix + 1, len(params)):
            points.append([params[ix], params[iy]])
    return points

e2ohist = EFT2ObsHist(
    terms=initTerms([X['name'] for X in pars]),
    sumW=[[sumW(hist, ib, overflow) for ib in bin_range] for hist in hists],
    sumW2=[[sumW2(hist, ib, overflow) for ib in bin_range] for hist in hists],
    numEntries=[[numEntries(hist, ib, overflow, aos) for ib in bin_range] for hist in hists],
    bin_edges=edges,
    bin_labels=bin_labels)

e2ohist.printToScreen(style=args.print_style, colorAbove=args.color_above)
e2oscaling = EFTScaling.fromEFT2ObsHist(e2ohist, filter=filter)

if args.save_raw:
    print('>> Saving EFT2ObsHist as %s_raw.json' % args.output)
    e2ohist.writeToJSON('%s_raw.json' % args.output)

if 'json' in save_formats:
    print('>> Saving histogram parametrisation to %s.json' % args.output)
    e2oscaling.writeToJSON('%s.json' % args.output, legacy=args.legacy)

if 'yaml' in save_formats:
    print('>> Saving histogram parametrisation to %s.yaml' % args.output)
    e2oscaling.writeToYAML('%s.yaml' % args.output)

if 'txt' in save_formats:
    print('>> Saving histogram parametrisation to %s.txt' % args.output)
    e2oscaling.writeToTxt('%s.txt' % args.output, translate_txt)

if 'tex' in save_formats:
    print('>> Saving histogram parametrisation to %s.tex' % args.output)
    e2oscaling.writeToTex('%s.tex' % args.output, translate_tex)
