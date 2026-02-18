import subprocess
import argparse
import os
import glob

def merge(dir):
    yoda_file = os.path.join(dir, "Rivet.yoda")
    yoda_files = glob.glob(os.path.join(dir, "Rivet_*.yoda"))
    if not os.path.isfile(yoda_file):
        subprocess.check_call(['yodamerge', '-o', yoda_file] + yoda_files)
    return yoda_file

os.chdir('/grid_mnt/data__data.polcms/cms/song/CMSSW_14_1_0_pre4/src/EFT2Obs')
x_titles = {'p_pt': 'p_{T}(\\gamma) [GeV]', 'ptg': 'p_{T}(\\gamma) [GeV]', 'ptg_ajet': 'p_{T}(\\gamma) [GeV]', 'phi_sjet': '\\phi', 'phi_quark': '\\phi_{quark}', 
            'phi_sjet_matched': '\\phi', 'phi_quark_matched': '\\phi_{quark}', 'phi_sjet_noajet': '\\phi', 'phi_sjet_ajet': '\\phi',
            'w_pt_matched': 'p_{T}(W) [GeV]', 'w_pt': 'p_{T}(W) [GeV]', 'wp_mass': 'M_{W\\gamma} [GeV]'}
ranges = {'p_pt': '0.5,1.5', 'ptg_ajet': '0.5,1.5', 'ptg': '0.5,1.5', 'phi_sjet': '0.5,1.5', 'phi_quark': '0.0,2.0', 
          'phi_sjet_matched': '0.0,5.0', 'phi_quark_matched': '0.0,5.0', 'phi_sjet_noajet': '0.5,1.5',
          "w_pt_matched": '0.0,5.0', 'w_pt': '0.0,2.0', 'wp_mass': '0.0,5.0'}
cw_vals = ['cwtil=0.5:2', 'cw=0.5:4']

parser = argparse.ArgumentParser()
parser.add_argument('dir', type=str, default='inclusive',
                    help='the directory with the yoda files')
parser.add_argument('var', type=str, default='sjet_phi',
                    help='the variable to plot')
parser.add_argument('--browser', '-b', action='store_true')
parser.add_argument('--square', '-s', action='store_true')
parser.add_argument('--cross', '-c', action='store_true')
parser.add_argument('--to-scaling', '-t', default=False, action='store_true')
parser.add_argument('--title', '-title', type=str, default='')
args = parser.parse_args()

dir = os.path.join("./WGQQ", args.dir)
cfg = os.path.join('./cards', 'WGQQ', 'config.json')
var = args.var
pta_bins = ['200', '400', '600', '800', '1000']
if 'inclusive' in dir or args.dir == "ptg400":
    yoda_file = os.path.join(dir, "Rivet.yoda")
    yoda_files = []
    if not os.path.isfile(yoda_file):
        subprocess.check_call(['mkdir', '-p', dir])
        for i in range(len(pta_bins)):
            if i < len(pta_bins) - 1:
                subdir = os.path.join("./WGQQ", os.path.dirname(args.dir), 'ptg' + pta_bins[i] + 'to' + pta_bins[i + 1])
            else:
                subdir = os.path.join("./WGQQ", os.path.dirname(args.dir), 'ptg' + pta_bins[i])
            print(subdir)
            yoda_files.append(merge(subdir))
        subprocess.check_call(['yodastack', '-o', yoda_file] + yoda_files)
else:
    yoda_file = merge(dir)

rivet = 'WGQQ'
if args.to_scaling:
    subprocess.check_call(['python3', 'scripts/get_scaling.py', '-c', cfg, '-i', yoda_file, '--hist', '/' + rivet + '/' +
                      var, '--save', 'json,txt,tex', '--translate-tex', 'resources/translate_tex.json', '--bin-labels', 'WGQQ/bin_labels.json', '--dir', dir])

else:
    subprocess.check_call(['python3', 'scripts/get_scaling.py', '-c', cfg, '-i', yoda_file, '--hist', '/' + rivet + '/' +
                    var, '--save', 'json,txt,tex', '--translate-tex', 'resources/translate_tex.json', '--dir', dir])
    if args.var not in ranges:
        range_v = '0.5,1.5'
    else:
        range_v = ranges[args.var]
    if args.var not in x_titles:
        if 'ptg' in var:
            x_titles[var] = 'p_{T}(\\gamma) [GeV]'
        if 'mwg' in var:
            x_titles[var] = 'M_{W\\gamma} [GeV]'
    plot_args = []
    if not args.square:
        plot_args.append('--no-square')
    if not args.cross:
        plot_args.append('--no-cross')
    # plot_args.append('--logy')

    subprocess.check_call(['python3', 'scripts/makePlot.py', '--hist', os.path.join(dir, rivet + '_' + var + '.json'), '-c', cfg, '--x-title', x_titles[var], '--title-left',
                          'W^{\\pm}\\gamma \\rightarrow qq\\gamma', '--title-right', args.title, '--ratio', range_v, '--draw'] + cw_vals + ['--show-unc', '--y-min', '1E-9', '--translate', 'resources/translate_root.json'] + plot_args)

    if args.browser:
        subprocess.check_call(['mkdir', '-p', os.path.join('plots_browser/', args.dir)])

        subprocess.check_call(['cp', os.path.join(
            dir, rivet + '_' + var + '.pdf'), os.path.join('plots_browser/', args.dir)])
        plots = glob.glob(os.path.join(dir, 'WGQQ_' + var + '.p*'))
        subprocess.check_call(['cp'] + plots + [os.path.join('plots_browser/', dir)])
