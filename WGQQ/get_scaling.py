import subprocess
import os
from pathlib import Path
from storage import output_directory
import argparse
import math


parser = argparse.ArgumentParser()
parser.add_argument('--rebin', type=str, default='200,350,500,700,1000,2000',
                    help='Comma separated ptg bin edges (must match existing histogram edges)')
parser.add_argument('--output', type=str, default=None,
                    help='Output directory under WGQQ/output (default: scaling_<bin edges>)')
args = parser.parse_args()
try:
    pt_edges = [float(edge) for edge in args.rebin.split(',')]
except ValueError:
    parser.error('--rebin must contain comma separated numbers')
if len(pt_edges) < 2 or not all(math.isfinite(edge) for edge in pt_edges) or any(
        left >= right for left, right in zip(pt_edges, pt_edges[1:])):
    parser.error('--rebin needs at least two finite, strictly increasing edges')

os.chdir(Path(__file__).resolve().parent.parent)
output_dir = Path(output_directory())
bin_tag = '_'.join(format(edge, 'g') for edge in pt_edges)
scaling_dir = output_dir / (args.output or 'scaling_' + bin_tag)
inclusive_dir = output_dir / 'inclusive'

n_phi_bins = 10
pt_rebin = args.rebin

subprocess.check_call(['mkdir', '-p', str(scaling_dir)])
bin_labels = str(scaling_dir / 'bin_labels.json')
subprocess.check_call(['python3', 'WGQQ/set_bin_labels.py', '--nbins', str(len(pt_edges) - 1),
                       '--output', bin_labels])
subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'ptg', '-t', '--bin-labels', bin_labels, '--rebin', pt_rebin])
subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'ptg_ajet', '-t', '--bin-labels', bin_labels, '--rebin', pt_rebin])
subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_ptg.json'), str(scaling_dir / 'ptg.json')])
subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_ptg_ajet.json'), str(scaling_dir / 'ptg_ajet.json')])
for i in range(n_phi_bins):
    subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'ptg_phi_{}'.format(i), '-t', '--bin-labels', bin_labels, '--rebin', pt_rebin])
    subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'mwp_phi_{}'.format(i), '-t', '--bin-labels', bin_labels])
    subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'ptg_phi_{}_ajet'.format(i), '-t', '--bin-labels', bin_labels, '--rebin', pt_rebin])
    subprocess.check_call(['python3', 'WGQQ/makePlot.py', 'inclusive', 'mwp_phi_{}_ajet'.format(i), '-t', '--bin-labels', bin_labels])
    subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_ptg_phi_{}.json'.format(i)), str(scaling_dir / 'ptg_phi_{}.json'.format(i))])
    subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_mwp_phi_{}.json'.format(i)), str(scaling_dir / 'mwp_phi_{}.json'.format(i))])
    subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_ptg_phi_{}_ajet.json'.format(i)), str(scaling_dir / 'ptg_phi_{}_ajet.json'.format(i))])
    subprocess.check_call(['mv', str(inclusive_dir / 'WGQQ_mwp_phi_{}_ajet.json'.format(i)), str(scaling_dir / 'mwp_phi_{}_ajet.json'.format(i))])
