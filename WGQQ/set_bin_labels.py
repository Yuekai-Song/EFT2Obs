import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--nbins', type=int, default=6)
parser.add_argument('--output', default=None, help='Destination for the bin labels JSON')
args = parser.parse_args()

n_ptg_bins = args.nbins
n_mwp_bins = 5
n_phi_bins = 10

bin_labels = {}
bin_labels['/WGQQ/ptg'] = []
bin_labels['/WGQQ/ptg_ajet'] = []

for j in range(n_ptg_bins):
    bin_labels['/WGQQ/ptg'].append('ptg_{}'.format(j))
    bin_labels['/WGQQ/ptg_ajet'].append('ptg_{}_ajet'.format(j))
for i in range(n_phi_bins):
    bin_labels['/WGQQ/ptg_phi_{}'.format(i)] = []
    bin_labels['/WGQQ/mwp_phi_{}'.format(i)] = []
    bin_labels['/WGQQ/ptg_phi_{}_ajet'.format(i)] = []
    bin_labels['/WGQQ/mwp_phi_{}_ajet'.format(i)] = []
    for j in range(n_ptg_bins):
        bin_labels['/WGQQ/ptg_phi_{}'.format(i)].append('ptg_{}'.format(j) + '_phi_{}'.format(i))
        bin_labels['/WGQQ/ptg_phi_{}_ajet'.format(i)].append('ptg_{}'.format(j) + '_phi_{}_ajet'.format(i))
    for j in range(n_mwp_bins):
        bin_labels['/WGQQ/mwp_phi_{}'.format(i)].append('mwp_{}'.format(j) + '_phi_{}'.format(i))
        bin_labels['/WGQQ/mwp_phi_{}_ajet'.format(i)].append('mwp_{}'.format(j) + '_phi_{}_ajet'.format(i))
with open(args.output or f'./bin_labels_{n_ptg_bins}.json', 'w') as json_file:
    json.dump(bin_labels, json_file, indent=4)
