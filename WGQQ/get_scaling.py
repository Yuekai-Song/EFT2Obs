import subprocess

n_phi_bins = 10
subprocess.check_call(['rm', '-rf', 'scaling'])
subprocess.check_call(['mkdir', '-p', 'scaling'])
subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'ptg', '-t'])
subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'ptg_ajet', '-t'])
subprocess.check_call(['mv', 'inclusive/WGQQ_ptg.json', 'scaling/ptg.json'])
subprocess.check_call(['mv', 'inclusive/WGQQ_ptg_ajet.json', 'scaling/ptg_ajet.json'])
for i in range(n_phi_bins):
    subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'ptg_phi_{}'.format(i), '-t'])
    subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'mwp_phi_{}'.format(i), '-t'])
    subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'ptg_phi_{}_ajet'.format(i), '-t'])
    subprocess.check_call(['python3', 'makePlot.py', 'inclusive', 'mwp_phi_{}_ajet'.format(i), '-t'])
    subprocess.check_call(['mv', 'inclusive/WGQQ_ptg_phi_{}.json'.format(i), 'scaling/ptg_phi_{}.json'.format(i)])
    subprocess.check_call(['mv', 'inclusive/WGQQ_mwp_phi_{}.json'.format(i), 'scaling/mwp_phi_{}.json'.format(i)])
    subprocess.check_call(['mv', 'inclusive/WGQQ_ptg_phi_{}_ajet.json'.format(i), 'scaling/ptg_phi_{}_ajet.json'.format(i)])
    subprocess.check_call(['mv', 'inclusive/WGQQ_mwp_phi_{}_ajet.json'.format(i), 'scaling/mwp_phi_{}_ajet.json'.format(i)])