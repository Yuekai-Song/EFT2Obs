from builtins import range
import sys
sys.path.append('./scripts')
from jobs import Jobs
import argparse
import os
os.chdir('/grid_mnt/data__data.polcms/cms/song/EFT2Obs')

job_mgr = Jobs()
parser = argparse.ArgumentParser()
parser.add_argument('--step', default='make', choices=['make', 'run'], help="which step for the gridpack")   
parser.add_argument('-p', '--pt-bins', type=str, default='200.0,1000.0', help="Initial random number seed")

job_mgr.attach_job_args(parser)
(args, unknown) = parser.parse_known_args()
job_mgr.set_args(args)

if args.step == 'make':
    job_mgr.job_mode = 'condor'
    job_mgr.task_name = 'divide'
    job_mgr.task_dir = 'jobs'
    job_mgr.bopts = '+JobFlavour = "longlunch"\nT3Queue = long\ngetenv = False\nuse_x509userproxy = True\nWNTag=el9\n+SingularityCmd = ""\ninclude : /opt/exp_soft/cms/t3_tst/t3queue |'
else:
    job_mgr.job_mode = 'interactive'
    job_mgr.task_name = 'run'

iwd = os.getcwd()
pt_bins = args.pt_bins.split(',')
pt_bins.append('-1.0')
print(iwd)
for i in range(len(pt_bins) - 1):
    cmd = []
    pt_min = pt_bins[i]
    pt_max = pt_bins[i + 1]
    if args.step == 'make':
        gp_cmd = '%s/WGQQ/make_gridpack.sh %s %s' % (iwd, pt_min, pt_max)
    else:
        gp_cmd = '%s/WGQQ/run_gridpack.sh %s %s' % (iwd, pt_min, pt_max)
    if len(unknown) > 0:
        gp_cmd += ' %s' % (' '.join(unknown))
    cmd.append(gp_cmd)
    job_mgr.job_queue.append('; '.join(cmd))

job_mgr.flush_queue()
