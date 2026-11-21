from builtins import range
import sys
import argparse
import os
from pathlib import Path

repo_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_dir / 'scripts'))
from jobs import Jobs

os.chdir(repo_dir)
os.environ['PWD'] = str(repo_dir)  # Used by the job script generator.

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
    job_mgr.task_dir = 'jobs/make'
    job_mgr.bopts = 'request_cpus=8\n+MaxRuntime = 28800\ngetenv = False\nuse_x509userproxy = True\nshould_transfer_files = NO'
    os.makedirs(job_mgr.task_dir, exist_ok=True)
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
