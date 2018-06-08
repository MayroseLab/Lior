# Utility functions for interacting with the SGE queue system
from __future__ import print_function
import re
from datetime import datetime
import subprocess
from time import sleep
from IPython.lib import backgroundjobs as bg
from IPython import get_ipython
import os

def read_config(config_path):
  """
  read a config file and return the respective dict.
  Configuration lines expected format is: <key> = <value>
  """
  conf_dict = {}
  spaces_regex = re.compile(r'^[ \t]*\n?$')
  valid_regex = re.compile(r'^([^=\s]+)\s*=\s*(.+)$')
  n = 0
  with open(config_path) as f:
    for line in f:
      n += 1
      if line.startswith('#') or spaces_regex.match(line):
        continue
      line = line.strip()
      res = valid_regex.search(line)
      if not res:
        print("WARNING: Invalid configuration line %s (line #%s) in file %s. Skipping line" % (line,n,config_path) )
        continue
      key, val = res.groups()
      conf_dict[key] = val
  return conf_dict

def prep_job_for_qsub(job_name, commands_list, config_file, n_cpu = None):
  """
  Prepare a file to be sent to queue, given a list of commands.
  """
  # read config and assert all required config params are given
  config_dict = read_config(config_file)
  required_params = ['RESOURCE','OUT_DIR','QUEUE']
  for param in required_params:
    assert (param in config_dict and config_dict[param]), "Missing or None value for param %s" % param

  #prepare q file
  # SGE-related commands
  sge_commands = []
  sge_commands.append("#!/bin/bash")
  sge_commands.append("#$ -N %s" % job_name)
  sge_commands.append("#$ -S /bin/bash")
  sge_commands.append("#$ -l %s" % config_dict['RESOURCE'])
  if 'HOSTS' in config_dict:
    sge_commands.append("#$ -l h=%s" % config_dict['HOSTS'])
  sge_commands.append("#$ -e %s/$JOB_NAME.$JOB_ID.ER" % config_dict['OUT_DIR'])
  sge_commands.append("#$ -o %s/$JOB_NAME.$JOB_ID.OU" % config_dict['OUT_DIR'])
  sge_commands.append("#$ -q %s" % config_dict['QUEUE'])
  if n_cpu:
    sge_commands.append("#$ -pe mpi %s" % n_cpu)
  # job commands
  commands_list.insert(0,'set -e')
  commands_list.insert(1,"echo JOB START")
  commands_list.append("echo JOB END")
  # print commands
  all_commands = sge_commands + commands_list
  time_stamp =  str(datetime.now()).replace(' ','_')
  q_file_path = "%s/%s_%s_commands.q" % ( config_dict['OUT_DIR'], job_name, time_stamp )
  with open(q_file_path,'w') as fo:
    print('\n'.join(all_commands), file=fo)
  return(q_file_path)

def get_job_exit_code(job_id):
  """
  Returns exit status (int) of job id.
  Returns None if job is still running.
  Returns False if job id is not found.
  """
  # check if job is currently running
  try:
    qstat_j_res = subprocess.check_output("qstat -j %s" % job_id, shell=True, stderr=subprocess.PIPE).decode('ascii')
    return None
  except:
    # if it's not running, check if can get exit code
    try:
      job_stat = subprocess.check_output("qacct -j %s" % job_id, shell=True, stderr=subprocess.PIPE).decode('ascii')
      exit_status_regex = re.compile(r'exit_status +(\d+)')
      res = exit_status_regex.search(job_stat)
      if res:
        exit_status = res.group(1)
        return int(exit_status)
    # if not, it means that job id doesn't exist
    except subprocess.CalledProcessError:
      return False

def send_commands_to_queue(job_name, commands_list, config_file, n_cpu = None, block = True, verbose = True):
  """
  Send a list of commands to queue as a single job.
  Will block until job is complete and report exit status.
  """
  # convert single command to list, if given as string
  if type(commands_list) == str:
    commands_list = [commands_list]
  # prepare required files
  job_file_path = prep_job_for_qsub(job_name, commands_list, config_file, n_cpu)
  # send command and fetch job id
  qsub_command = "qsub %s" % job_file_path
  qsub_res = subprocess.check_output(qsub_command, shell=True).decode('ascii')
  job_id = qsub_res.split()[2]
  exit_status = None
  if not block:
    if verbose:
      print("Job %s (job id %s) sent to queue" % (job_name, job_id) )
    return job_id, exit_status
  # wait until job is done or fails
  exit_status = get_job_exit_code(job_id)
  while exit_status is None:
    sleep(10)
    exit_status = get_job_exit_code(job_id)
  sleep(10)	# this help avoid losing jobs in the 'twilight zone'
  exit_status = get_job_exit_code(job_id)
  if verbose:
    if type(exit_status) is int and exit_status == 0:
      print("Job %s (job id %s) completed successfully" % (job_name, job_id) )
    elif type(exit_status) is int:
      print("Job %s (job id %s) failed with exit code %s" % (job_name, job_id, exit_status) )
    else:
      print("Can't get status for Job %s (job id %s). Something went wrong!" % (job_name, job_id))
  return job_id, exit_status

def send_jobs_to_queue(jobs_dict, config_file, n_cpu = None):
  """
  Send a batch of jobs to queue, each with its own name and commands list.
  Input is a dict with the format {job_id: [commands]}
  Does not block and returns a dict {job_name: job id}
  """
  res_dict = {}
  for job_name, commands in jobs_dict.items():
    j = send_commands_to_queue(job_name, commands, config_file, n_cpu, block = False, verbose = False)
    job_id = j[0]
    res_dict[job_name] = job_id
  return res_dict

def run_script_in_bg(script_path,*args):
  """
  Not actually an interaction with the queue, but useful for running
  python scripts that send commands to queue (such as pipelines) without
  blocking the notebook.
  """
  def run_script(sp,args_tup):
    magic_statement = 'run ' + sp + ' ' + ' '.join(args_tup)
    get_ipython().magic(magic_statement)
  jobs = bg.BackgroundJobManager()
  jobs.new(run_script,script_path,args)

if __name__ == "__main__":
  test_commands = ['echo \"test started...\"', 'sleep 30', 'echo \"test complete\"']
  config_file = os.path.dirname(os.path.realpath(__file__)) + '/queue.conf'
  job_name = 'queue_utils_test'
  send_commands_to_queue(job_name, test_commands, config_file)
