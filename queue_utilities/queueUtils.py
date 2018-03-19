from __future__ import print_function
import re
from datetime import datetime
import subprocess
from time import sleep

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
  Prepare a file to be sent to queue, given a list of commands
  and the required queue configuration file.
  """
  # read config and assert all required config params are given
  config_dict = read_config(config_file)
  required_params = ['RESOURCE','OUT_DIR','QUEUE']
  for param in required_params:
    assert (param in config_dict and config_dict[param]), "Missing or None value for param %s" % param

  #prepare bash file with commands
  commands_list.insert(0,'set -e')
  commands_list.insert(1,"echo JOB START")
  commands_list.append("echo JOB END")
  time_stamp =  str(datetime.now()).replace(' ','_')
  bash_file_path = "%s/%s_%s_commands.sh" % ( config_dict['OUT_DIR'], job_name, time_stamp )
  with open(bash_file_path,'w') as fo:
    print('\n'.join(commands_list), file=fo)

  # prepare queue file
  out_lines = []
  out_lines.append("#!/bin/tcsh")
  out_lines.append("#$ -N %s" % job_name)
  out_lines.append("#$ -S /bin/tcsh")
  out_lines.append("#$ -l %s" % config_dict['RESOURCE'])
  if 'HOSTS' in config_dict:
    out_lines.append("#$ -l h=%s" % config_dict['HOSTS'])
  out_lines.append("#$ -e %s/$JOB_NAME.$JOB_ID.ER" % config_dict['OUT_DIR'])
  out_lines.append("#$ -o %s/$JOB_NAME.$JOB_ID.OU" % config_dict['OUT_DIR'])
  out_lines.append("#$ -q %s" % config_dict['QUEUE'])
  if n_cpu:
    out_lines.append("#$ -pe mpi %s" % n_cpu)
  out_lines.append("bash %s" % bash_file_path)

  # print to job file
  job_file_path = "%s/%s_%s.q" %(config_dict['OUT_DIR'], job_name, time_stamp)
  with open(job_file_path,'w') as fo:
    print('\n'.join(out_lines), file=fo)
  return(job_file_path)

def send_commands_to_queue(job_name, commands_list, config_file, n_cpu = None, verbose = True):
  """
  Send a list of commands to queue as a single job.
  Will block until job is complete and report exit status.
  """
  # prepare required files
  job_file_path = prep_job_for_qsub(job_name, commands_list, config_file, n_cpu)
  # send command and fetch job id
  qsub_command = "qsub %s" % job_file_path
  qsub_res = subprocess.check_output(qsub_command, shell=True).decode('ascii')
  job_id = qsub_res.split()[2]
  # wait until job is done or fails
  exit_status_regex = re.compile(r'exit_status +(\d+)')
  exit_status = None
  while not exit_status:
    sleep(10)
    try:
      job_stat = subprocess.check_output("qacct -j %s" % job_id, shell=True, stderr=subprocess.PIPE).decode('ascii')
    except subprocess.CalledProcessError:
      continue
    res = exit_status_regex.search(job_stat)
    if res:
      exit_status = res.group(1)
      if verbose:
        if exit_status == '0':
          print("Job %s (job id %s) completed successfully" % (job_name, job_id) )
        else:
          print("Job %s (job id %s) failed with exit error %s" % (job_name, job_id, exit_status) )
      return exit_status

if __name__ == "__main__":
  test_commands = ['echo \"test started...\"', 'sleep 30', 'echo \"test complete\"']
  config_file = 'queue.conf'
  job_name = 'queue_utils_test'
  send_commands_to_queue(job_name, test_commands, config_file)
