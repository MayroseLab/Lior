#!python
import os
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

if "threads" in job_properties:
  ppn = job_properties['threads']
else:
  ppn = 1
if "mem_gb" in job_properties["resources"]:
  mem = ',mem=%sg' % job_properties["resources"]['mem_gb']
else:
  mem = ',mem=1g'
# job name
rule_name = job_properties["rule"]
pref = []
for x in ['sample', 'phenotype', 'genome']:
  if x in job_properties["wildcards"]:
    pref.append(job_properties["wildcards"][x])
if not pref:
  base = "all_samples"
else:
  base = '_'.join(pref)
queue = job_properties['params']["queue"]
priority = job_properties['params']["priority"]
logs_dir = os.path.dirname(job_properties["log"][0])
os.system("qsub -N {base}_{rule} -p {priority} -q {queue} -l nodes=1:ppn={ppn}{mem} -o {logs_dir}/{base}_{rule}.out -e {logs_dir}/{base}_{rule}.err {script}".format(base=base, rule=rule_name, priority=priority, queue=queue, ppn=ppn, mem=mem, logs_dir=logs_dir, script=jobscript))
