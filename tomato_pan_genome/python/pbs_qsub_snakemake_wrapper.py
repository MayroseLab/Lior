#!python
#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

if "nodes" in job_properties["params"]:
  nodes = job_properties["params"]["nodes"]
else:
  nodes = 1
if "ppn" in job_properties["params"]:
  ppn = job_properties["params"]["ppn"]
else:
  ppn = 1
rule_name = job_properties["rule"]
if 'sample' in job_properties["wildcards"]:
  base = job_properties["wildcards"]['sample']
elif "chunk" in job_properties["wildcards"]:
  base = job_properties["params"]["sample"]
  base += '_' + job_properties["wildcards"]['chunk']
else:
  base = "X"
queue = job_properties["params"]["queue"]
priority = job_properties["params"]["priority"]
logs_dir = job_properties["params"]["logs_dir"]
os.system("qsub -N {base}_{rule} -p {priority} -q {queue} -l nodes={nodes}:ppn={ppn} -o {logs_dir}/{base}_{rule}.out -e {logs_dir}/{base}_{rule}.err {script}".format(base=base, rule=rule_name, priority=priority, queue=queue, nodes=nodes, ppn=ppn, logs_dir=logs_dir, script=jobscript))
