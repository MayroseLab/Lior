#!python
#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

nodes = job_properties["params"]["nodes"]
ppn = job_properties["params"]["ppn"]
rule_name = job_properties["rule"]
sample_name = job_properties["wildcards"]['sample']

os.system("qsub -N {sample}_{rule} -p 0 -q itaym2 -l nodes={nodes}:ppn={ppn} {script}".format(sample=sample_name, rule=rule_name, nodes=nodes, ppn=ppn, script=jobscript))
