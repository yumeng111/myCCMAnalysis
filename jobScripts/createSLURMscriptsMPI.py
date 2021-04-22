#! /usr/bin/python3
# To run you must pass it a configuration file in the .ini format
# An example of the configuration file is shown in $CCMINSTALL/jobScripts/*.ini
# 
# The following is a copy of one of those files (exampleCreateSLURMConfig.ini) for completeness sake
#
# ; I defined the following enviromental variabels in my .bash_profile to the lines do not have to be so long below
# ; $CCMPROCESSED = /lustre/scratch3/turquoise/rtthorn/processedFiles
# ; $CCMRAW = /lustre/scratch3/turquoise/rtthorn/rawFiles
# ; $CCMINSTALL = /users/rtthorn/CCMCode/CCM_analysis_sw_dev_xmlConfig
#
# ; These are the defults. You only need to include them below if you want to change them
# [DEFAULT]
# Rewrite = false
# Submit = false
# MaxNumJobs = 300
# MultiFileRegex = 
# JobNamePrefix = findEvents
# InputEndString = _physics
# NewFilePostfix = events
# SourceFile = $HOME/.bash_profile
# InputDataRegex = *.root
#
# ; Have to include: InputDir, OutputDir, LogDir, ScriptDir, SLURMScript, XMLScript
# ; Can include any of the ones that are listed under DEFAULT to change them from their
# ; default values
# [Params]
# InputDir = $CCMRAW
# OutputDir = $CCMPROCESSED/beam/test
# LogDir = $CCMPROCESSED/beam/test/log
# ScriptDir = $CCMPROCESSED/beam/test/scripts
# SLURMScript = $CCMINSTALL/jobScripts/findEvents_slrum.script
# XMLScript = $CCMINSTALL/configFiles/findEventsConfig.xml
# InputDataRegex = *run000188-*.root

import sys, os
import subprocess
from subprocess import PIPE
import re
import time
import glob
import datetime
import argparse
import sys
import codecs
import configparser

def unescaped_str(arg_str):
    return codecs.decode(str(arg_str), 'unicode_escape')

#Creates slurm and xml config file for the find events module. Submits the slurm script if submit = true.
config = configparser.ConfigParser()
config.read(sys.argv[1])
options = config['Params']

# -------- default settings ------------------------
#Directory location of the input files
input_loc      = os.path.expandvars(options['InputDir'])
#Directory location of the output files
output_loc     = os.path.expandvars(options['OutputDir'])
#Path to the directory to save the log file. Default is a directory called log in the directory where the output files are saved
log_loc     = os.path.expandvars(options['LogDir'])
#Path to the directory to save the script files. Default is a directory called log in the directory where the output files are saved
script_loc     = os.path.expandvars(options['ScriptDir'])
#Name of the default slurm script file
default_slurm  = os.path.expandvars(options['SLURMScript'])
#Name of the default xml script file
default_xml    = os.path.expandvars(options['XMLScript'])
#Path to the source file to source in slurm script 
source_file     = os.path.expandvars(options['SourceFile'])
#The prefix for the slurm job name
job_name_prefix  = options['JobNamePrefix']
#The postfix string to add to the output file
run_type        = options['NewFilePostfix']
#Rewrite all files. Default is false.
rewrite        = options.getboolean('Rewrite')
#Submit the scripts to the queue. Default is false.
submit         = options.getboolean('Submit')
#If string is non empty, it is passed to re.sub to change the format of filenames to allow more than one file per call of CCMAnalysis
multi_file_regex = options['MultiFileRegex']
#Default file format in regular expression format
default_file_regex = options['DefaultFileRegex']
#The maximum number of jobs that can be submitted. Default is 1000 (HPC max)
max_num_jobs     = options.getint('MaxNumJobs')
#The parameters to pass to ls to select the input files.
input_data_regex = options['InputDataRegex']
#The string that ends the input files (before the .root), used to get a unique name.
input_end_string = options['InputEndString']
#Flag to say if the files are from the MC or not
is_mc = options.getboolean("IsMC")

user    = os.getenv("USER")

# -------- main code ------------------------
#for key in options:
#  print(key)
print({section: dict(config[section]) for section in config.sections()})

input_loc   = os.path.abspath(input_loc)
output_loc  = os.path.abspath(output_loc)
log_loc  = os.path.abspath(log_loc)
script_loc  = os.path.abspath(script_loc)
default_slurm  = os.path.abspath(default_slurm)
default_xml  = os.path.abspath(default_xml)

#  will need to change conditions on what to look for
raw_list  = glob.glob(os.path.join(input_loc, input_data_regex))

#if user supplied a multi_file_regex string use it in re.sub to replace the filenames in raw_list and return a list of unique results
if multi_file_regex:
  pattern = re.compile(r""+default_file_regex,re.VERBOSE)
  replaced_filenames = [pattern.sub(multi_file_regex,file) for file in raw_list]             
  raw_list = list(set(replaced_filenames))

slurm_scripts = []
current_xml_list = []

xml_script = script_loc+"/"+os.path.basename(default_xml)

if not os.path.exists(xml_script):
  fin = open(default_xml,"rt")
  fout = open(xml_script,"wt")
  for line in fin:
    copyLine = line
    fout.write(copyLine)

  fin.close()
  fout.close()

run = ""

script_num = 0

count = 0

for raw_file in raw_list:
  baseName = os.path.basename(raw_file.replace(".root", "_%s.root"%(run_type)))
  if multi_file_regex:
    #remove regex characters from filenames
    baseName = re.sub('[\?\*]','x',baseName)

  out_name = output_loc+"/"+baseName

  if os.path.exists(out_name) and not rewrite:
    continue

  run = ""
  if not is_mc:
    start = baseName.find("run")
    end = baseName.find(input_end_string)
    run = baseName[start:end]
    job_name = job_name_prefix+'_'+run

  if run == "":
    run = raw_file.replace(".root","")
    run = run.replace(input_loc+"/","")

  print("run ",run, " out ",out_name)

  log_script = log_loc+"/"+job_name+".log"

  current_xml_list.append("CCMAnalysis -i "+raw_file+" -o "+out_name+" -c "+xml_script+" -l "+log_script)

  print("len(current_xml_list) ",len(current_xml_list))

  if len(current_xml_list) == 34:
    script_num = script_num + 1
    slurm_script = script_loc+"/"+run_type+"_slurm_"+run+".script"
    app_script = script_loc+"/"+run_type+"_mpiapp_"+run
    slurm_scripts.append(slurm_script)

    if not os.path.exists(slurm_script) or rewrite:

      print("Create and Write SLURM script file")
      fin = open(default_slurm,"rt")
      fout = open(slurm_script,"wt")
      for line in fin:
        copyLine = line
        copyLine = copyLine.replace('SOURCEFILE',source_file)
        copyLine = copyLine.replace('JOBNAME',job_name)
        copyLine = copyLine.replace('LOGFILELOCATION',log_loc)
        copyLine = copyLine.replace('MPIAPP',app_script)
        fout.write(copyLine)
      fin.close()
      fout.close()

      print("Create and Write MPI APP File")
      fout = open(app_script,"wt")
      for line in current_xml_list:
        fout.write("--map-by node -n 1 "+line+"\n")
      fout.close()

    current_xml_list = []

  if len(slurm_scripts) == max_num_jobs:
    current_xml_list = []
    break

if len(current_xml_list) != 0:
    slurm_script = script_loc+"/"+run_type+"_slurm_"+run+".script"
    app_script = script_loc+"/"+run_type+"_mpiapp_"+run
    slurm_scripts.append(slurm_script)

    if not os.path.exists(slurm_script) or rewrite:

      print("Create and Write SLURM script file")
      fin = open(default_slurm,"rt")
      fout = open(slurm_script,"wt")
      for line in fin:
        copyLine = line
        copyLine = copyLine.replace('SOURCEFILE',source_file)
        copyLine = copyLine.replace('JOBNAME',job_name)
        copyLine = copyLine.replace('LOGFILELOCATION',log_loc)
        copyLine = copyLine.replace('MPIAPP',app_script)
        fout.write(copyLine)
      fin.close()
      fout.close()

      print("Create and Write MPI APP File")
      fout = open(app_script,"wt")
      for line in current_xml_list:
        fout.write("--map-by node -n 1 "+line+"\n")
      fout.close()

    current_xml_list = []


out = subprocess.run(["squeue","-l", "-u", user], stdout=PIPE, universal_newlines=True)
num_lines = len(re.findall("\w*"+user+"\w*",out.stdout))

if submit:
  n = min(max_num_jobs,len(slurm_scripts))
  if n+num_lines > max_num_jobs:
    n = max_num_jobs-num_lines
  for slurm in slurm_scripts[:n]:
    subprocess.call(["sbatch",slurm])

print("MaxNumJobs ",max_num_jobs," num_lines ",num_lines," numScripts ",len(slurm_scripts))

