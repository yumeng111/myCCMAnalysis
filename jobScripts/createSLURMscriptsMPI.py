#! /usr/bin/python3

import sys, os
import subprocess
from subprocess import PIPE
import re
import time
import glob
import datetime
import argparse
import sys

# -------- default settings ------------------------
inputDir      = os.getenv("CCMRAW")+"/"
outputDir     = os.getenv("CCMPROCESSED")+"/na22_20191221/"
defaultSLURM  = "findEvents_slurm.script"
defaultXML    = "findEventsConfig.xml"
hvOffList     = os.getenv("CCMINSTALL")+"/mappings/hv_off_2019.csv"
calibrationFile = os.getenv("CCMINSTALL")+"/mappings/root_out_2019ledData_run179_legFix_integral_all_round4_.root"
sourceFile     = os.getenv("HOME")+"/.bash_profile"
logFileLoc     = outputDir+"/log"
jobNamePrefix  = "findEvents_"
runType        = "events"
rewrite        = False
submit         = False
currentUser    = os.getenv("USER")
maxNumJobs     = 1000 # this is a hardcoded limit from LANL HPC
inputDataRegex = "*run000193*.root"
inputEndString = "_physics"

def getOptions(args=sys.argv[1:]):
  parser = argparse.ArgumentParser(description="Creates slurm and xml config file for the find events module. Submits the slurm script if --submit option is given.")
  parser.add_argument("--input-loc", default=inputDir,
      help="Directory location of the input files")
  parser.add_argument("--output-loc", default=outputDir,
      help="Directory location of the output files")
  parser.add_argument("--default-slurm", default=defaultSLURM,
      help="Name of the default slurm script file")
  parser.add_argument("--default-xml", default=defaultXML,
      help="Name of the default xml script file")
  parser.add_argument("--hv-off-list", default=hvOffList,
      help="Path to the HV off csv file to use") 
  parser.add_argument("--calibration-file", default=calibrationFile,
      help="Path to the calibration ROOT file to use") 
  parser.add_argument("--source-file", default=sourceFile,
      help="Path to the source file to source in slurm script") 
  parser.add_argument("--log-file", default=logFileLoc,
      help="Path to the directory to save the log file. Default is a directory called log in the directory where the output files are saved")
  parser.add_argument("--job-name-prefix", default=jobNamePrefix,
      help="The prefix for the slurm job name")
  parser.add_argument("--run-type", default=runType,
      help="The postfix string to add to the output file")
  parser.add_argument("--rewrite", default=rewrite, dest='rewrite', action='store_true',
      help="Rewrite all files. Default is false.")
  parser.add_argument("--submit", default=submit, dest='submit', action='store_true',
      help="Submit the scripts to the queue. Default is false.")
  parser.add_argument("--user", default=currentUser,
      help="Name of user submitting the jobs")
  parser.add_argument("--max-jobs", default=maxNumJobs, type=int,
      help="The maximum number of jobs that can be submitted. Default is 1000 (HPC max)")
  parser.add_argument("--input-data-regex", default=inputDataRegex,
      help="The parameters to pass to ls to select the input files.")
  parser.add_argument("--input-end", default=inputEndString,
      help="The string that ends the input files (before the .root), used to get a unique name.")
  parser.add_argument("--trigger-type", default="BEAM",
      help="Which triggers to look at. See xml configuration file for possible options. Default = BEAM")
  options = parser.parse_args(args)
  return options

# -------- main code ------------------------
options = getOptions(sys.argv[1:])
print(options)

inputDir   = os.path.abspath(options.input_loc)
outputDir  = os.path.abspath(options.output_loc)

#  will need to change conditions on what to look for
rawList  = glob.glob(os.path.join(inputDir, options.input_data_regex))

slurmScripts = []
currentXMLList = []

xmlScript = outputDir+"/scripts/"+"findEventsConfig.xml"

if not os.path.exists(xmlScript):
  fin = open(options.default_xml,"rt")
  fout = open(xmlScript,"wt")
  for line in fin:
    copyLine = line
    copyLine = copyLine.replace('HVOFFLIST',options.hv_off_list)
    copyLine = copyLine.replace('CALIBRATIONFILE',options.calibration_file)
    copyLine = copyLine.replace('TRIGGER',options.trigger_type)
    copyLine = copyLine.replace('THRESHOLD',"0.4")
    fout.write(copyLine)

  fin.close()
  fout.close()


count = 0
for rawFile in rawList:
  baseName = os.path.basename(rawFile.replace(".root", "_%s.root"%(options.run_type)))
  outName = outputDir+"/"+baseName

  if os.path.exists(outName) and not options.rewrite:
    continue

  start = baseName.find("run")
  end = baseName.find(inputEndString)
  run = baseName[start:end]
  jobName = options.job_name_prefix+run

  print("run ",run, " out ",outName)

  logScript = outputDir+"/log/"+"findEvents_"+run+".log"

  currentXMLList.append("CCMAnalysis -i "+rawFile+" -o "+outName+" -c "+xmlScript+" -l "+logScript)

  print("len(currentXMLList) ",len(currentXMLList))

  if len(currentXMLList) == 12:
    slurmScript = outputDir+"/scripts/"+"findEvents_slurm_"+run+".script"
    appScript = outputDir+"/scripts/"+"findEvents_mpiapp_"+run
    slurmScripts.append(slurmScript)

    if not os.path.exists(slurmScript) or options.rewrite:

      if options.log_file == logFileLoc:
        options.log_file = options.output_loc+"/log"

      print("Create and Write SLURM script file")
      fin = open(options.default_slurm,"rt")
      fout = open(slurmScript,"wt")
      for line in fin:
        copyLine = line
        copyLIne = copyLine.replace('SOURCEFILE',options.source_file)
        copyLine = copyLine.replace('CONFIGFILE',xmlScript)
        copyLine = copyLine.replace('JOBNAME',jobName)
        copyLine = copyLine.replace('LOGFILELOCATION',options.log_file)
        copyLine = copyLine.replace('MPIAPP',appScript)
        fout.write(copyLine)
      fin.close()
      fout.close()

      print("Create and Write MPI APP File")
      fout = open(appScript,"wt")
      for line in currentXMLList:
        fout.write("--map-by node -n 1 "+line+"\n")
      fout.close()

    currentXMLList = []

  if len(slurmScripts) == options.max_jobs:
    break


out = subprocess.run(["squeue","-l", "-u", options.user], stdout=PIPE, universal_newlines=True)
numLines = len(re.findall("\w*"+options.user+"\w*",out.stdout))

if options.submit:
  n = min(options.max_jobs,len(slurmScripts))
  if n+numLines > maxNumJobs:
    n = maxNumJobs-numLines
  for slurm in slurmScripts[:n]:
    subprocess.call(["sbatch",slurm])

print("MaxNumJobs ",options.max_jobs," numLines ",numLines," numScripts ",len(slurmScripts))

