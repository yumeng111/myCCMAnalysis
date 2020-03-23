#! /usr/bin/python3

import sys, os
import subprocess
import time
import glob
import datetime

# -------- settings ------------------------
inputDir      = os.getenv("CCMRAW")+"/"
outputDir     = os.getenv("CCMPROCESSED")+"/"
defaultSLRUM  = os.getenv("LUSTRETMP")+"/findEvents_slrum.script"
defaultXML    = os.getenv("CCMINSTALL")+"/configFiles/local/findEventsConfig.xml"
hvOffList     = os.getenv("CCMINSTALL")+"/mappings/hv_off_2019.csv"
calibrationFile = os.getenv("CCMINSTALL")+"/mappings/root_out_2019ledData_run179_legFix_integral_all_round4_.root"
sourceFile     = os.getenv("HOME")+"/.bash_profile"
logFileLoc     = os.getenv("CCMINSTALL")+"/log"
runType       = "events"

# -------- main code ------------------------
inputDir   = os.path.abspath(inputDir)
outputDir  = os.path.abspath(outputDir)

#  will need to change conditions on what to look for
rawList  = glob.glob(os.path.join(inputDir, "*run000167*.root"))

slrumScripts = []

for rawFile in rawList:
  baseName = os.path.basename(rawFile.replace(".root", "_%s.root"%(runType)))
  outName = outputDir+"/"+baseName
  start = baseName.find("run")
  end = baseName.find("_physics")
  run = baseName[start:end]
  print("baseName",baseName," run ",run, " out ",outName)

  slrumScript = outputDir+"/"+"findEvents_slrum_"+run+".script"
  xmlScript = outputDir+"/"+"findEventsConfig_"+run+".xml"

  fin = open(defaultSLRUM,"rt")
  fout = open(slrumScript,"wt")
  for line in fin:
    copyLine = line
    copyLine = copyLine.replace('CONFIGFILE',xmlScript)
    copyLine = copyLine.replace('RUN',run)
    copyLine = copyLine.replace('LOGFILELOCATION',logFileLoc)
    copyLine = copyLine.replace('SOURCEFILE',sourceFile)
    fout.write(copyLine)

  fin.close()
  fout.close()

  fin = open(defaultXML,"rt")
  fout = open(xmlScript,"wt")
  for line in fin:
    copyLine = line
    copyLine = copyLine.replace('INFILENAME',rawFile)
    copyLine = copyLine.replace('OUTFILENAME',outName)
    copyLine = copyLine.replace('HVOFFLIST',hvOffList)
    copyLine = copyLine.replace('CALIBRATIONFILE',calibrationFile)
    copyLine = copyLine.replace('TRIGGER',"BEAM")
    fout.write(copyLine)

  fin.close()
  fout.close()

  slrumScripts.append(slrumScript)

for slrum in slrumScripts:
  subprocess.call(["sbatch",slrum])


"""
fin = open("/lustre/scratch3/turquoise/rtthorn/findEvents_slrum.script", "rt")
fout = open("/lustre/scratch3/turquoise/rtthorn/findEvents_slrum_specific.script", "wt")

for line in fin:
	fout.write(line.replace('CONFIGFILE', '$CCMINSTALL/configFiles/local/findEventsConfig.xml'))
	
fin.close()
fout.close()
"""
