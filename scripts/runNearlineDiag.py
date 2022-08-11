#! /bin/python

# Written by Remington T. Thornton on 08/09/2019
# Based on convertBinary2ROOT.py En-Chuan Huang on 05/02/2019
# Function:
#	Run NearlineDaig every 10 min on the previous 10 min worth of data
#
# Note:
# 1. If executable requires env variables, please set it first.
# 2. The program will check files in the [outputDir] that are less than 10 min old

import sys, os
import subprocess, shlex
import time
import glob
import datetime

# -------- settings ------------------------
commandFormat = "CCMAnalysis -c {configFileName} -I {inputFileName} -o {outputFileName} -l {logFileName}"#command for CCMAnalysis(note -I takes list of files

configFileName    = "/home/daquser/CCMCode/CCM_analysis_sw_dev/configFiles/speEventFinderRateCalc.xml"
inputFileName = 'files_for_nearlineDiag.txt'
#location for output file and log, change to make legitimate
outputFileName  = "/dev/null/out.out"
logFileName = "/dev/null/log.log"
logFile = "logFile"#log file for process.Popen

#old paths from CCM2019
#executable    = "/home/daquser/CCMCode/build/NearlineDaig"
#outputDir     = "/home/daquser/CCM_DAQ_test/data/compressed"
#outputDirA     = "/mnt/CCM_Data_A/CCM_Data_A"
#outputDirA     = "/mnt/CCM_Data_A/CCM_Data_A"
#outputDirB     = "/mnt/CCM_Data_B/CCM_Data_B"
#new paths
outputDirB      = "/mnt/CCM_Data_B/2022_data/root_files/led_calibration/"
outputDirA = outputDirB
outputDir = outputDirB
runType       = "waveforms"
env_setup     = None  # if environment is needed, please uncomment below as well.
delayTime     = 10 * 60. # binary file will only be converted after last modified time is more than 10 mins ago
# ------------------------------------------

def printProcesses(myList):
  # each item inside list = [name, Popen, starttime, logFile]
  if len(myList)==0:
    print "No process running"
  for item in myList:
    print item[0]," Status: ", item[1].poll(), " running for %.1fs"%(time.time()-item[2])

def cleanProcesses(myList): # clean up finished processes and print exit code
  newList = []
  for item in myList:
    if item[1].poll()==None:
      newList.append(item)
    else:
      print "%s finished with exit code %d in %.0fs "%(item[0], item[1].poll(),time.time()-item[2])
      item[3].close()
  return newList

def countProcesses(myList): # count number of running processes
  if len(myList)==0:
    return 0
  return sum([item[1].poll()==None for item in myList])

# If environment

#os.system("source %s"%(env_setup))


while True: # continue searching for new files
    # avoid miscommunication, use absolute path
    #executable = os.path.abspath(executable)
    outputDirA  = os.path.abspath(outputDirA)
    outputDirB  = os.path.abspath(outputDirB)

    # List files only contain the base name with .bin/.root or runType
    timenow = time.time()
    #rootListA    = glob.glob(os.path.join(outputDirA, "*_%s.root"%(runType)))
    rootListB    = glob.glob(os.path.join(outputDirB, "*_%s.root"%(runType)))
    rootList = []

    """
    #removed because diskA is broken
    for fileName in rootListA:
    if os.path.exists(fileName):
      if timenow - os.path.getmtime(fileName) < delayTime and timenow-os.path.getmtime(fileName) > 60:
        rootList.append(fileName)
    """
    for fileName in rootListB:
        if os.path.exists(fileName):
          if timenow - os.path.getmtime(fileName) < delayTime and timenow-os.path.getmtime(fileName) > 60:
            rootList.append(fileName)

    if len(rootList)!=0:
        print "New Files to be processed: ", rootList
    else:
        print "No files to process for NearlineDiag sleeping for",delayTime,"s",time.asctime(time.localtime())
        time.sleep(delayTime)
        continue

    with open(inputFileName, 'w') as f:
        for item in  rootList:
            f.write("%s\n" % item)

    print "Starting ",time.asctime(time.localtime())

    command = commandFormat.format(configFileName=configFileName, inputFileName = inputFileName,outputFileName = outputFileName,logFileName = logFileName)
    command = shlex.split(command)
    print(command)
    #process = subprocess.Popen(command, shell=False, stdout=logFile,stderr=logFile)
    process = subprocess.Popen(command, shell=False)

    #  subprocess.call(["NearlineDiag","files_for_nearlineDiag.txt"])

    print "Sleeping for",(delayTime-300),"s",time.asctime(time.localtime())
    time.sleep(delayTime-300)

