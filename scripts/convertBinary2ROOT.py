#! /bin/python

# Written by En-Chuan Huang on 06/03/2019
# Function:
#    convert *.bin in [inputDir] to *_[runType].root in [outputDir]
#    with with [executable] and the [compressLevel] at maximum [MAX_NTHREADS].
#    Output and error message is saved to [outputDir]/*_[runType].log 
#
# Note:
# 1. If executable requires env variables, please set it first.
# 2. The program will check new files in the [inputDir] to convert
# 3. Currently, the program does NOT remove/move the .bin file
# 4. One can change the order in [commandFormat] to suit the executable

import sys, os
import subprocess, shlex
import time
import glob
import datetime

# -------- settings ------------------------
MAX_NTHREADS  = 4
commandFormat = "CCMAnalysis -c {configFileName} -r {inputFileName} -o {outputFileName} -l {logFileName}"

configFileName    = "/home/daquser/CCMCode/daq_machine_software/sources/CCM_analysis_200/configFiles/findPulsesConfig.xml"
logDir    = "/home/daquser/CCM_DAQ_Nboard/data/root_files/find_pulses_log_files"

#inputDir      = "/home/daquser/CCM_DAQ_Nboard/data/"


inputDir      = "/home/daquser/CCM_DAQ_Nboard/data/"
#inputDir =     "/mnt/CCM_Data_B/LED_2021/"
#outputDir     = "/home/daquser/CCMCode/ROOT_files/"
#outputDir     = "/home/daquser/CCM_DAQ_test/data/compressed/"
#outputDirA    = "/mnt/CCM_Data_A/CCM_Data_A/"
#outputDirB    = "/mnt/CCM_Data_B/CCM_Data_B/"
#outputDirB      = "/mnt/CCM_Data_B/LED_2021/root_files"#2021 output directory
#outputDirB      ="/mnt/CCM_Data_B/2022_data/root_files/led_calibration"#2022 run led callibration
outputDirB      = "/home/daquser/CCM_DAQ_Nboard/data/root_files"
#outputDir     = "/home/daquser/CCM_DAQ_1board/data/"
outputDir     = outputDirB
runType       = "waveforms"
env_setup     = None  # if environment is needed, please uncomment below as well.
delayTime     = 1 * 60. # binary file will only be converted after last modified time is more than 10 mins ago

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
      timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
      print "%s: %s finished with exit code %d in %.0fs "%(timestamp,item[0], item[1].poll(),time.time()-item[2])
      item[3].close()
  return newList

def countProcesses(myList): # count number of running processes
  if len(myList)==0:
    return 0
  return sum([item[1].poll()==None for item in myList])

# If environment

#os.system("source %s"%(env_setup))

processingList = [] # each item [name, Popen, starttime]
workingList = [] # contain working file names

truth = True
while (truth): # continue searching for new files
  # avoid miscommunication, use absolute path
  inputDir   = os.path.abspath(inputDir)
  outputDir  = os.path.abspath(outputDir)
  configFileName = os.path.abspath(configFileName)
  logDir = os.path.abspath(logDir)

  # List files only contain the base name with .bin/.root or runType
  timenow = time.time()
  binaryList  = glob.glob(os.path.join(inputDir, "*.bin"))
  rootListA   = glob.glob(os.path.join(outputDir, "*_%s.root"%(runType)))
  rootListB   = glob.glob(os.path.join(outputDirB, "*_%s.root"%(runType)))
  logListA    = glob.glob(os.path.join(outputDir, "*_%s.log"%(runType)))
  logListB    = glob.glob(os.path.join(outputDirB, "*_%s.log"%(runType)))
  binaryList  = [os.path.basename(fileName.replace(".bin", "")) for fileName in binaryList \
 if timenow-os.path.getmtime(fileName)>delayTime]

  rootListCheck = []
  for fileName in rootListA:
    if os.path.exists(fileName):
      if timenow - os.path.getmtime(fileName) > delayTime:
        f = fileName.replace("_%s.root"%(runType),"")
        rootListCheck.append(os.path.basename(f))

  for fileName in rootListB:
    if os.path.exists(fileName):
      if timenow - os.path.getmtime(fileName) > delayTime:
        f = fileName.replace("_%s.root"%(runType),"")
        rootListCheck.append(os.path.basename(f))

  rootListA    = [os.path.basename(fileName.replace("_%s.root"%(runType), "")) for fileName in rootListA]
  logListA     = [os.path.basename(fileName.replace("_%s.log"%(runType), "")) for fileName in rootListA]
  rootListB    = [os.path.basename(fileName.replace("_%s.root"%(runType), "")) for fileName in rootListB]
  logListB     = [os.path.basename(fileName.replace("_%s.log"%(runType), "")) for fileName in rootListB]
  rootList    = rootListA + rootListB + logListA + logListB # logList to make sure no one is double counted
  convertList = list(set(binaryList) - set(rootList))
  convertList.sort()
  convertListCheck = [value for value in binaryList if value in rootListCheck]
  convertListCheck.sort()

  processingList = cleanProcesses(processingList)
  # before running double check previous processesm
  if countProcesses(processingList)!=0:
    print "Running Threads:"
    printProcesses(processingList)
    print ""

  if len(convertList)!=0:
    print "New Files to be processed: ", convertList
  else:
    print "No new file to process at",time.asctime(time.localtime())
    time.sleep(delayTime)
  print "Searching for Files to remove..."
  for itemName in convertListCheck:
    print "Removed "+itemName+".bin"
    os.remove(os.path.join(inputDir,itemName+".bin"))

  # loop through convertList
  for itemName in convertList:
    # submit a new process for each file
    # first find out how much disk space is left
    #10/5 - WT commented out 
  #  out = subprocess.Popen(['parseDFOutput.sh', '/mnt/CCM_Data_A', '11'],
  #    stdout=subprocess.PIPE,
  #    stderr=subprocess.STDOUT)

  #  ccm_data_a,junk = out.communicate()
  #  ccm_data_a = str(ccm_data_a)
  #  ccm_data_a = ccm_data_a.replace("b'","")
  #  ccm_data_a = ccm_data_a.replace("\\n'","")
  #  ccm_data_a = int(ccm_data_a[:-2])
  #  print "Disk Space on /mnt/CCM_Data_A/ = ",ccm_data_a

#    out = subprocess.Popen(['parseDFOutput.sh', '/mnt/CCM_Data_B', '11'],
#      stdout=subprocess.PIPE,
#      stderr=subprocess.STDOUT)
#
#    ccm_data_b,junk = out.communicate()
#    ccm_data_b = str(ccm_data_b)
#    ccm_data_b = ccm_data_b.replace("b'","")
#    ccm_data_b = ccm_data_b.replace("\\n'","")
#    ccm_data_b = int(ccm_data_b[:-2])
#    print "Disk Space on /mnt/CCM_Data_B/ = ",ccm_data_b

    # choose disk based on available disk space
#    if ccm_data_a < 98:
#      outputDir = outputDirA
#    elif ccm_data_b < 98:
#      outputDir = outputDirB
#    outputDir = inputDir

    print "going to save to ",outputDir

    inputFileName = os.path.join(inputDir, itemName+".bin")
    outputFileName = os.path.join(outputDir, itemName + "_"+runType+ ".root")
    outputFileName = os.path.join(outputDir, itemName + "_"+runType+ ".root")
    logFileName = os.path.join(logDir, itemName + "_"+runType+ ".log")
    logFile = open(outputFileName.replace(".root", ".log"), "w")
    command = commandFormat.format(configFileName=configFileName, inputFileName = inputFileName, outputFileName = outputFileName,logFileName=logFileName)
    command = shlex.split(command)
    print(command)
    process = subprocess.Popen(command, shell=False, stdout=logFile,stderr=logFile)
    t = time.time()
    processingList.append([itemName, process, t, logFile])


    # no more than MAX_NTHREADS processes running 
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print "%s: Start Processing %s.bin while %d threads running"%(timestamp,itemName, countProcesses(processingList))
    while True:
      if countProcesses(processingList)<MAX_NTHREADS:
        processingList = cleanProcesses(processingList)
        break
      else:
        time.sleep(delayTime)

  #truth = False

  time.sleep(delayTime)


# finishing
while True:
  if countProcesses(processingList)>0:
    processingList = cleanProcesses(processingList)
    time.sleep(delayTime)
    continue
  else:
    break
