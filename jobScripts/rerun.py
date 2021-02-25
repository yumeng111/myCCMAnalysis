#! /bin/python3

import sys, os
import subprocess
from subprocess import PIPE
import re
import time
import glob
from datetime import datetime
import argparse
import sys

maxLines = 300

def get_current_num_jobs():
  out = subprocess.run(["squeue","-l", "-u", "rtthorn"], stdout=PIPE, universal_newlines=True)
  return len(re.findall("\w*rtthorn\w*",out.stdout))

with open('rerunDM.txt') as f:
  lines = [line.rstrip() for line in f]

lenLines = len(lines)
for index,line in enumerate(lines):
  while get_current_num_jobs() >= maxLines:
    print(datetime.now(),"\tSleeping")
    time.sleep(60)

  print(datetime.now(),"\tsubmitting index ",index," out of ",lenLines," ",line)
  subprocess.call(["sbatch",line])

