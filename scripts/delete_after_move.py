import os
import glob
import time

delete_location = "/data/to_delete/"

while(True):
    fnames = sorted(glob.glob(delete_location+"*.bin"))
    if len(fnames) > 10:
        print("Removing files:\n" + "\n".join(["\t" + s for s in fnames[:-10]]))
        for i in range(len(fnames) - 10):
            print("os.remove(" + fnames[i] + ")")
            os.remove(fnames[i])
            #os.rename(fnames[i],fnames[i]+"_DELETED")

    time.sleep(10)
