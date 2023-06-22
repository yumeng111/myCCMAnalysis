
import os
import glob
import numpy as np
import argparse

import convert_to_i3
    
parser = argparse.ArgumentParser(
        prog = "Convert_and_move",
        description = "Converts binary to i3.zst and moves binary to landing pad",
)

parser.add_argument("--machine",
        type=str,
        required=True,
        help="Machine: \"m\" (\"w\") for m(w)illstester")

parser.add_argument("--nchunks",
        type=int,
        required=True,
        help="number of chunks for input filelist")

parser.add_argument("--chunk",
        type=int,
        required=True,
        help="chunk to convert and move, 0 < chunk < nchunks-1")

args = parser.parse_args()

if args.machine=="m": 
    filelist_file = "millstester_filelist.txt"
    source = "/data/2022_run/binary_files/"
    comp_destination = "/data/2023/compressed_unmerged_files/millstester"
elif args.machine=="w": 
    filelist_file = "willstester_filelist.txt"
    source = "/data/2022_run_m2/binary_files/"
    comp_destination = "/data/2023/compressed_unmerged_files/willstester"
elif args.machine=="t": 
    filelist_file = "test_filelist.txt"
    source = "/data/test_source/"
    comp_destination = "/data/test_destination/"
else:
    print("Give a valid machine!")
    exit(0)

with open(filelist_file) as file:
    filelist = [line.strip() for line in file]

runs_to_compress = set(np.arange(12010, 12500))

delete_location = "/data/to_delete/"

def chunk(items, chunks, chunk_number):
    if chunks > 0:
        a = int(np.floor(float(len(items)) / float(chunks)))
        b = int(np.ceil(float(len(items)) / float(chunks)))
        x = len(items) - a*chunks
        if chunk_number < x:
            n = b
            n0 = n*chunk_number
        else:
            n = a
            n0 = b*x + a*(chunk_number-x)
        n1 = min(n0+n,len(items))
        if n0 >= len(items):
            return []
        items = items[n0:n1]
    return items

def get_files_to_convert(path):
    fnames = np.array(sorted(glob.glob(path + "/ccmdata_*.bin")))
    #print(fnames)
    run_numbers = np.array([int(fname.split("/")[-1].split("_")[1][3:]) for fname in fnames])
    run_mask = np.array([bool(run in runs_to_compress) for run in run_numbers]).astype(bool)
    fnames = fnames[run_mask]
    return sorted(fnames)

def get_ouput_file_name(fname, path):
    return path + "/" + fname.split("/")[-1][:-3] + "i3.zst"

def get_destination_file_name(fname, path):
    return path + "/" + fname.split("/")[-1]

def convert_file(fname, destination):
    output = get_ouput_file_name(fname, destination)
    convert_to_i3.convert([[fname]], output)

def move_file(fname, destination):
    output = get_destination_file_name(fname, destination)
    os.rename(fname, output)

def check_output(inputfile, outputfile):
    if not os.path.exists(outputfile): return False
    input_size = os.path.getsize(inputfile)
    output_size = os.path.getsize(outputfile)
    return output_size/input_size > 0.1


files_to_convert = chunk(filelist, args.nchunks, args.chunk)


for fname in files_to_convert:
    if not os.path.exists(fname):
        continue
    if os.path.getsize(fname) <= 8191:
        move_file(fname, delete_location)
        continue

    convert_file(fname,comp_destination)
    if check_output(fname, get_ouput_file_name(fname, comp_destination)):
        move_file(fname, delete_location)
    




