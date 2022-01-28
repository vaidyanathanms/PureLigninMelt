# Aux python codes for generating GROMACS inputs for analysis
# Ver: May-11-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
import datetime
#------------------------------------------------------------------

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        raise RuntimeError('ERROR: ' + srcfyl + ' not found')

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

# Set working directory
def set_working_dir(rundir,inp_type,solv_type = 'None'):
    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, "does not exist")
    
    return workdir1
#------------------------------------------------------------------

#Check pdb/psf/top files for the melt (or polymer)
def find_latest_file(dum_inpdir,ext):
    # check structure files (*.pdb/.gro)
    if glob.glob(dum_inpdir+'/*' + ext) == []:
        print("No files with " + ext + " found")
        return 'None'
    elif len(glob.glob(dum_inpdir+'/*' + ext)) == 1:
        out_fname = glob.glob(dum_inpdir+'/*' + ext)[0]
    elif len(glob.glob(dum_inpdir+'/*' + ext)) > 1:
        print("More than one file with " + ext + " found")
        fnames = glob.glob(dum_inpdir+'/*' + ext)
        out_fname = max(fnames, key=os.path.getmtime)

    return out_fname
#------------------------------------------------------------------

# Extract filename from full path
def extract_file_name(fname):
    outname = fname.split("/")
    return outname[len(outname)-1]
#------------------------------------------------------------------

# Read file in reverse order
#https://thispointer.com/python-read-a-file-in-reverse-order-line-by-line/
def read_reverse_order(file_name):
    # Open file for reading in binary mode
    with open(file_name, 'rb') as read_obj:
        # Move the cursor to the end of the file
        read_obj.seek(0, os.SEEK_END)
        # Get the current position of pointer i.e eof
        pointer_location = read_obj.tell()
        # Create a buffer to keep the last read line
        buffer = bytearray()
        # Loop till pointer reaches the top of the file
        while pointer_location >= 0:
            # Move the file pointer to the location pointed by pointer_location
            read_obj.seek(pointer_location)
            # Shift pointer location by -1
            pointer_location = pointer_location -1
            # read that byte / character
            new_byte = read_obj.read(1)
            # If the read byte is new line character then it means one line is read
            if new_byte == b'\n':
                # Fetch the line from buffer and yield it
                yield buffer.decode()[::-1]
                # Reinitialize the byte array to save next line
                buffer = bytearray()
            else:
                # If last read character is not eol then add it in buffer
                buffer.extend(new_byte)
        # As file is read completely, if there is still data in buffer, then its the first line.
        if len(buffer) > 0:
            # Yield the first line too
            yield buffer.decode()[::-1]

#--------------------------------------------------------------------

# if __name__ 
if __name__ == 'main':
   main()
#--------------------------------------------------------------------
