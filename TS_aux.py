import pandas as pd
import xarray as xr
from dataclasses import dataclass
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re


shot_data = pd.DataFrame({
    "shot" : [92522,92524,92525,92527,92528,92530,92531,92532,92533,92534,92535,92536,92537,92538,117830,117831,117832,117837,117838,117839],
    "pointing" : [4,2,4,6,8,5.2,4,4,5,3,6,5,4,4,8,8,5.5,5.5,5.5,5.5],
    "timing" : [4,2,4,6,8,6,4,4,4,4,6,6,6,4,9,9,7.6,7.3,6,5]
})

shot_data.set_index("shot", inplace=True)

def load_file(data, file_tag, file_path):
    """
    Load multiple CSV files for a given shot number into a nested dictionary.
    """
    try:
        s = file_tag.split('.')
        if s[-1]=="csv":
            data[file_tag] = pd.read_csv(file_path)
            print(f"Loaded {file_tag}")

        if s[-1]=="nc":
            data[file_tag] = xr.open_dataset(file_path)
            print(data[file_tag])
            
    except FileNotFoundError:
        print(f"⚠ File not found: {file_path}")

def select_shots(folder_path):
    """
    Select specific shots from the data dictionary.
    """
    ls = [f for f in os.listdir(folder_path) if os.path.isdir(Path(folder_path)/f)]
    data = {}

    print("Stuff in folder :\n")
    for i in range(len(ls)):
        print(f"[{i:g}] {ls[i]}\n")
        
    while True:
        n = int(input("Which shot would you like to see fits for?\n"))
        fp2 = os.listdir(Path(folder_path)/ls[n])
        if str(ls[n]) not in data:
            data[str(ls[n])] = {}
        lscsv = [f for f in fp2]
        print(f"Files for shot {ls[n]} :\n")
        for i in range(len(lscsv)):
            print(f"[{i}] {lscsv[i]}")
        
        fn = input("Which files do you want to import?\n")
        files  = re.split(r"[ ,]+", fn)
        files = [int(x) for x in files]
        for i in files:
            fp = Path(folder_path)/ls[n]/lscsv[i]
            if fp.is_dir():
                if str(lscsv[i]) not in data[str(ls[n])]:
                    data[str(ls[n])][str(lscsv[i])] = {}
                files2 = os.listdir(fp)
                for f in files2:
                    load_file(data[str(ls[n])][str(lscsv[i])],f,Path(fp)/f)
            else:
                load_file(data[str(ls[n])],lscsv[i],fp)
                
        flag = input("Load more shots?\n")
        if flag in ("n","no","nah","","nope"):
            break

    return data

