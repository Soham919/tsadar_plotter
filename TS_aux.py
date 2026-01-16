import pandas as pd
import xarray as xr
from dataclasses import dataclass
from datetime import date

shot_data = pd.DataFrame({
    "shot" : [92522,92524,92525,92527,92528,92530,92531,92532,92533,92534,92535,92536,92537,92538],
    "pointing" : [4,2,4,6,8,5.2,4,4,5,3,6,5,4,4],
    "timing" : [4,2,4,6,8,6,4,4,4,4,6,6,6,4]
})

shot_data.set_index("shot", inplace=True)

def load_file(data, shot_no, file_tag, file_path):
    """
    Load multiple CSV files for a given shot number into a nested dictionary.
    """
    try:
        s = file_tag.split('.')
        if s[-1]=="csv":
            data[shot_no][file_tag] = pd.read_csv(file_path)
            print(f"Loaded {file_tag}")

        if s[-1]=="nc":
            data[shot_no][file_tag] = xr.open_dataset(file_path)
            print(data[shot_no][file_tag])
            
    except FileNotFoundError:
        print(f"⚠ File not found: {file_path}")