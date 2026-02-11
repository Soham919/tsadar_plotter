import h5py
import os 
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def print_h5(name, obj):
    print(name)

def h5_to_dict(node):
    """Recursively converts an HDF5 group/file into a nested dictionary."""
    data_dict = {}
    
    # Iterate through all items in the current group
    for key, item in node.items():
        if isinstance(item, h5py.Group):
            # If it's a group, recurse deeper
            #data_dict[key] = h5_to_dict(item)
            continue
        elif isinstance(item, h5py.Dataset):
            print(f"Loading dataset: {key} ({item.attrs['unit']}) with shape {item.shape} and dtype {item.dtype}")
            # If it's a dataset, read it into memory as a NumPy array [:]
            data_dict[key] = item[:]
            
    # Optional: Include attributes (metadata) if your file uses them
    if node.attrs.keys():
        data_dict['_attrs'] = dict(node.attrs)
        
    return data_dict