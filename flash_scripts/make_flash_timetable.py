import os
from pathlib import Path
from flash_helpers import make_flash_time_table

# --------- Mac ----------- #
flash_folder = Path("/Users/soham/Documents/Flash/test_runs/1mmSpot_8mm_offset")

# ---- Windows ---- #
#flash_folder = r"\\profiles\Users$\sban\Documents\FLASH\2D\10um_spot_2"
#flash_folder = Path(r"C:\Simulation_data\FLASH\2D\1mm_spot_Si_2")

make_flash_time_table(
    fp=flash_folder,
    file_pattern="*hdf5_plt_cnt_*",
    output_txt="flash_plotfile_times.txt"
)