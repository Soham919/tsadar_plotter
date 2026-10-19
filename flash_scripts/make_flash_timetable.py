import os
from pathlib import Path
from flash_helpers import make_flash_time_table

# --------- Mac ----------- #
#flash_folder = Path("/Users/soham/Documents/Flash/2DCylindrical/Si3N4_cyl_5TW_1300psiGJS")

# ---- Windows ---- #
flash_folder = Path(r"C:\Simulation_data\FLASH\2D_Cylindrical\FV4.7.1_carbuncle_debugs\19A_hll1_2.5TW")  # directory containing the FLASH plot files


make_flash_time_table(
    fp=flash_folder,
    file_pattern="*hdf5_plt_cnt_*",
    output_txt="flash_plotfile_times.txt"
)