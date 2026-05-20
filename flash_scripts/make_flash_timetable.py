import os
from flash_helpers import make_flash_time_table

# ---- Windows ---- #
flash_folder = r"\\profiles\Users$\sban\Documents\FLASH\1D\1mm_spot\1D_Si_4.16TW_long"

make_flash_time_table(
    fp=flash_folder,
    file_pattern="4.16TW_hdf5_plt_cnt_*",
    output_txt="flash_plotfile_times.txt"
)