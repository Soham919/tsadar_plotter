import os
from flash_helpers import make_flash_time_table

# ---- Windows ---- #
flash_folder = r"\\profiles\Users$\sban\Documents\FLASH\1D\1mm_spot\1D_1um_resolve"

make_flash_time_table(
    fp=flash_folder,
    file_pattern="*hdf5_plt_cnt_*",
    output_txt="flash_plotfile_times.txt"
)