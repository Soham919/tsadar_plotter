import subprocess
from pathlib import Path

# -------- user settings --------
username = "sban"
host = "cl7head"

remote_dir = "/home/sban/Vacuum_6_EnergyBins"  # remote directory containing the files to copy
#remote_dir = "/g1/hdd/sban/FLASH_runs/1D/1mm_spot/1D2Dcompare_Si_3TW"
file_pattern = "*"          # examples: "*.h5", "*.png", "*.txt"
#file_pattern = "*LaserEnergyProfile.dat"        
# ---- Mac ---- #
# local_dir = "/Users/soham/Documents/Flash/test_runs/1mmSpot_8mm_offset"

# ---- Windows ---- #
local_dir = Path(r"C:\Simulation_data\PROPACEOS\Vacuum_6_EnergyBins")  # local directory to copy files to
# -------------------------------

def scp_copy_files(username, host, remote_dir, file_pattern, local_dir):
    """
    Copy matching files from a remote host to a local directory using scp.

    Parameters
    ----------
    username : str
        Remote username.
    host : str
        Remote hostname.
    remote_dir : str
        Remote directory containing the files.
    file_pattern : str
        File wildcard pattern, e.g. '*.h5' or '*.png'.
    local_dir : str
        Local destination directory.
    """
    local_path = Path(local_dir).expanduser()
    local_path.mkdir(parents=True, exist_ok=True)

    remote_path = f'{username}@{host}:{remote_dir.rstrip("/")}/{file_pattern}'

    cmd = ["scp", remote_path, str(local_path)]

    print("Running command:")
    print(" ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
        print(f"\nTransfer complete to: {local_path}")
    except subprocess.CalledProcessError as e:
        print("\nscp failed.")
        print(e)

if __name__ == "__main__":
    scp_copy_files(username, host, remote_dir, file_pattern, local_dir)