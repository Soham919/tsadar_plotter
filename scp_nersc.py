import subprocess
from pathlib import Path

# -------- user settings --------
username = "soham_b"
host = "perlmutter.nersc.gov"
remote_dir = "/global/homes/s/soham_b/tsadar_copy/tsadar/external/data/"
file_pattern = "*92534*"          # examples: "*.h5", "*.png", "*.txt"
local_dir = "/Users/soham/Documents/nersc"
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