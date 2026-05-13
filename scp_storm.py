import subprocess
from pathlib import Path

# -------- user settings --------
username = "sban"
host = "cl4head"
remote_dir = "/home/sban/Si"
file_pattern = "*"
local_dir = "//profiles/Users$/sban/Documents/Kinshock/LILAC_simulations/Si_5"
# -------------------------------

def scp_copy_files(username, host, remote_dir, file_pattern, local_dir):
    """
    Copy matching files from a remote host to a local directory using scp.
    """
    local_path = Path(local_dir).expanduser()
    local_path.mkdir(parents=True, exist_ok=True)

    remote_path = f'{username}@{host}:{remote_dir.rstrip("/")}/{file_pattern}'
    cmd = ["scp","-r", remote_path, str(local_path)]

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