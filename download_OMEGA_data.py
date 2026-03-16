from pathlib import Path
import re
import csv
import requests
from bs4 import BeautifulSoup
import os

# =========================
# User settings
# =========================
USERNAME = os.environ["OMEGA_USER"]
PASSWORD = os.environ["OMEGA_PASSWORD"]

SHOTNUMBER = "117828"   # The shot data you want
DIAG = ["EPW","IAW","P9TBD_CCD","XRPHC-CID","SRS_STREAK","SBS_STREAK","TPDI"]   # The diagnostic images you want
PAGE_URL = (
    "https://omegaops.lle.rochester.edu/lir"
    f"?goSingle=goSingle&singleReport=Admin_Summary&shotnumber={SHOTNUMBER}"
)   # This is the url for the Omega shot images and reports
Download_url = "https://omegaops.lle.rochester.edu/lirdir/archiveDownload.py"  # This is the url for the download button
baseDir = Path("/Users/soham/Documents/Kinshock/Kinshock-26A/data")   # where you want to save teh data
shot = Path(f"{SHOTNUMBER}")
shot = baseDir / shot
if not shot.exists():
    shot.mkdir(parents = True, exist_ok = True)
    
OUTPUT_CSV = f"shot_{SHOTNUMBER}_file_lookup.csv"
SAVE_HTML_COPY = False


# =========================
# Fetch page - see if the OMEGA images page can be accessed 
# =========================
session = requests.Session()
session.auth = (USERNAME, PASSWORD) # This will enter the username and password
session.headers.update({"User-Agent": "Mozilla/5.0"})

response = session.get(PAGE_URL, timeout=60)
response.raise_for_status()

html = response.text  # This is html text for the page

if SAVE_HTML_COPY:  # Keep false unless you want to see and understand the parsing
    Path(f"shot_{SHOTNUMBER}_page.html").write_text(html, encoding="utf-8")

print(f"Fetched page successfully: status = {response.status_code}")


# =========================
# Parse page - Look through the page with BeautifulSoup and find the download/ image IDs for your desired diagnostic images 
# =========================
soup = BeautifulSoup(html, "html.parser")

records = []  # This is the list of dictionaries with all the parsed html information

for a in soup.find_all("a"):
    filename = a.get_text(strip=True)
    onclick = a.get("onclick", "").strip()

    # We only want file entries that call DoImage(...)
    if not onclick.startswith("DoImage("):
        continue

    # Pull out everything inside DoImage(...)
    m = re.search(r"DoImage\((.*)\)\s*;?\s*$", onclick)
    if not m:
        continue

    args_str = m.group(1)

    # Extract all quoted arguments
    quoted = re.findall(r"'([^']*)'", args_str)

    # Based on the page structure, these quoted args typically look like:
    # 0: shotnumber
    # 1: date string
    # 2: srcfile
    # 3: archive_file
    # 4: format
    # 5: diagType
    # 6: viewer
    # ...
    # last: image_id
    if len(quoted) < 7:
        continue

    image_id = quoted[-1]
    if not image_id.isdigit():
        continue
    
    # save as a dictionary
    record = {
        "shotnumber": quoted[0] if len(quoted) > 0 else "",
        "date": quoted[1] if len(quoted) > 1 else "",
        "srcfile": quoted[2] if len(quoted) > 2 else "",
        "archive_file": quoted[3] if len(quoted) > 3 else "",
        "format": quoted[4] if len(quoted) > 4 else "",
        "diag": quoted[5] if len(quoted) > 5 else "",
        "viewer": quoted[6] if len(quoted) > 6 else "",
        "filename": filename,
        "image_id": image_id,
        "onclick": onclick,
    }

    # append to make list of dictionaries
    records.append(record)

print(f"Parsed {len(records)} file records")


# =========================
# Optional: remove exact duplicates
# =========================
unique_records = []
seen = set()

for r in records:
    key = (
        r["filename"],
        r["diag"],
        r["image_id"],
        r["srcfile"],
        r["archive_file"],
    )
    if key not in seen:
        seen.add(key)
        unique_records.append(r)

records = unique_records
print(f"{len(records)} unique records after deduplication")


# =========================
# Save to CSV
# =========================
fieldnames = [
    "shotnumber",
    "date",
    "diag",
    "filename",
    "image_id",
    "format",
    "viewer",
    "srcfile",
    "archive_file",
    "onclick",
]


with open(shot / OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(records)

print(f"Saved CSV: {OUTPUT_CSV}")

# =========================
# Download desired files
# =========================

for diag in DIAG: # step through the diagnostics you want

    matches = [r for r in records if r["diag"] == diag]

    for match in matches:
        image_id = match["image_id"]
        params = {
            "convert": "y",
            "image_id": image_id,
        }

        try:  # This time you need to use the Download_url
            with session.get(Download_url, params=params, stream=True, timeout=60) as r:
                print(f"{diag}: status {r.status_code}")
                r.raise_for_status()

                content_disp = r.headers.get("Content-Disposition")
                filename = f"{image_id}.h5"
                if content_disp and "filename=" in content_disp:  # take out the file name from the html headers
                    filename = content_disp.split("filename=")[-1].strip('"')
                
                out_path = shot / filename

                with open(out_path, "wb") as f:  # download the file by creating and writing it in your out_path
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

                print(f"Downloaded: {filename}")

        except requests.exceptions.RequestException as e:
            print(f"Failed for image_id={image_id}: {e}")
