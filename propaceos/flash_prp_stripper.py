from pathlib import Path


def strip_propaceos_ion_fractions(infile, outfile=None):
    infile = Path(infile)

    if outfile is None:
        outfile = infile.with_name(infile.stem + "_stripped" + infile.suffix)
    else:
        outfile = Path(outfile)

    lines = infile.read_text(errors="replace").splitlines(keepends=True)

    start = None
    end = None

    for i, line in enumerate(lines):
        if "Ionization Fractions" in line:
            start = i
        if start is not None and "Zbar" in line:
            end = i
            break

    if start is None:
        raise ValueError("Could not find 'Ionization Fractions' section.")

    if end is None:
        raise ValueError("Could not find 'Zbar' section after Ionization Fractions.")

    # Remove ion fractions block, but keep Zbar and everything after it
    new_lines = lines[:start] + lines[end:]

    outfile.write_text("".join(new_lines))

    print(f"Input file:  {infile}")
    print(f"Output file: {outfile}")
    print(f"Removed lines {start+1} through {end}")
    print(f"Kept Zbar starting at original line {end+1}")

    return outfile


# Example
fp = Path(r"C:\Simulation_data\PROPACEOS\Si3N4_2\Si3N4_2.prp")
op = Path(r"C:\Simulation_data\PROPACEOS\Si3N4_2\Si3N4_2_stripped.prp")
strip_propaceos_ion_fractions(fp)