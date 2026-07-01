import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LogNorm


float_re = re.compile(r"[-+]?\d*\.\d+E[-+]?\d+|[-+]?\d+E[-+]?\d+|[-+]?\d+\.\d+|[-+]?\d+")

class PropaceosTables:
    def __init__(self, filename):
        print("Creating PROPACEOS table")

        self.filename = filename

    def read(self):
        print(f"Reading {self.filename}")

        with open(self.filename) as f:
            self.lines = f.readlines()
        
        print(f"Read {len(self.lines)} lines")
    
    def find_line(self, text):
        text = text.lower()
        for i, line in enumerate(self.lines):
            if text in line.lower():
                return i
        raise ValueError(f"Could not find: {text}")
    
    def nums(self, line):
        return [float(x) for x in float_re.findall(line)]
    
    def read_n_values(self, start, n):
        vals = []
        i = start

        while len(vals) < n and i < len(self.lines):
            line = self.lines[i]

            if "****" not in line:
                vals.extend(self.nums(line))

            i += 1

        if len(vals) < n:
            raise ValueError(f"Only found {len(vals)} values, expected {n}")

        return np.array(vals[:n]), i
    
    def read_mesh(self):
        i = self.find_line("mesh parameters for EoS")

        nT = int(self.nums(self.lines[i+5])[0])
        T, j = self.read_n_values(i+6, nT)

        nD = int(self.nums(self.lines[j])[0])
        D, j = self.read_n_values(j+1, nD)

        self.nT = nT
        self.nD = nD
        self.T = T
        self.ni = D

        print(f"Grid size [nT, nD] = [{self.nT}, {self.nD}]")

    def read_quantity(self, name, attr_name=None):
        if attr_name is None:
            attr_name = name

        i = self.find_line(name)
        vals, j = self.read_n_values(i + 1, self.nT * self.nD)
        arr = vals.reshape(self.nD, self.nT)

        setattr(self, attr_name, arr)

        print(f"Read {name} -> self.{attr_name}: shape = {arr.shape}")
        return arr
    
    def load_all(self):
        self.read()
        self.read_mesh()

        quantities = {
            "Zbar": "Zbar",
            "Int. Rosseland Mean Opacity": "rosseland_int",
            "Int. emis. Planck Mean Opacity": "planck_emis_int",
            "Int. abs. Planck Mean Opacity": "planck_abs_int",
            "Eint": "Eint",
            "Eion": "Eion",
            "Eele": "Eele",
            "Pion": "Pion",
            "Pele": "Pele",
        }

        for name, attr in quantities.items():
            try:
                self.read_quantity(name, attr)
            except ValueError:
                print(f"Skipping {name}: not found")
                
    def plot(self, quantity, logz=True):
        T = self.T
        ni = self.ni
        if not hasattr(self, quantity):
            raise ValueError(f"{quantity} has not been read.")
        
        Z = getattr(self, quantity)

        TT, NN = np.meshgrid(T, ni)

        fig, ax = plt.subplots(figsize=(7,5))

        if logz:
            Zplot = np.where(Z > 0, Z, np.nan)
            pcm = ax.pcolormesh(TT, NN, Zplot, shading="auto", norm=LogNorm())
        else:
            pcm = ax.pcolormesh(TT, NN, Z, shading="auto")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Temperature [eV]")
        ax.set_ylabel(r"Ion density [cm$^{-3}$]")
        ax.set_title(quantity)

        cbar = fig.colorbar(pcm, ax=ax)
        cbar.set_label(quantity)

        plt.tight_layout()
        return fig, ax
    
    from pathlib import Path

    def strip_for_flash(self, outpath):
        """
        1. Print FLASH rt_mgdNumGroups / rt_mgdBounds from group structure.
        2. Remove only the Ionization Fractions data section.
        3. Write stripped file to outpath.
        """

        # ---- Read opacity group bounds ----
        bounds = []
        reading_groups = False

        for line in self.lines:
            s = line.strip()
            low = s.lower()

            if "group structure" in low and "(ev)" in low:
                reading_groups = True
                continue

            if reading_groups:
                # stop only at the section header, not header metadata
                if "********************" in line and "ionization fractions" in low:
                    break

                if not s or "*" in s:
                    continue

                for token in s.replace("D", "E").split():
                    try:
                        bounds.append(float(token))
                    except ValueError:
                        pass

        num_groups = len(bounds) - 1

        print(f"rt_mgdNumGroups = {num_groups}")
        for i, b in enumerate(bounds, start=1):
            print(f"rt_mgdBounds_{i:<2} = {b:.6e}")

        # ---- Strip only Ionization Fractions section ----
        stripped_lines = []
        skipping = False

        for line in self.lines:
            low = line.lower()

            ion_frac_section = (
                "********************" in line
                and "ionization fractions" in low
            )

            zbar_section = (
                "********************" in line
                and "zbar" in low
            )

            if ion_frac_section:
                skipping = True
                continue

            if skipping and zbar_section:
                skipping = False
                stripped_lines.append(line)
                continue

            if skipping:
                continue

            stripped_lines.append(line)

        self.lines = stripped_lines

        with open(outpath, "w") as f:
            for line in self.lines:
                f.write(line if line.endswith("\n") else line + "\n")

        print(f"\nWrote stripped FLASH table to: {outpath}")
        return self.lines
    


# usage
fp = Path(r"C:\Simulation_data\PROPACEOS\gas_6_EnergyBins\gas_6_EnergyBins.prp")
table = PropaceosTables(fp)
table.load_all()
fig, ax = table.plot("rosseland_int", logz = False)
plt.show()
op = Path(r"C:\Simulation_data\PROPACEOS\gas_6_EnergyBins\gas_6_EnergyBins_stripped.prp")
table.strip_for_flash(op)
