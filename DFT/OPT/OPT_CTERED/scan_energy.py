import os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import shutil

template_file = "mos2.fdf"

a_values = np.arange(3.16, 3.19, 0.002)

siesta_command = "siesta < input.fdf > output.out"

# ==========================================================
# MODIFICAR LATTICE CONSTANT
# ==========================================================
def modify_latticeconstant(content, a):
    pattern = r"LatticeConstant\s+.*"
    replacement = f"LatticeConstant {a:.4f} Ang"
    return re.sub(pattern, replacement, content)

# ==========================================================
# EXTRAER ENERGÍA
# ==========================================================
def extract_energy(outfile):
    energy = None
    with open(outfile) as f:
        for line in f:
            if "siesta: E_KS(eV)" in line:
                try:
                    energy = float(line.split()[-1])
                except:
                    pass
    return energy

# ==========================================================
# MAIN
# ==========================================================
with open(template_file) as f:
    original_content = f.read()

energies = []

for a in a_values:

    dirname = f"a_{a:.3f}"
    os.makedirs(dirname, exist_ok=True)

    # input.fdf
    modified_content = modify_latticeconstant(original_content, a)
    with open(os.path.join(dirname, "input.fdf"), "w") as f:
        f.write(modified_content)

    # ======================================================
    # COPIAR PSEUDOPOTENCIALES (PSML)
    # ======================================================
    for p in ["Mo.psml", "S.psml"]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Falta {p} en el directorio base")
        shutil.copy(p, dirname)

    # ======================================================
    # EJECUTAR SIESTA
    # ======================================================
    print(f"Ejecutando a = {a:.3f}")

    result = subprocess.run(
        siesta_command,
        shell=True,
        cwd=dirname,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        print("ERROR en SIESTA para a =", a)
        print(result.stderr)

    # ======================================================
    # LEER ENERGÍA
    # ======================================================
    energy = extract_energy(os.path.join(dirname, "output.out"))
    energies.append(energy)

    print(f"E = {energy} eV")

# ==========================================================
# RESULTADOS
# ==========================================================
energies = np.array(energies, dtype=float)

valid = ~np.isnan(energies)

a_valid = a_values[valid]
e_valid = energies[valid]

idx = np.argmin(e_valid)

print("\n==========================")
print(f"a óptimo = {a_valid[idx]:.4f} Å")
print(f"E mínima = {e_valid[idx]:.6f} eV")
print("==========================")

# ==========================================================
# GRÁFICA
# ==========================================================
plt.figure(figsize=(8,6))
plt.plot(a_valid, e_valid, "o-")
plt.xlabel("Lattice constant a (Å)")
plt.ylabel("Energía total (eV)")
plt.grid()
plt.tight_layout()
plt.savefig("energy_curve.png", dpi=300)
plt.show()
