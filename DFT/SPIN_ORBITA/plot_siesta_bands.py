"""
plot_siesta_bands.py
--------------------
Parse and plot SIESTA .bands output files.

Features:
  - Automatic spin-polarised detection (spin weight = 2)
  - Fermi-level shift to 0 eV (optional, on by default)
  - Band-gap identification (optional)
  - High-symmetry k-point labels read from the file footer
  - Per-spin plot or combined plot for spin-polarised calculations

Usage examples
--------------
# Basic plot (Fermi level shifted to 0):
    python plot_siesta_bands.py my_calc.bands

# Keep absolute energies (no Fermi shift):
    python plot_siesta_bands.py my_calc.bands --no-shift

# Enable gap detection:
    python plot_siesta_bands.py my_calc.bands --gap

# Spin-polarised: plot only spin-up channel:
    python plot_siesta_bands.py my_calc.bands --spin up

# Spin-polarised: plot both spins side by side:
    python plot_siesta_bands.py my_calc.bands --spin both

# Set energy window around Fermi level:
    python plot_siesta_bands.py my_calc.bands --emin -10 --emax 10

# Save figure instead of showing it:
    python plot_siesta_bands.py my_calc.bands --save bands.png
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_bands_file(filepath):
    """
    Parse a SIESTA .bands file.

    Returns
    -------
    dict with keys:
        fermi        : float  – Fermi energy (eV)
        kmin, kmax   : float  – k-path range (Å⁻¹)
        emin, emax   : float  – energy range (eV)
        n_bands      : int    – number of bands per spin channel
        spin_weight  : int    – 1 (non-polarised) or 2 (spin-polarised)
        n_kpoints    : int    – number of k-points
        kpoints      : ndarray (n_kpoints,)
        bands        : ndarray
                       shape (n_kpoints, n_bands)        if spin_weight == 1
                       shape (2, n_kpoints, n_bands)     if spin_weight == 2
        klabels      : list of (k_position, label) tuples  (may be empty)
    """
    with open(filepath, "r") as fh:
        lines = fh.readlines()

    idx = 0

    # Line 1 – Fermi level
    fermi = float(lines[idx].split()[0]); idx += 1

    # Line 2 – k-path range
    kmin, kmax = map(float, lines[idx].split()); idx += 1

    # Line 3 – energy range
    emin_file, emax_file = map(float, lines[idx].split()); idx += 1

    # Line 4 – n_bands, spin_weight, n_kpoints
    n_bands, spin_weight, n_kpoints = map(int, lines[idx].split()); idx += 1

    # -----------------------------------------------------------------------
    # Read band data
    # n_bands   = number of bands PER SPIN CHANNEL
    # spin_weight = 1 (non-polarised) or 2 (spin-polarised)
    # Total energy values per k-point = n_bands * spin_weight
    #   → non-polarised: k  e1 e2 ... e_n_bands
    #   → spin-polarised: k  e1_up...e_n_up  e1_dn...e_n_dn
    # -----------------------------------------------------------------------
    total_per_k = n_bands * spin_weight
    kpoints     = np.empty(n_kpoints)
    all_energies = np.empty((n_kpoints, total_per_k))

    i = idx
    for ik in range(n_kpoints):
        tokens = []
        while len(tokens) < total_per_k + 1:
            tokens += lines[i].split()
            i += 1
        kpoints[ik]       = float(tokens[0])
        all_energies[ik]  = list(map(float, tokens[1:total_per_k + 1]))
    idx = i

    if spin_weight == 1:
        bands = all_energies                          # shape (nk, n_bands)
    else:
        bands = np.array([
            all_energies[:, :n_bands],                # spin-up   (nk, n_bands)
            all_energies[:, n_bands:],                # spin-down (nk, n_bands)
        ])                                            # shape (2, nk, n_bands)

    # -----------------------------------------------------------------------
    # Read high-symmetry k-point labels (footer)
    # Format:  <n_labels>
    #          <k_pos>   '<Label>'
    # -----------------------------------------------------------------------
    klabels = []
    try:
        n_labels = int(lines[idx].strip())
        idx += 1
        for _ in range(n_labels):
            parts = lines[idx].split()
            k_pos = float(parts[0])
            label = parts[1].strip("'\"")
            # Convert Gamma → proper symbol
            if label.lower() in ("gamma", "g"):
                label = r"$\Gamma$"
            klabels.append((k_pos, label))
            idx += 1
    except (IndexError, ValueError):
        pass  # footer absent or malformed – silently ignore

    return dict(
        fermi=fermi,
        kmin=kmin, kmax=kmax,
        emin_file=emin_file, emax_file=emax_file,
        n_bands=n_bands,
        spin_weight=spin_weight,
        n_kpoints=n_kpoints,
        kpoints=kpoints,
        bands=bands,
        klabels=klabels,
    )


# ---------------------------------------------------------------------------
# Gap finder
# ---------------------------------------------------------------------------

def find_gap(kpoints, bands_2d, fermi_ref=0.0):
    """
    Find the band gap for one spin channel.

    Parameters
    ----------
    bands_2d  : ndarray (n_kpoints, n_bands) – energies already Fermi-shifted
    fermi_ref : float – reference Fermi level in the shifted frame (usually 0)

    Returns
    -------
    dict with keys:
        type       : 'metallic' | 'direct' | 'indirect'
        gap        : float (eV)  – 0 for metallic
        vbm        : float (eV)
        cbm        : float (eV)
        vbm_k      : float – k-position of VBM
        cbm_k      : float – k-position of CBM
    """
    # Occupied bands: highest energy <= fermi_ref at each k
    # Unoccupied bands: lowest energy > fermi_ref at each k
    vbm_per_k = np.max(bands_2d[bands_2d <= fermi_ref].reshape(
        bands_2d.shape[0], -1
    ) if False else  # placeholder – computed properly below
        np.where(bands_2d <= fermi_ref, bands_2d, -np.inf), axis=1)

    cbm_per_k = np.min(
        np.where(bands_2d > fermi_ref, bands_2d, np.inf), axis=1)

    vbm = np.max(vbm_per_k)
    cbm = np.min(cbm_per_k)

    if cbm <= vbm:
        return dict(type="metallic", gap=0.0,
                    vbm=vbm, cbm=cbm, vbm_k=None, cbm_k=None)

    gap = cbm - vbm
    vbm_k = kpoints[np.argmax(vbm_per_k)]
    cbm_k = kpoints[np.argmin(cbm_per_k)]

    gap_type = "direct" if np.isclose(vbm_k, cbm_k, atol=1e-4) else "indirect"

    return dict(type=gap_type, gap=gap,
                vbm=vbm, cbm=cbm, vbm_k=vbm_k, cbm_k=cbm_k)


# ---------------------------------------------------------------------------
# Plotter
# ---------------------------------------------------------------------------

SPIN_COLORS = {
    "up":   "#d62728",   # red
    "down": "#1f77b4",   # blue
    "none": "#2c2c2c",   # near-black for non-polarised
}

def _plot_one_channel(ax, kpoints, bands_2d, color, label=None, lw=0.8, alpha=1.0):
    """Plot all bands for one spin channel onto ax."""
    for ib in range(bands_2d.shape[1]):
        ax.plot(kpoints, bands_2d[:, ib],
                color=color, lw=lw, alpha=alpha,
                label=label if ib == 0 else None)


def _add_klabels(ax, klabels):
    """Draw vertical lines and tick labels at high-symmetry points."""
    if not klabels:
        return
    positions = [k for k, _ in klabels]
    labels    = [l for _, l in klabels]
    for pos in positions:
        ax.axvline(pos, color="black", lw=0.8, ls="-")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=12)


def _annotate_gap(ax, gap_info, kpoints):
    """Draw VBM/CBM markers and a gap bracket."""
    if gap_info["type"] == "metallic":
        ax.set_title(ax.get_title() + "  [Metallic]", fontsize=11)
        return

    vbm, cbm = gap_info["vbm"], gap_info["cbm"]
    vbm_k, cbm_k = gap_info["vbm_k"], gap_info["cbm_k"]

    ax.axhline(vbm, color="#e67e22", lw=1.0, ls="--", alpha=0.8)
    ax.axhline(cbm, color="#27ae60", lw=1.0, ls="--", alpha=0.8)

    ax.scatter([vbm_k], [vbm], color="#e67e22", s=40, zorder=5)
    ax.scatter([cbm_k], [cbm], color="#27ae60", s=40, zorder=5)

    mid_k = (kpoints[0] + kpoints[-1]) / 2
    ax.annotate(
        "", xy=(mid_k, cbm), xytext=(mid_k, vbm),
        arrowprops=dict(arrowstyle="<->", color="purple", lw=1.5)
    )
    ax.text(mid_k, (vbm + cbm) / 2,
            f"  {gap_info['gap']:.3f} eV\n  ({gap_info['type']})",
            color="purple", fontsize=9, va="center")

    title_suffix = (f"Gap = {gap_info['gap']:.3f} eV  [{gap_info['type']}]  "
                    f"VBM={vbm:.3f} eV  CBM={cbm:.3f} eV")
    ax.set_title(ax.get_title() + "\n" + title_suffix, fontsize=9)


def plot_bands(data, shift_fermi=True, show_gap=False,
               spin_choice="both", emin=None, emax=None, save=None):
    """
    Main plotting routine.

    Parameters
    ----------
    data         : dict returned by parse_bands_file
    shift_fermi  : bool  – shift energies so E_F = 0
    show_gap     : bool  – detect and annotate the band gap
    spin_choice  : 'up' | 'down' | 'both' | 'all'
                   For non-polarised calculations this parameter is ignored.
    emin, emax   : float or None – energy window for the plot
    save         : str or None – filename to save the figure
    """
    fermi       = data["fermi"]
    kpoints     = data["kpoints"]
    bands_raw   = data["bands"]
    spin_weight = data["spin_weight"]
    klabels     = data["klabels"]

    shift = fermi if shift_fermi else 0.0
    fermi_ref = 0.0 if shift_fermi else fermi

    # ---- determine energy window ----
    e_lo = emin if emin is not None else -20.0
    e_hi = emax if emax is not None else  20.0

    # -----------------------------------------------------------------------
    # Non-polarised case
    # -----------------------------------------------------------------------
    if spin_weight == 1:
        bands = bands_raw - shift        # (nk, nb)

        fig, ax = plt.subplots(figsize=(7, 6))
        _plot_one_channel(ax, kpoints, bands, color=SPIN_COLORS["none"])

        ax.axhline(fermi_ref, color="gray", lw=1.0, ls="--", alpha=0.7,
                   label=f"E$_F$ = {fermi:.4f} eV" if not shift_fermi else "E$_F$")

        if show_gap:
            gap_info = find_gap(kpoints, bands, fermi_ref)
            _annotate_gap(ax, gap_info, kpoints)
            print(_gap_summary("", gap_info))

        _add_klabels(ax, klabels)
        ax.set_xlim(kpoints[0], kpoints[-1])
        ax.set_ylim(e_lo, e_hi)
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_xlabel("k-path", fontsize=12)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set_title("SIESTA Band Structure", fontsize=13)
        ax.legend(fontsize=9, loc="upper right")
        fig.tight_layout()

    # -----------------------------------------------------------------------
    # Spin-polarised case
    # -----------------------------------------------------------------------
    else:
        bands_up = bands_raw[0] - shift   # (nk, nb)
        bands_dn = bands_raw[1] - shift

        spin_choice = spin_choice.lower()
        fig, ax = plt.subplots(figsize=(7, 6))

        if spin_choice == "up":
            channels = [(bands_up, "Spin ↑", SPIN_COLORS["up"])]
        elif spin_choice in ("down", "dn"):
            channels = [(bands_dn, "Spin ↓", SPIN_COLORS["down"])]
        else:  # "both" or "all" – overlay on the same axes
            channels = [
                (bands_up, "Spin ↑", SPIN_COLORS["up"]),
                (bands_dn, "Spin ↓", SPIN_COLORS["down"]),
            ]

        for bnd, lbl, color in channels:
            _plot_one_channel(ax, kpoints, bnd, color=color, label=lbl, alpha=0.75)
            if show_gap:
                gap_info = find_gap(kpoints, bnd, fermi_ref)
                _annotate_gap(ax, gap_info, kpoints)
                print(_gap_summary(lbl, gap_info))

        ax.axhline(fermi_ref, color="gray", lw=1.0, ls="--", alpha=0.7, label="E$_F$")
        _add_klabels(ax, klabels)
        ax.set_xlim(kpoints[0], kpoints[-1])
        ax.set_ylim(e_lo, e_hi)
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_xlabel("k-path", fontsize=12)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.legend(fontsize=9, loc="upper right")
        ax.set_title("SIESTA Band Structure (spin-polarised)", fontsize=13)
        fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save}")
    else:
        plt.show()


def _gap_summary(label, gap_info):
    prefix = f"[{label}] " if label else ""
    if gap_info["type"] == "metallic":
        return f"{prefix}System is METALLIC (bands cross E_F)."
    return (f"{prefix}Band gap: {gap_info['gap']:.4f} eV  "
            f"({gap_info['type']})  |  "
            f"VBM = {gap_info['vbm']:.4f} eV at k = {gap_info['vbm_k']:.5f}  |  "
            f"CBM = {gap_info['cbm']:.4f} eV at k = {gap_info['cbm_k']:.5f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot SIESTA band structure from a .bands file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("bands_file", help="Path to the .bands file")
    parser.add_argument("--no-shift", action="store_true",
                        help="Do NOT shift energies to E_F = 0 (plot absolute energies)")
    parser.add_argument("--gap", action="store_true",
                        help="Detect and annotate the band gap")
    parser.add_argument("--spin", default="both",
                        choices=["up", "down", "dn", "both", "all"],
                        help="For spin-polarised: which channel(s) to plot "
                             "[up | down | both (default)]")
    parser.add_argument("--emin", type=float, default=None,
                        help="Lower energy bound for plot (eV, relative to E_F if shifted)")
    parser.add_argument("--emax", type=float, default=None,
                        help="Upper energy bound for plot (eV, relative to E_F if shifted)")
    parser.add_argument("--save", type=str, default=None,
                        help="Save figure to this filename (e.g. bands.png) instead of showing")

    args = parser.parse_args()

    print(f"Reading: {args.bands_file}")
    data = parse_bands_file(args.bands_file)

    print(f"  Fermi level   : {data['fermi']:.6f} eV")
    print(f"  Spin channels : {data['spin_weight']}")
    print(f"  Bands per spin: {data['n_bands']}"
          + (f"  (total {data['n_bands'] * data['spin_weight']} per k-point)" if data['spin_weight'] == 2 else ""))
    print(f"  k-points      : {data['n_kpoints']}")
    if data["klabels"]:
        kl_str = "  →  ".join(f"{l} ({k:.4f})" for k, l in data["klabels"])
        print(f"  k-path        : {kl_str}")

    plot_bands(
        data,
        shift_fermi=not args.no_shift,
        show_gap=args.gap,
        spin_choice=args.spin,
        emin=args.emin,
        emax=args.emax,
        save=args.save,
    )


if __name__ == "__main__":
    main()
