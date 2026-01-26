#!/usr/bin/env python3
"""
Sanity plots for CCM gamma simulation H5:
Truth:
  1) Primary gamma energy (GeV)
  2) Primary gamma vertex x/y/z (m)
Detector-physics:
  3) Heatmap (2D hist): TotalDepositedEnergyLAr vs primary gamma vertex z

cmd:
  cd ~/workspaces/CCM/sources/CCMAnalysis
  python3 -m venv venv_plot_gamma
  source venv_plot_gamma/bin/activate
  python -m pip install --upgrade pip
  python -m pip install matplotlib h5py numpy

  python plot_gamma.py gamma_100evt.h5

  deactivate
"""

import os
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


def outdir_from_input(fn: str) -> str:
    base = os.path.basename(fn)
    stem = base[:-3] if base.endswith(".h5") else os.path.splitext(base)[0]
    outdir = os.path.join("plots", stem)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def read_primary_gamma_per_event(h5: h5py.File):
    """
    Robust: get one primary gamma per event using the ragged-table index.
    Returns arrays: E_primary, x, y, z, zenith, azimuth (length = n_events where found)
    """
    mct = h5["I3MCTree"]
    idx = h5["__I3Index__/I3MCTree"]

    n_evt = idx.shape[0]
    E_list, x_list, y_list, z_list, zen_list, azi_list = [], [], [], [], [], []

    for i in range(n_evt):
        if idx["exists"][i] != 1:
            continue
        start = int(idx["start"][i])
        stop  = int(idx["stop"][i])
        rows = mct[start:stop]
        if rows.size == 0:
            continue

        # Primary gamma: exists==1, type==22, depth==0
        sel = (rows["exists"] == 1) & (rows["type"] == 22) & (rows["depth"] == 0)
        if not np.any(sel):
            continue

        r0 = rows[sel][0]  # take the first one if multiple
        E_list.append(float(r0["energy"]))
        x_list.append(float(r0["x"]))
        y_list.append(float(r0["y"]))
        z_list.append(float(r0["z"]))
        zen_list.append(float(r0["zenith"]))
        azi_list.append(float(r0["azimuth"]))

    return (np.array(E_list), np.array(x_list), np.array(y_list),
            np.array(z_list), np.array(zen_list), np.array(azi_list))


def read_edep_per_event(h5: h5py.File):
    """
    TotalDepositedEnergyLAr is per event (not ragged), but still has exists flag.
    Returns edep array (GeV) aligned by event row index; also returns a mask exists.
    """
    d = h5["TotalDepositedEnergyLAr"][:]
    exists = (d["exists"] == 1)
    edep = d["value"].astype(float)
    return edep, exists


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else "gamma_100evt.h5"
    outdir = outdir_from_input(fn)

    with h5py.File(fn, "r") as f:
        Eprim, x, y, z, zen, azi = read_primary_gamma_per_event(f)
        edep, edep_exists = read_edep_per_event(f)

    if Eprim.size == 0:
        raise RuntimeError("No primary gammas found (type==22, depth==0).")

    # -------------------------
    # Plot 1: primary energy
    # -------------------------
    plt.figure()
    plt.hist(Eprim, bins=50)
    plt.xlabel("Primary gamma energy [GeV]")
    plt.ylabel("Counts")
    plt.title("Primary gamma energy (truth)")
    plt.grid(True)
    plt.tight_layout()
    p1 = os.path.join(outdir, "truth_primary_gamma_energy.png")
    plt.savefig(p1, dpi=150)
    # -----------------------------------------
    # Plot 2: primary vertex x/y/z distributions
    # -----------------------------------------
    plt.figure()
    plt.hist(x, bins=50, alpha=0.7, label="x")
    plt.hist(y, bins=50, alpha=0.7, label="y")
    plt.hist(z, bins=50, alpha=0.7, label="z")
    plt.xlabel("Primary gamma vertex coordinate [m]")
    plt.ylabel("Counts")
    plt.title("Primary gamma vertex (truth)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    p2 = os.path.join(outdir, "truth_primary_gamma_vertex_xyz.png")
    plt.savefig(p2, dpi=150)
    # ---------------------------------------------------------
    # Plot 3: HEATMAP Edep vs primary vertex z
    # edep is stored per event row. primary arrays are per-event too,
    # but extracted only events where primary gamma exists.
    #
    # For strict alignment, map by (Run,Event,SubEvent), but in
    # sim output it’s typically 1-to-1 with event row. This plot is for sanity checks.
    #
    # Easiest: use same length by truncating to min length where both exist.
   # Plot 3: Heatmap (2D hist) Edep vs primary vertex z
    # ---------------------------------------------------------
    n = min(len(z), np.count_nonzero(edep_exists))
    if n > 0:
        z_use = z[:n]
        edep_use = edep[edep_exists][:n]

        plt.figure()
        plt.hist2d(z_use, edep_use, bins=(50, 50))
        plt.xlabel("Primary gamma vertex z [m]")
        plt.ylabel("Total deposited energy in LAr [GeV]")
        plt.title("Heatmap: Edep(LAr) vs primary vertex z")
        plt.colorbar(label="Counts")
        plt.tight_layout()
        p3 = os.path.join(outdir, "heatmap_edep_vs_vertex_z.png")
        plt.savefig(p3, dpi=150)
        print(f"Saved plot to {p3}")


if __name__ == "__main__":
    main()

