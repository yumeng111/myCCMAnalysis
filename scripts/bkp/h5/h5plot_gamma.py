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
    One primary gamma per event using ragged index.
    Primary gamma definition: exists==1, type==22, depth==0
    Returns: E [GeV], x,y,z [m]
    """
    mct = h5["I3MCTree"]
    idx = h5["__I3Index__/I3MCTree"]

    E, x, y, z = [], [], [], []
    for i in range(idx.shape[0]):
        if idx["exists"][i] != 1:
            continue
        start = int(idx["start"][i])
        stop = int(idx["stop"][i])
        rows = mct[start:stop]
        if rows.size == 0:
            continue

        sel = (rows["exists"] == 1) & (rows["type"] == 22) & (rows["depth"] == 0)
        if not np.any(sel):
            continue

        r0 = rows[sel][0]
        E.append(float(r0["energy"]))
        x.append(float(r0["x"]))
        y.append(float(r0["y"]))
        z.append(float(r0["z"]))

    return np.array(E), np.array(x), np.array(y), np.array(z)


def read_total_edep(h5: h5py.File):
    """
    TotalDepositedEnergyLAr per event with exists flag.
    Returns: edep [GeV], mask_exists [bool]
    """
    d = h5["TotalDepositedEnergyLAr"][:]
    exists = (d["exists"] == 1)
    edep = d["value"].astype(float)
    return edep, exists


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else "gamma_100evt.h5"
    outdir = outdir_from_input(fn)

    with h5py.File(fn, "r") as f:
        Eprim_GeV, x, y, z = read_primary_gamma_per_event(f)
        edep_GeV, edep_exists = read_total_edep(f)

    # Convert to MeV
    Eprim_MeV = 1e3 * Eprim_GeV
    edep_MeV = 1e3 * edep_GeV

    # -------------------------
    # Plot 1: primary energy (MeV)
    plt.figure()
    plt.hist(Eprim_MeV, bins=50)
    plt.xlabel("Primary gamma energy [MeV]")
    plt.ylabel("Counts")
    plt.title("Primary gamma energy (truth)")
    plt.grid(True)
    plt.tight_layout()
    p1 = os.path.join(outdir, "truth_primary_gamma_energy.png")
    plt.savefig(p1, dpi=150)

    # -------------------------
    # Plot 2: primary vertex xyz (m)
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

    # -------------------------
    # Plot 3: heatmap Edep vs vertex z (MeV)
    n = min(len(z), np.count_nonzero(edep_exists))
    if n > 0:
        z_use = z[:n]
        edep_use_MeV = edep_MeV[edep_exists][:n]

        plt.figure()
        plt.hist2d(z_use, edep_use_MeV, bins=(50, 50))
        plt.xlabel("Primary gamma vertex z [m]")
        plt.ylabel("Total deposited energy [MeV]")
        plt.title("Edep vs primary vertex z")
        plt.colorbar(label="Counts")
        plt.tight_layout()
        p3 = os.path.join(outdir, "heatmap_edep_vs_vertex_z.png")
        plt.savefig(p3, dpi=150)

    # -------------------------
    # Plot 4: detector deposited energy distribution (MeV)
    edep_only = edep_MeV[edep_exists]
    if edep_only.size > 0:
        plt.figure()
        plt.hist(edep_only, bins=50)
        plt.xlabel("Total deposited energy [MeV]")
        plt.ylabel("Events")
        plt.title("Total deposited energy")
        plt.grid(True)
        plt.tight_layout()
        p4 = os.path.join(outdir, "detector_edep.png")
        plt.savefig(p4, dpi=150)

if __name__ == "__main__":
    main()
