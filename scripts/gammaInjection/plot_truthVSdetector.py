#   python3 plot_truthVSdetector.py --i3 output/gamma_1000evt_charge.i3.zst --count --out_prefix gammas_2d 

#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from icecube import dataio, icetray


def get_i3double(frame, key):
    if not frame.Has(key):
        return None
    obj = frame[key]
    if hasattr(obj, "value"):
        try:
            return float(obj.value)
        #if value not numeric
        except Exception:
            return None
    #if no .value
    try:
        return float(obj)
    except Exception:
        return None


def sum_pulse_charge(frame, use_per_particle=True):
    total_q = 0.0
    found = False

    # MultiParticle: pid -> pmt -> pulses
    if use_per_particle and frame.Has("MCRecoPulsesMultiParticle"):
        mp = frame["MCRecoPulsesMultiParticle"]
        for _, pmt_map in mp.items():
            for _, pulse_series in pmt_map.items():
                for pulse in pulse_series:
                    total_q += float(pulse.charge)
                    found = True
        return total_q if found else None

    # Merged: pmt -> pulses
    if frame.Has("MCRecoPulses"):
        merged = frame["MCRecoPulses"]
        for _, pulse_series in merged.items():
            for pulse in pulse_series:
                total_q += float(pulse.charge)
                found = True
        return total_q if found else None

    return None


def has_nonempty_pulses(frame, use_per_particle=True):
    """Return True if the frame has a pulses container with at least one pulse."""
    # MultiParticle: pid -> pmt -> pulses
    if use_per_particle and frame.Has("MCRecoPulsesMultiParticle"):
        mp = frame["MCRecoPulsesMultiParticle"]
        for _, pmt_map in mp.items():
            for _, pulse_series in pmt_map.items():
                try:
                    if len(pulse_series) > 0:
                        return True
                except TypeError:
                    # some pulse series may typeerror(not support len()); iterate instead
                    for _ in pulse_series:
                        return True
        return False

    # Merged: pmt -> pulses
    if frame.Has("MCRecoPulses"):
        merged = frame["MCRecoPulses"]
        for _, pulse_series in merged.items():
            try:
                if len(pulse_series) > 0:
                    return True
            except TypeError:
                for _ in pulse_series:
                    return True
        return False

    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--i3", required=True, help="Input .i3/.i3.zst file")
    ap.add_argument("--out_dir", default="output/plots", help="Output directory (default: output/plots)")
    ap.add_argument("--out_prefix", default=None, help="Output filename prefix (default: input basename)")
    ap.add_argument("--bins", type=int, default=120)
    ap.add_argument("--bins1d", type=int, default=100)
    ap.add_argument("--count", action="store_true", help="Compute Q->P counts and efficiency")
    args = ap.parse_args()

    in_path = Path(args.i3)
    base_name = in_path.name
    if base_name.endswith(".i3.zst"):
        base_name = base_name[:-len(".i3.zst")]
    else:
        base_name = in_path.stem

    out_dir = Path(args.out_dir) / base_name
    out_dir.mkdir(parents=True, exist_ok=True)

    out_prefix = args.out_prefix if args.out_prefix else base_name

    n_q_total = 0
    n_q_with_pulses_key = 0
    n_q_nonempty_frames = 0
    q_nonempty_evt = set()     # denominator set: unique Q events with non-empty pulses
    p_evt_all = set()          # numerator set: events that have at least one P frame

    edep_by_evt = {}                # event_id -> truth energy
    qtot_by_evt = {}                # event_id -> total pulse charge from Q frame
    q90sum_by_evt = {}              # event_id -> sum of Charge90NS over P frames
   
    E_truth = []                    # truth energy (stored in GeV)

    GeV_to_MeV = 1000.0
    x_min, x_max = 9.5, 10.5       
    f = dataio.I3File(args.i3)
    try:
        while f.more():
            fr = f.pop_frame()

            # Q frames: counting, truth energy, qtot_by_evt
            if fr.Stop == icetray.I3Frame.DAQ:
                n_q_total += 1

                if fr.Has("MCRecoPulsesMultiParticle") or fr.Has("MCRecoPulses"):
                    n_q_with_pulses_key += 1

                if has_nonempty_pulses(fr, use_per_particle=True):
                    n_q_nonempty_frames += 1
                    if fr.Has("CCMEventHeader"):
                        evt = int(fr["CCMEventHeader"].event_id)
                        q_nonempty_evt.add(evt)

                # Store DAQ pulse charge per event
                if fr.Has("CCMEventHeader"):
                    evt = int(fr["CCMEventHeader"].event_id)
                    e = get_i3double(fr, "TotalDepositedEnergyLAr")
                    if e is not None and evt not in edep_by_evt:
                        edep_by_evt[evt] = float(e)
                        E_truth.append(float(e))
                    
                    qtot = sum_pulse_charge(fr, use_per_particle=True)
                    if qtot is not None:
                        # there should be one DAQ per event; overwrite is fine
                        qtot_by_evt[evt] = float(qtot)
                continue

            #P frames: counting, EventCharge90NS
            if fr.Stop == icetray.I3Frame.Physics:
                if not fr.Has("CCMEventHeader"):
                    continue

                evt = int(fr["CCMEventHeader"].event_id)
                p_evt_all.add(evt)

                q90 = get_i3double(fr, "UncoatedAndCoated8inEventCharge90NS")
                if q90 is not None:
                    q90sum_by_evt[evt] = q90sum_by_evt.get(evt, 0.0) + float(q90)

    finally:
        f.close()

    #efficiency:
    if args.count:
        denom = len(q_nonempty_evt)
        numer = len(q_nonempty_evt & p_evt_all)
        eff = numer / denom if denom > 0 else float("nan")

        print(f"DAQ(Q) frames total: {n_q_total}")
        print(f"DAQ(Q) frames with pulses key present: {n_q_with_pulses_key}")
        print(f"DAQ(Q) frames with non-empty pulses: {n_q_nonempty_frames}")
        print(f"Unique non-empty DAQ events (denominator): {denom}")
        print(f"Unique non-empty DAQ events that have ≥1 Physics frame (numerator): {numer}")
        print(f"Efficiency = (Q non-empty & has P) / (Q non-empty) = {eff:.6f}")

    # Plot 1: 1DTruth deposited energy distribution
    if len(E_truth) > 0:
        E_truth_MeV = np.asarray(E_truth, dtype=float) * GeV_to_MeV
        plt.figure()
        plt.hist(E_truth_MeV, bins=80)
        plt.xlabel("TotalDepositedEnergyLAr(Truth)[MeV]")
        plt.ylabel("Counts")
        plt.title("Truth deposited energy distribution")
        plt.tight_layout()
        plt.savefig(out_dir / f"{out_prefix}_1dtruthE_dist.png", dpi=200)

        #DEBUG
        print("[DEBUG] N =", len(E_truth_MeV))
        print("[DEBUG] min/max E_truth_MeV =", np.min(E_truth_MeV), np.max(E_truth_MeV))
        print("[DEBUG] N(E>15 MeV) =", np.sum(E_truth_MeV > 15))

    # Plot 2: 1D distribution of total pulse charge
    common_evt_pulse = sorted(set(edep_by_evt) & set(qtot_by_evt))
    if len(common_evt_pulse) > 0:
        E_arr = np.asarray([edep_by_evt[evt] for evt in common_evt_pulse], dtype=float) * GeV_to_MeV
        Q_arr = np.asarray([qtot_by_evt[evt] for evt in common_evt_pulse], dtype=float)

        plt.figure()
        plt.hist(Q_arr, bins=args.bins1d)
        plt.xlabel("TotalPulseCharge")
        plt.ylabel("Counts")
        plt.title("Total pulse charge distribution")
        plt.tight_layout()
        plt.savefig(out_dir / f"{out_prefix}_1dPulseQ_dist.png", dpi=200)
        
        # Plot 3: TotalPulseCharge vs Edep
        USE_ZOOM = True

        for zoom in [False, USE_ZOOM]:
            suffix = "_zoom" if zoom else "_full"
            plt.figure()
            plt.hist2d(E_arr, Q_arr, bins=args.bins)
            if zoom:
                plt.xlim(x_min, x_max)
                title = "Pulse charge vs deposited energy (zoomed)"
            else:
                title = "Pulse charge vs deposited energy"
            plt.xlabel("TotalDepositedEnergyLAr(Truth)[MeV]")
            plt.ylabel("TotalPulseCharge")
            plt.title(title)
            cb = plt.colorbar()
            cb.set_label("Counts")
            plt.tight_layout()
            plt.savefig(out_dir / f"{out_prefix}_PulseQ_vs_truthE{suffix}.png", dpi=200)
            plt.close()
        
    # Plot 5: 1D distribution of EventCharge90NS (summed per event_id)
    common_evt_q90 = sorted(set(edep_by_evt) & set(q90sum_by_evt))
    if len(common_evt_q90) > 0:
        E2 = np.asarray([edep_by_evt[evt] for evt in common_evt_q90], dtype=float) * GeV_to_MeV
        Q2 = np.asarray([q90sum_by_evt[evt] for evt in common_evt_q90], dtype=float)

        plt.figure()
        plt.hist(Q2, bins=args.bins1d)
        plt.xlabel("UncoatedAndCoated8inEventCharge90NS")
        plt.ylabel("Counts")
        plt.title("Event charge (90ns) distribution")
        plt.tight_layout()
        plt.savefig(out_dir / f"{out_prefix}_1deventQ90_dist.png", dpi=200)
        
        # Plot 6: Charge90NS vs Edep
        USE_ZOOM = True
        
        for zoom in [False, USE_ZOOM]:
            suffix = "_zoom" if zoom else "_full"
            plt.figure()
            plt.hist2d(E2, Q2, bins=args.bins)
            if zoom:
                plt.xlim(x_min, x_max)
                title = "Event charge(90ns) vs deposited energy (zoomed)"
            else:
                title = "Event charge(90ns) vs deposited energy"
            plt.xlabel("TotalDepositedEnergyLAr[MeV]")
            plt.ylabel("UncoatedAndCoated8inEventCharge90NS")
            plt.title(title)
            cb = plt.colorbar()
            cb.set_label("Counts")
            plt.tight_layout()
            plt.savefig(out_dir / f"{out_prefix}_eventQ90_vs_truthE{suffix}.png", dpi=200)
            plt.close()
    
    print("Done.")
    if len(common_evt_pulse) == 0:
        print("No events. check pulses exist in DAQ frames and event_id matches Physics frames.")
    if len(common_evt_q90) == 0:
        print("No events. check IntervalChargeSum enabled or key missing.")


if __name__ == "__main__":
    main()
