#!/usr/bin/env python3
"""
Main entry point for running the biomass pyrolysis simulation.

Usage
-----
    python run_simulation.py                  # default settings
    python run_simulation.py --mass 0.5       # 500 g biomass
    python run_simulation.py --T-max 800      # hold at 800 °C
    python run_simulation.py --heating-rate 20 # 20 K/min
    python run_simulation.py --no-plot         # skip plot generation
"""

import argparse
import sys

import numpy as np

from pyrolysis_simulation.reactor import FixedBedReactor
from pyrolysis_simulation.visualization import plot_all, print_summary


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Simulate biomass pyrolysis in a fixed bed reactor.",
    )
    p.add_argument(
        "--mass",
        type=float,
        default=1.0,
        help="Initial dry biomass charge [kg] (default: 1.0)",
    )
    p.add_argument(
        "--cellulose",
        type=float,
        default=0.45,
        help="Mass fraction of cellulose (default: 0.45)",
    )
    p.add_argument(
        "--hemicellulose",
        type=float,
        default=0.25,
        help="Mass fraction of hemicellulose (default: 0.25)",
    )
    p.add_argument(
        "--lignin",
        type=float,
        default=0.25,
        help="Mass fraction of lignin (default: 0.25)",
    )
    p.add_argument(
        "--moisture",
        type=float,
        default=0.05,
        help="Moisture mass fraction (default: 0.05)",
    )
    p.add_argument(
        "--ash",
        type=float,
        default=0.05,
        help="Ash mass fraction (default: 0.05)",
    )
    p.add_argument(
        "--T-start",
        type=float,
        default=25.0,
        help="Starting temperature [°C] (default: 25)",
    )
    p.add_argument(
        "--T-max",
        type=float,
        default=600.0,
        help="Maximum (hold) temperature [°C] (default: 600)",
    )
    p.add_argument(
        "--heating-rate",
        type=float,
        default=10.0,
        help="Heating rate [K/min] (default: 10)",
    )
    p.add_argument(
        "--hold-time",
        type=float,
        default=30.0,
        help="Hold time at T_max [min] (default: 30)",
    )
    p.add_argument(
        "--dt",
        type=float,
        default=5.0,
        help="Output time step [s] (default: 5)",
    )
    p.add_argument(
        "--output",
        type=str,
        default="pyrolysis_results.png",
        help="Output plot filename (default: pyrolysis_results.png)",
    )
    p.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plot generation",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # --- Set up reactor ---
    composition = {
        "cellulose": args.cellulose,
        "hemicellulose": args.hemicellulose,
        "lignin": args.lignin,
    }

    reactor = FixedBedReactor(
        biomass_mass=args.mass,
        composition=composition,
        moisture=args.moisture,
        ash=args.ash,
    )

    # --- Temperature programme ---
    T_start_K = args.T_start + 273.15
    T_max_K = args.T_max + 273.15
    heating_rate_Ks = args.heating_rate / 60.0  # K/min -> K/s

    temp_func = FixedBedReactor.linear_ramp(T_start_K, heating_rate_Ks, T_max_K)

    # Total time = ramp time + hold time
    ramp_time = (T_max_K - T_start_K) / heating_rate_Ks if heating_rate_Ks > 0 else 0
    hold_time_s = args.hold_time * 60.0
    t_end = ramp_time + hold_time_s

    print(f"\n  Heating from {args.T_start:.0f} °C to {args.T_max:.0f} °C "
          f"at {args.heating_rate:.1f} K/min")
    print(f"  Ramp duration: {ramp_time/60:.1f} min | Hold: {args.hold_time:.1f} min")
    print(f"  Total simulation time: {t_end/60:.1f} min\n")

    # --- Run simulation ---
    results = reactor.run(t_end=t_end, temp_func=temp_func, dt=args.dt)

    # --- Output ---
    print_summary(results)

    if not args.no_plot:
        fig = plot_all(results, save_path=args.output)
        print(f"\n  Plots saved to: {args.output}")

    return results


if __name__ == "__main__":
    main()
