"""
Visualisation utilities for pyrolysis simulation results.
"""

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI / headless environments
import matplotlib.pyplot as plt
import numpy as np


def plot_temperature_profile(results, ax=None, celsius=True):
    """Plot reactor temperature vs time."""
    if ax is None:
        _, ax = plt.subplots()
    t = results["time"]
    T = results["temperature"]
    if celsius:
        T = T - 273.15
    ax.plot(t, T, "r-", linewidth=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]" if celsius else "Temperature [K]")
    ax.set_title("Reactor Temperature Profile")
    ax.grid(True, alpha=0.3)
    return ax


def plot_product_yields(results, ax=None):
    """Plot cumulative product yields (wt%) vs time."""
    if ax is None:
        _, ax = plt.subplots()
    t = results["time"]
    for prod, vals in results["product_yield_pct"].items():
        ax.plot(t, vals, linewidth=2, label=prod.capitalize())
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Yield [wt%]")
    ax.set_title("Product Yields vs Time")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


def plot_component_mass(results, ax=None):
    """Plot remaining biomass component mass vs time."""
    if ax is None:
        _, ax = plt.subplots()
    t = results["time"]
    for comp, mass in results["component_mass"].items():
        ax.plot(t, mass * 1000, linewidth=2, label=comp.capitalize())
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Remaining Mass [g]")
    ax.set_title("Biomass Component Decomposition")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


def plot_gas_species(results, ax=None):
    """Plot cumulative gas species mass vs time."""
    if ax is None:
        _, ax = plt.subplots()
    t = results["time"]
    for sp, mass in results["gas_species_mass"].items():
        ax.plot(t, mass * 1000, linewidth=2, label=sp)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Cumulative Mass [g]")
    ax.set_title("Gas Species Evolution")
    ax.legend()
    ax.grid(True, alpha=0.3)
    return ax


def plot_conversion(results, ax=None):
    """Plot biomass conversion vs time."""
    if ax is None:
        _, ax = plt.subplots()
    t = results["time"]
    ax.plot(t, results["conversion"] * 100, "k-", linewidth=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Conversion [%]")
    ax.set_title("Overall Biomass Conversion")
    ax.grid(True, alpha=0.3)
    return ax


def plot_final_product_pie(results, ax=None):
    """Pie chart of final product distribution."""
    if ax is None:
        _, ax = plt.subplots()
    labels = []
    sizes = []
    for prod in ["gas", "tar", "char"]:
        labels.append(prod.capitalize() if prod != "tar" else "Bio-oil (tar)")
        sizes.append(results["product_yield_pct"][prod][-1])
    colors = ["#ff9999", "#66b3ff", "#99ff99"]
    ax.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=140)
    ax.set_title("Final Product Distribution [wt%]")
    return ax


def plot_all(results, save_path=None):
    """Generate a comprehensive dashboard of all plots."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        "Biomass Pyrolysis in Fixed Bed Reactor – Simulation Results",
        fontsize=14,
        fontweight="bold",
    )

    plot_temperature_profile(results, ax=axes[0, 0])
    plot_product_yields(results, ax=axes[0, 1])
    plot_component_mass(results, ax=axes[0, 2])
    plot_gas_species(results, ax=axes[1, 0])
    plot_conversion(results, ax=axes[1, 1])
    plot_final_product_pie(results, ax=axes[1, 2])

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def print_summary(results):
    """Print a tabular summary of the simulation results."""
    t = results["time"]
    T = results["temperature"]
    m0 = results["initial_biomass_mass"]

    print("=" * 70)
    print("  BIOMASS PYROLYSIS IN FIXED BED REACTOR – SIMULATION SUMMARY")
    print("=" * 70)
    print(f"  Initial biomass charge    : {m0*1000:.1f} g")
    print(f"  Ash content               : {results['ash_mass']*1000:.1f} g")
    print(f"  Moisture content          : {results['moisture_mass']*1000:.1f} g")
    print(f"  Simulation time           : {t[-1]:.0f} s")
    print(f"  Final temperature         : {T[-1]-273.15:.1f} °C")
    print(f"  Final conversion (DAF)    : {results['conversion'][-1]*100:.1f} %")
    print("-" * 70)
    print("  FINAL PRODUCT YIELDS (wt% of initial biomass):")
    print("-" * 70)
    for prod in ["gas", "tar", "char"]:
        label = "Bio-oil (tar)" if prod == "tar" else prod.capitalize()
        pct = results["product_yield_pct"][prod][-1]
        mass_g = results["product_mass"][prod][-1] * 1000
        print(f"    {label:<18s}: {pct:6.2f} %  ({mass_g:8.2f} g)")
    print("-" * 70)
    total_yield = sum(results["product_yield_pct"][p][-1] for p in ["gas", "tar", "char"])
    print(f"    {'Total':<18s}: {total_yield:6.2f} %")
    print("-" * 70)
    print("  GAS SPECIES COMPOSITION (final, mass basis):")
    print("-" * 70)
    total_gas = results["product_mass"]["gas"][-1]
    for sp, mass_arr in results["gas_species_mass"].items():
        sp_mass = mass_arr[-1]
        sp_pct = sp_mass / total_gas * 100 if total_gas > 0 else 0
        print(f"    {sp:<6s}: {sp_mass*1000:8.2f} g  ({sp_pct:5.1f} % of gas)")
    print("-" * 70)
    print(f"  Mass balance check: {results['total_output_mass'][-1]*1000:.2f} g "
          f"out of {m0*1000:.1f} g input")
    print("=" * 70)
