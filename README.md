# Biomass Pyrolysis in Fixed Bed Reactor – Simulation

A Python-based simulation of **biomass pyrolysis in a fixed bed reactor**, modelling the thermal decomposition of biomass into **gas, bio-oil (tar), and char** as functions of **time and temperature**.

## Features

| Feature | Description |
|---|---|
| **Multi-component kinetics** | Cellulose, hemicellulose, and lignin decompose via parallel first-order Arrhenius reactions |
| **Product yields** | Cumulative mass of gas, bio-oil (tar), and char tracked over time |
| **Gas species** | Individual gas species: CO, CO₂, CH₄, H₂, H₂O |
| **Temperature programmes** | Constant temperature, linear ramp + hold, or custom multi-stage profiles |
| **Moisture & ash** | Moisture evaporation and inert ash accounted for |
| **Mass balance** | Verified mass balance across all products |
| **Visualisation** | Dashboard with 6 plots: temperature profile, product yields, component decomposition, gas species, conversion, and final product pie chart |

## Quick Start

### Install dependencies

```bash
pip install -r requirements.txt
```

### Run the simulation (default settings)

```bash
python run_simulation.py
```

This runs a 1 kg biomass charge heated from 25 °C to 600 °C at 10 K/min, then held for 30 min. Results are printed to the console and plots saved to `pyrolysis_results.png`.

### Customise parameters

```bash
# 500 g biomass, heated to 800 °C at 20 K/min
python run_simulation.py --mass 0.5 --T-max 800 --heating-rate 20

# Change biomass composition
python run_simulation.py --cellulose 0.50 --hemicellulose 0.20 --lignin 0.25

# Skip plot generation (console output only)
python run_simulation.py --no-plot
```

### All CLI options

| Option | Default | Description |
|---|---|---|
| `--mass` | 1.0 | Initial dry biomass [kg] |
| `--cellulose` | 0.45 | Cellulose mass fraction |
| `--hemicellulose` | 0.25 | Hemicellulose mass fraction |
| `--lignin` | 0.25 | Lignin mass fraction |
| `--moisture` | 0.05 | Moisture mass fraction |
| `--ash` | 0.05 | Ash mass fraction |
| `--T-start` | 25 | Starting temperature [°C] |
| `--T-max` | 600 | Maximum temperature [°C] |
| `--heating-rate` | 10 | Heating rate [K/min] |
| `--hold-time` | 30 | Hold time at T_max [min] |
| `--dt` | 5 | Output time step [s] |
| `--output` | pyrolysis_results.png | Plot output file |
| `--no-plot` | — | Skip plot generation |

## Project Structure

```
├── run_simulation.py              # Main entry point (CLI)
├── requirements.txt               # Python dependencies
├── pyrolysis_simulation/
│   ├── __init__.py
│   ├── kinetics.py                # Arrhenius kinetic model
│   ├── reactor.py                 # Fixed bed reactor ODE model
│   └── visualization.py           # Plotting and summary output
└── tests/
    └── test_simulation.py         # Unit tests
```

## Model Description

The simulation uses a **competitive parallel reaction scheme**:

```
                    ┌─→ Gas   (CO, CO₂, CH₄, H₂, H₂O)
Cellulose      ────┼─→ Tar   (bio-oil)
                    └─→ Char

                    ┌─→ Gas
Hemicellulose  ────┼─→ Tar
                    └─→ Char

                    ┌─→ Gas
Lignin         ────┼─→ Tar
                    └─→ Char
```

Each reaction follows first-order Arrhenius kinetics: **k = A · exp(−Eₐ / RT)**.

Kinetic parameters are based on published literature (Ranzi et al.; Di Blasi, *Prog. Energy Combust. Sci.*, 2008).

## Run Tests

```bash
pip install pytest
pytest tests/ -v
```

## Example Output

```
  BIOMASS PYROLYSIS IN FIXED BED REACTOR – SIMULATION SUMMARY
======================================================================
  Initial biomass charge    : 1000.0 g
  Simulation time           : 5250 s
  Final temperature         : 600.0 °C
  Final conversion (DAF)    : 98.7 %
----------------------------------------------------------------------
  FINAL PRODUCT YIELDS (wt% of initial biomass):
----------------------------------------------------------------------
    Gas               :  35.42 %  (  354.20 g)
    Bio-oil (tar)     :  38.15 %  (  381.50 g)
    Char              :  21.33 %  (  213.30 g)
----------------------------------------------------------------------
```
