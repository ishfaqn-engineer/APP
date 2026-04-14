"""
Kinetic model for biomass pyrolysis in a fixed bed reactor.

This module implements a multi-component competitive kinetic scheme
for the primary pyrolysis of biomass. Biomass is modelled as three
pseudo-components (cellulose, hemicellulose, lignin), each
decomposing independently via parallel first-order Arrhenius
reactions into gas, tar (bio-oil), and char.

Reference kinetic parameters are based on the widely used scheme by
Ranzi et al. and Di Blasi (Prog. Energy Combust. Sci., 2008).
"""

import numpy as np


# ---------------------------------------------------------------------------
# Kinetic parameters: (A [1/s], Ea [J/mol]) for each component -> product
# ---------------------------------------------------------------------------
# Each pseudo-component produces gas, tar (bio-oil), and char.
KINETICS = {
    "cellulose": {
        "gas":  {"A": 2.80e19, "Ea": 242.4e3},
        "tar":  {"A": 3.28e14, "Ea": 196.5e3},
        "char": {"A": 1.30e10, "Ea": 150.5e3},
    },
    "hemicellulose": {
        "gas":  {"A": 9.20e16, "Ea": 186.7e3},
        "tar":  {"A": 2.10e16, "Ea": 186.7e3},
        "char": {"A": 3.05e1,  "Ea": 46.0e3},
    },
    "lignin": {
        "gas":  {"A": 7.70e6,  "Ea": 111.4e3},
        "tar":  {"A": 1.50e9,  "Ea": 143.8e3},
        "char": {"A": 7.70e6,  "Ea": 111.4e3},
    },
}

# Gas-phase universal constant
R_GAS = 8.314  # J/(mol·K)


def rate_constant(A: float, Ea: float, T: float) -> float:
    """Arrhenius rate constant k = A * exp(-Ea / (R*T))."""
    return A * np.exp(-Ea / (R_GAS * T))


def decomposition_rates(component: str, T: float) -> dict:
    """Return rate constants for gas, tar, and char for *component* at *T* [K]."""
    params = KINETICS[component]
    return {
        product: rate_constant(p["A"], p["Ea"], T)
        for product, p in params.items()
    }
