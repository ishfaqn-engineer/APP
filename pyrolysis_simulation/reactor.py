"""
Reactor model for a fixed bed pyrolysis reactor.

This module solves the coupled ODE system that tracks the mass of
each biomass pseudo-component (cellulose, hemicellulose, lignin) and
the cumulative yields of gas, tar (bio-oil), and char as functions
of time, subject to a prescribed temperature programme.
"""

import numpy as np
from scipy.integrate import solve_ivp

from .kinetics import KINETICS, decomposition_rates

# -----------------------------------------------------------------------
# Default biomass composition (dry, ash-free basis, mass fractions)
# -----------------------------------------------------------------------
DEFAULT_COMPOSITION = {
    "cellulose": 0.45,
    "hemicellulose": 0.25,
    "lignin": 0.25,
    # The remainder (0.05) is treated as inert ash
}

# Typical proximate-analysis moisture & ash for reference
DEFAULT_MOISTURE = 0.05   # 5 wt%
DEFAULT_ASH = 0.05        # 5 wt%

COMPONENTS = list(KINETICS.keys())  # cellulose, hemicellulose, lignin
PRODUCTS = ["gas", "tar", "char"]

# -----------------------------------------------------------------------
# Gas species split factors (mass basis, approximate)
# These describe how the *gas* fraction from each component is split
# into individual gas species.  Values are indicative and can be
# replaced with experimental data.
# -----------------------------------------------------------------------
GAS_SPECIES_SPLIT = {
    "cellulose": {"CO": 0.40, "CO2": 0.35, "CH4": 0.05, "H2": 0.02, "H2O": 0.18},
    "hemicellulose": {"CO": 0.30, "CO2": 0.45, "CH4": 0.04, "H2": 0.01, "H2O": 0.20},
    "lignin": {"CO": 0.25, "CO2": 0.20, "CH4": 0.20, "H2": 0.10, "H2O": 0.25},
}


def _build_ode_rhs(temp_func):
    """
    Build the right-hand side function for the ODE system.

    State vector layout (12 elements):
        [0:3]  m_cell, m_hemi, m_lign   – remaining mass of each component
        [3:6]  gas_cell, gas_hemi, gas_lign
        [6:9]  tar_cell, tar_hemi, tar_lign
        [9:12] char_cell, char_hemi, char_lign
    """

    def rhs(t, y):
        T = temp_func(t)
        dydt = np.zeros_like(y)

        for i, comp in enumerate(COMPONENTS):
            m = y[i]
            if m < 1e-15:
                continue
            rates = decomposition_rates(comp, T)
            k_total = sum(rates.values())
            dm_dt = -k_total * m
            dydt[i] = dm_dt

            # Product accumulation
            for j, prod in enumerate(PRODUCTS):
                dydt[3 + j * 3 + i] = rates[prod] * m

        return dydt

    return rhs


class FixedBedReactor:
    """
    Fixed bed pyrolysis reactor simulation.

    Parameters
    ----------
    biomass_mass : float
        Initial dry biomass charge [kg].
    composition : dict, optional
        Mass fractions of cellulose, hemicellulose, lignin (dry, ash-free).
    moisture : float, optional
        Initial moisture mass fraction (released as H2O during drying).
    ash : float, optional
        Ash mass fraction (inert).
    """

    def __init__(
        self,
        biomass_mass: float = 1.0,
        composition: dict | None = None,
        moisture: float = DEFAULT_MOISTURE,
        ash: float = DEFAULT_ASH,
    ):
        self.biomass_mass = biomass_mass
        self.composition = composition or dict(DEFAULT_COMPOSITION)
        self.moisture = moisture
        self.ash = ash

        # Compute initial component masses (dry, ash-free)
        daf_mass = biomass_mass * (1.0 - moisture - ash)
        self.initial_masses = {
            comp: daf_mass * self.composition.get(comp, 0.0)
            for comp in COMPONENTS
        }

    # ----- temperature programme helpers --------------------------------

    @staticmethod
    def constant_temperature(T: float):
        """Return a callable that gives constant temperature *T* [K]."""
        def _f(_t):
            return T
        return _f

    @staticmethod
    def linear_ramp(T_start: float, heating_rate: float, T_max: float):
        """
        Return a temperature callable: linear ramp then hold.

        Parameters
        ----------
        T_start : float – initial temperature [K]
        heating_rate : float – heating rate [K/s]
        T_max : float – maximum (hold) temperature [K]
        """
        def _f(t):
            T = T_start + heating_rate * t
            return min(T, T_max)
        return _f

    @staticmethod
    def multi_stage_ramp(stages: list):
        """
        Return a temperature callable for a multi-stage programme.

        Parameters
        ----------
        stages : list of tuples
            Each tuple is (heating_rate [K/s], target_T [K], hold_time [s]).
            The reactor starts at the first stage target minus the ramp.
        """
        # pre-compute transition times
        transitions = []
        t_cursor = 0.0
        T_cursor = stages[0][2] if len(stages) > 0 else 300.0
        # Actually let caller set T_start separately; simplify:
        T_cursor = 300.0  # ambient start

        for rate, target, hold in stages:
            if rate > 0:
                dt_ramp = (target - T_cursor) / rate
            else:
                dt_ramp = 0.0
            transitions.append((t_cursor, T_cursor, rate, target, dt_ramp, hold))
            t_cursor += max(dt_ramp, 0) + hold
            T_cursor = target

        def _f(t):
            for (t0, T0, rate, target, dt_ramp, hold) in transitions:
                if t < t0 + dt_ramp:
                    return T0 + rate * (t - t0)
                if t < t0 + dt_ramp + hold:
                    return target
            # After all stages, hold last temperature
            return transitions[-1][3]

        return _f

    # ----- solve ---------------------------------------------------------

    def run(self, t_end: float, temp_func, dt: float = 1.0):
        """
        Run the simulation.

        Parameters
        ----------
        t_end : float
            Total simulation time [s].
        temp_func : callable
            Temperature as a function of time, T(t) [K].
        dt : float
            Output time step [s].

        Returns
        -------
        results : dict
            Dictionary containing time array, temperature profile,
            component masses, product yields, gas species, etc.
        """
        # Initial conditions
        y0 = np.zeros(12)
        for i, comp in enumerate(COMPONENTS):
            y0[i] = self.initial_masses[comp]

        t_eval = np.arange(0, t_end + dt, dt)

        sol = solve_ivp(
            _build_ode_rhs(temp_func),
            [0, t_end],
            y0,
            method="RK45",
            t_eval=t_eval,
            rtol=1e-8,
            atol=1e-10,
            max_step=dt,
        )

        if not sol.success:
            raise RuntimeError(f"ODE solver failed: {sol.message}")

        time = sol.t
        T_profile = np.array([temp_func(t) for t in time])

        # Extract solution
        m_components = {comp: sol.y[i] for i, comp in enumerate(COMPONENTS)}

        product_yields = {}
        for j, prod in enumerate(PRODUCTS):
            total = np.zeros(len(time), dtype=float)
            for i in range(len(COMPONENTS)):
                total += sol.y[3 + j * 3 + i]
            product_yields[prod] = total

        # Add moisture as water in the gas phase (simple flash at ~373 K)
        moisture_mass = self.biomass_mass * self.moisture
        # Simple model: moisture evaporates linearly between 373–393 K
        moisture_released = np.zeros(len(time), dtype=float)
        for idx, T in enumerate(T_profile):
            if T >= 393:
                moisture_released[idx] = moisture_mass
            elif T > 373:
                moisture_released[idx] = moisture_mass * (T - 373) / 20.0
        product_yields["gas"] += moisture_released

        # Ash remains as part of the solid residue
        ash_mass = self.biomass_mass * self.ash
        product_yields["char"] += ash_mass

        # ---- Gas species breakdown ----
        gas_species = {sp: np.zeros(len(time), dtype=float) for sp in ["CO", "CO2", "CH4", "H2", "H2O"]}
        for i, comp in enumerate(COMPONENTS):
            gas_from_comp = sol.y[3 + i]  # gas from this component
            for sp, frac in GAS_SPECIES_SPLIT[comp].items():
                gas_species[sp] += gas_from_comp * frac
        # Add evaporated moisture to H2O
        gas_species["H2O"] += moisture_released

        # ---- Overall mass balance check ----
        total_out = product_yields["gas"] + product_yields["tar"] + product_yields["char"]

        # ---- Conversion ----
        daf_mass = self.biomass_mass * (1.0 - self.moisture - self.ash)
        remaining = sum(m_components[c] for c in COMPONENTS)
        conversion = 1.0 - remaining / daf_mass

        # ---- Yield percentages (on initial biomass basis) ----
        yield_pct = {
            prod: vals / self.biomass_mass * 100.0
            for prod, vals in product_yields.items()
        }

        results = {
            "time": time,                      # [s]
            "temperature": T_profile,          # [K]
            "component_mass": m_components,    # [kg] remaining
            "product_mass": product_yields,    # [kg] cumulative
            "product_yield_pct": yield_pct,    # [wt%]
            "gas_species_mass": gas_species,   # [kg]
            "conversion": conversion,          # [-] array
            "total_output_mass": total_out,    # [kg] mass balance
            "initial_biomass_mass": self.biomass_mass,
            "ash_mass": ash_mass,
            "moisture_mass": moisture_mass,
        }
        return results
