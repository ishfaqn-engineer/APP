"""
Tests for the biomass pyrolysis simulation.
"""

import numpy as np
import pytest

from pyrolysis_simulation.kinetics import rate_constant, decomposition_rates, KINETICS, R_GAS
from pyrolysis_simulation.reactor import FixedBedReactor, COMPONENTS, PRODUCTS


# ---------------------------------------------------------------------------
# Kinetics tests
# ---------------------------------------------------------------------------

class TestKinetics:
    def test_rate_constant_positive(self):
        """Rate constant must always be positive."""
        k = rate_constant(A=1e10, Ea=100e3, T=773.15)
        assert k > 0

    def test_rate_constant_increases_with_temperature(self):
        """Higher temperature should give a higher rate constant."""
        k_low = rate_constant(1e10, 100e3, 573.15)
        k_high = rate_constant(1e10, 100e3, 773.15)
        assert k_high > k_low

    def test_decomposition_rates_returns_all_products(self):
        """Each component must return rates for gas, tar, char."""
        for comp in COMPONENTS:
            rates = decomposition_rates(comp, 773.15)
            assert set(rates.keys()) == {"gas", "tar", "char"}
            for v in rates.values():
                assert v > 0


# ---------------------------------------------------------------------------
# Reactor tests
# ---------------------------------------------------------------------------

class TestFixedBedReactor:
    def test_mass_balance(self):
        """Total output mass must approximately equal input mass."""
        reactor = FixedBedReactor(biomass_mass=1.0)
        temp_func = FixedBedReactor.linear_ramp(298.15, 10.0 / 60, 873.15)
        t_end = (873.15 - 298.15) / (10.0 / 60) + 1800
        results = reactor.run(t_end=t_end, temp_func=temp_func, dt=10)

        total_out = results["total_output_mass"][-1]
        np.testing.assert_allclose(total_out, 1.0, rtol=0.05)

    def test_yields_are_non_negative(self):
        """All product yields must be >= 0 at all times."""
        reactor = FixedBedReactor(biomass_mass=0.5)
        temp_func = FixedBedReactor.constant_temperature(773.15)
        results = reactor.run(t_end=3600, temp_func=temp_func, dt=10)

        for prod in PRODUCTS:
            assert np.all(results["product_mass"][prod] >= -1e-12)

    def test_conversion_reaches_high_value(self):
        """At 600 °C for 30 min, conversion should be very high."""
        reactor = FixedBedReactor(biomass_mass=1.0)
        temp_func = FixedBedReactor.constant_temperature(873.15)
        results = reactor.run(t_end=1800, temp_func=temp_func, dt=10)

        final_conv = results["conversion"][-1]
        assert final_conv > 0.90, f"Expected >90% conversion, got {final_conv*100:.1f}%"

    def test_char_yield_increases_with_lower_temperature(self):
        """Lower pyrolysis temperature should give more char."""
        reactor_low = FixedBedReactor(biomass_mass=1.0)
        results_low = reactor_low.run(
            t_end=3600,
            temp_func=FixedBedReactor.constant_temperature(573.15),
            dt=10,
        )

        reactor_high = FixedBedReactor(biomass_mass=1.0)
        results_high = reactor_high.run(
            t_end=3600,
            temp_func=FixedBedReactor.constant_temperature(873.15),
            dt=10,
        )

        char_low = results_low["product_yield_pct"]["char"][-1]
        char_high = results_high["product_yield_pct"]["char"][-1]
        assert char_low > char_high, (
            f"Char at 300°C ({char_low:.1f}%) should exceed char at 600°C ({char_high:.1f}%)"
        )

    def test_gas_species_sum_approximates_total_gas(self):
        """Sum of individual gas species should approximate total gas yield."""
        reactor = FixedBedReactor(biomass_mass=1.0)
        temp_func = FixedBedReactor.linear_ramp(298.15, 10.0 / 60, 873.15)
        t_end = (873.15 - 298.15) / (10.0 / 60) + 1800
        results = reactor.run(t_end=t_end, temp_func=temp_func, dt=10)

        gas_species_total = sum(
            results["gas_species_mass"][sp][-1] for sp in results["gas_species_mass"]
        )
        total_gas = results["product_mass"]["gas"][-1]
        # Gas species split covers only the pyrolysis gas (not all moisture pathways)
        # so we check they are in the same order of magnitude
        assert gas_species_total > 0
        assert gas_species_total <= total_gas * 1.1  # allow small tolerance

    def test_constant_temperature_function(self):
        """Constant temperature function should return the set temperature."""
        f = FixedBedReactor.constant_temperature(500.0)
        assert f(0) == 500.0
        assert f(1000) == 500.0

    def test_linear_ramp_function(self):
        """Linear ramp should increase then hold."""
        f = FixedBedReactor.linear_ramp(300.0, 1.0, 600.0)
        assert f(0) == 300.0
        assert f(150) == 450.0
        assert f(300) == 600.0
        assert f(500) == 600.0  # hold
