"""Is lithium plating behind the extra ageing of the cold cells? And where does fast charging start to plate?

Criterion: lithium metal can plate once the graphite surface potential vs Li/Li+ at the separator side, where it is
lowest during charging, drops below 0 V. Run fit_sei.py first; this script reuses its fitted SEI parameters.
Outputs: printed hypothesis check, results/fast_charge_limits.csv, results/fast_charge_map.png
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pybamm
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).parent
PHI = "Negative electrode surface potential difference at separator interface [V]"
T_KEY = "Ambient temperature [K]"
fit = json.load(open(HERE / "results" / "sei_fit.json"))


def plating_check(T_celsius, cycles=20):
    """Cycle between 70 and 85 % SoC as in the experiment. Returns the lowest anode potential [V] and the
    lithium lost to plating over the whole 258-day test [% of the lithium inventory]."""
    model = pybamm.lithium_ion.DFN({"SEI": "solvent-diffusion limited", "lithium plating": "partially reversible"})
    params = pybamm.ParameterValues("OKane2022")
    params.update(
        {
            "SEI solvent diffusivity [m2.s-1]": 10 ** fit["log10_D"],
            "SEI growth activation energy [J.mol-1]": 1e3 * fit["Ea_kJ_mol"],
            T_KEY: 273.15 + T_celsius,
            "Initial temperature [K]": 273.15 + T_celsius,
        }
    )
    experiment = pybamm.Experiment([("Charge at 0.3C for 30 minutes", "Discharge at 1C for 9 minutes")] * cycles)
    sol = pybamm.Simulation(model, parameter_values=params, experiment=experiment).solve(initial_soc=0.70)
    plated = [c["Loss of lithium to negative lithium plating [mol]"].entries[-1] for c in sol.cycles]
    # Assumes the settled per-cycle loss (second half of the run) stays constant over all cycles of the test;
    # simulating every cycle would be needed if plating turned out to matter.
    half = cycles // 2
    per_cycle = (plated[-1] - plated[half - 1]) / (cycles - half)
    lost = 100 * per_cycle * fit["cycles"] / sol["Total lithium in particles [mol]"].entries[0]
    return sol[PHI].entries.min(), lost


def fast_charge_limits(temps, c_rates):
    """SoC [%] at which a constant-current charge from 10 % SoC first drives the anode below 0 V vs Li.
    100 means the 4.2 V limit is reached first, i.e. no plating in the constant-current phase."""
    params = pybamm.ParameterValues("OKane2022")
    params.update({"Current function [A]": "[input]", T_KEY: "[input]"})
    sim = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=params)  # fresh cell, isothermal
    capacity = params["Nominal cell capacity [A.h]"]
    limits = np.empty((len(temps), len(c_rates)))
    for i, T in enumerate(temps):
        for j, c in enumerate(c_rates):
            sol = sim.solve([0, 3600 / c], inputs={"Current function [A]": -c * capacity, T_KEY: 273.15 + T},
                            initial_soc=0.1)
            phi = sol[PHI].entries
            soc = 10 + 100 * c * sol["Time [s]"].entries / 3600
            below = np.flatnonzero(phi < 0)
            if below.size == 0:
                limits[i, j] = 100
            elif below[0] == 0:
                limits[i, j] = soc[0]
            else:  # interpolate the 0 V crossing between the last safe and the first plating point
                k = below[0]
                limits[i, j] = np.interp(0, [phi[k], phi[k - 1]], [soc[k], soc[k - 1]])
    return limits


def plot_map(temps, c_rates, limits):
    ramp = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#0d366b"])
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.imshow(limits, cmap=ramp, vmin=0, vmax=100, origin="lower", aspect="auto")
    for i in range(len(temps)):
        for j in range(len(c_rates)):
            value = limits[i, j]
            ax.text(j, i, "no plating" if value == 100 else f"{value:.0f} %", ha="center", va="center",
                    color="white" if value > 55 else "#0b0b0b")
    # LG M50T datasheet, max continuous charge current: 0.3C at 0-25 degC, 0.7C at 25-45 degC
    x_cold, x_warm, y_25 = c_rates.index(0.3) + 0.5, c_rates.index(0.7) + 0.5, temps.index(25) - 0.5
    ax.plot([x_cold, x_cold, x_warm, x_warm], [-0.5, y_25, y_25, len(temps) - 0.5], color="#eb6834", lw=3,
            label="max. charge current allowed by the LG datasheet")
    ax.set_xticks(range(len(c_rates)), [f"{c:g}C" for c in c_rates])
    ax.set_yticks(range(len(temps)), [f"{T} °C" for T in temps])
    ax.set_xticks(np.arange(len(c_rates) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(temps) + 1) - 0.5, minor=True)
    ax.grid(False)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.set(xlabel="Charging current (C-rate)", ylabel="Cell temperature")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2))
    ax.set_title("Plating-free fast charging: max SoC before the anode reaches 0 V vs Li")
    fig.savefig(HERE / "results" / "fast_charge_map.png")


if __name__ == "__main__":
    print("Hypothesis: lithium plating causes the extra lithium loss of the cold cells")
    checks = {chamber: plating_check(T) for chamber, T in fit["chamber_T"].items()}
    for chamber, (phi_min, lost) in checks.items():
        print(f"  {chamber} degC chamber ({fit['chamber_T'][chamber]:.1f} degC cell): lowest anode potential "
              f"{1e3 * phi_min:.0f} mV vs Li, lithium lost to plating over the test {lost:.2f} %")
    phi_cold, lost_cold = checks["10"]
    gap = fit["cold_gap_pct"]
    verdict = "can" if phi_cold < 0 or lost_cold > gap / 2 else "cannot"
    print(f"  gap to explain: {gap:.1f} %-points -> plating {verdict} explain it with these parameters")

    # 0.3C and 0.7C are the datasheet's max charge currents below and above 25 degC
    temps, c_rates = [0, 10, 25, 40], [0.3, 0.5, 0.7, 1, 1.5, 2]
    limits = fast_charge_limits(temps, c_rates)
    # physics sanity check: plating starts earlier at higher current and at lower temperature
    assert np.all(np.diff(limits, axis=1) <= 1e-9) and np.all(np.diff(limits, axis=0) >= -1e-9), limits
    with open(HERE / "results" / "fast_charge_limits.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["T_degC"] + [f"{c:g}C" for c in c_rates])
        writer.writerows([T] + [round(v, 1) for v in row] for T, row in zip(temps, limits))
    print("Plating-free SoC limit [%] (rows: T in degC, columns: C-rate):")
    print(open(HERE / "results" / "fast_charge_limits.csv").read())
    plot_map(temps, c_rates, limits)
