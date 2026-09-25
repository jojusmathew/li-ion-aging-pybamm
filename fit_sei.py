"""Fit PyBaMM's SEI growth model to measured loss of lithium inventory (LLI), then test it on unseen cells.

Train: the cells cycled in the 25 and 40 degC chambers. Test: the cells in the 10 degC chamber.
Outputs: printed parameters and errors, results/lli_fit.png, results/lli_vs_temperature.png, results/sei_fit.json
"""

import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pybamm
from scipy.optimize import least_squares
from scipy.stats import t as student_t

HERE = Path(__file__).parent
D_KEY = "SEI solvent diffusivity [m2.s-1]"
E_KEY = "SEI growth activation energy [J.mol-1]"
T_KEY = "Ambient temperature [K]"
COLORS = {10: "#2a78d6", 25: "#1baf7a", 40: "#eb6834"}  # chamber temperature -> colour, same in every figure

# SEI growth limited by solvent diffusion through the SEI: rate = D * c_solvent / thickness, times an Arrhenius
# factor exp(E/R * (1/298.15 K - 1/T)). The rate does not depend on the current, so ageing by cycling equals
# storage for the same time at the same temperature; self_check() verifies this.
# This shortcut holds only for current-independent SEI models; potential-dependent ones need the real cycling
# protocol simulated.
model = pybamm.lithium_ion.SPM({"SEI": "solvent-diffusion limited"})
params = pybamm.ParameterValues("OKane2022")  # LG M50 cell, O'Kane et al. 2022
params.update({"Current function [A]": 0, D_KEY: "[input]", E_KEY: "[input]", T_KEY: "[input]"})
solver = pybamm.IDAKLUSolver(rtol=1e-8, atol=1e-10)  # tight, so finite-difference Jacobians are clean
storage = pybamm.Simulation(model, parameter_values=params, solver=solver)


def load_cells():
    cells = []
    for path in sorted((HERE / "data").glob("*.csv")):
        name, chamber = re.search(r"cell (\w) \((\d+)degC\)", path.name).groups()
        rows = list(csv.DictReader(open(path, encoding="utf-8")))[1:]  # skip the begin-of-life row (LLI = 0)
        cells.append(
            {
                "name": name,
                "chamber": int(chamber),
                # measured cell surface temperature: self-heating puts it several degC above the chamber
                "T": np.mean([float(r["Age set av. temperature [°C]"]) for r in rows]),
                "days": np.array([float(r["Days of degradation"]) for r in rows]),
                "lli": np.array([100 * float(r["LLI"]) for r in rows]),
                "cycles": int(rows[-1]["Ageing Cycles"]),
            }
        )
    return cells


def inputs(x, T_celsius):
    """x = [log10 of D at 25 degC in m2/s, activation energy in kJ/mol]"""
    return {D_KEY: 10 ** x[0], E_KEY: 1e3 * x[1], T_KEY: 273.15 + T_celsius}


def simulate_lli(x, T_celsius, days):
    sol = storage.solve([0, 86400 * max(days)], inputs=inputs(x, T_celsius), initial_soc=0.775)
    return np.atleast_1d(sol["Loss of lithium inventory [%]"](86400 * np.asarray(days)))


def residuals(x, cells):
    return np.concatenate([simulate_lli(x, c["T"], c["days"]) - c["lli"] for c in cells])


def rmse(x, cells):
    return np.sqrt(np.mean(residuals(x, cells) ** 2))


def self_check(x):
    """20 real cycles in the 70-85 % window must lose the same lithium as storage for the same time."""
    cycling = pybamm.Experiment([("Discharge at 1C for 9 minutes", "Charge at 0.3C for 30 minutes")] * 20)
    fixed = params.copy()
    fixed.update(inputs(x, 25))  # experiments cannot take temperature as an input parameter
    sol = pybamm.Simulation(model, parameter_values=fixed, experiment=cycling, solver=solver).solve(initial_soc=0.85)
    lli_cycling = sol["Loss of lithium inventory [%]"].entries[-1]
    lli_storage = simulate_lli(x, 25, [sol["Time [s]"].entries[-1] / 86400])[0]
    assert abs(lli_cycling / lli_storage - 1) < 0.01, (lli_cycling, lli_storage)
    print(f"self-check passed: 20 cycles LLI {lli_cycling:.4f} % vs storage {lli_storage:.4f} %")


def plot_fit(x, cells):
    fig, ax = plt.subplots(figsize=(6.4, 4))
    days = np.linspace(0, 260, 131)
    for chamber in (40, 25, 10):
        group = [c for c in cells if c["chamber"] == chamber]
        T = np.mean([c["T"] for c in group])
        tested = chamber == 10
        for c in group:
            ax.plot(c["days"], c["lli"], "o", color=COLORS[chamber], mec="white", mew=1)
        ax.plot(
            days,
            simulate_lli(x, T, days),
            "--" if tested else "-",
            color=COLORS[chamber],
            label=f"{chamber} °C chamber (cell ≈ {T:.0f} °C): " + ("predicted" if tested else "fitted"),
        )
    ax.plot([], [], "o", color="#898781", label="measured (dots), 2 cells per chamber")
    ax.set(xlabel="Time [days]", ylabel="Loss of lithium inventory [%]", xlim=(0, 265), ylim=(0, None))
    ax.set_title("SEI model fitted to the 25/40 °C cells, then tested on the 10 °C cells")
    ax.legend(loc="lower right")
    fig.savefig(HERE / "results" / "lli_fit.png")


def plot_temperature(cells, T_grid, lli_grid):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(T_grid, lli_grid, color="#52514e", label="SEI-only model (Arrhenius), fitted at 25/40 °C")
    for chamber in (10, 25, 40):
        group = [c for c in cells if c["chamber"] == chamber]
        ax.plot([c["T"] for c in group], [c["lli"][-1] for c in group], "o", ms=8, color=COLORS[chamber],
                mec="white", mew=1, label=f"measured, {chamber} °C chamber")
    cold = [c for c in cells if c["chamber"] == 10]
    T_cold = np.mean([c["T"] for c in cold])
    measured, predicted = np.mean([c["lli"][-1] for c in cold]), np.interp(T_cold, T_grid, lli_grid)
    x_arrow = T_cold - 1.5  # beside the two cold cells rather than on top of them
    ax.annotate("", xy=(x_arrow, measured), xytext=(x_arrow, np.interp(x_arrow, T_grid, lli_grid)),
                arrowprops={"arrowstyle": "->", "color": "#0b0b0b", "lw": 1.2})
    ax.text(x_arrow - 0.8, (measured + predicted) / 2, f"+{measured - predicted:.1f} %-points\nnot explained by SEI",
            ha="right", va="center", color="#0b0b0b")
    ax.set(xlabel="Cell temperature during cycling [°C]", ylabel="LLI after 258 days [%]")
    ax.set_title("Cold cells age faster than the SEI model predicts")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    fig.savefig(HERE / "results" / "lli_vs_temperature.png")


if __name__ == "__main__":
    cells = load_cells()
    train = [c for c in cells if c["chamber"] != 10]
    test = [c for c in cells if c["chamber"] == 10]

    fit = least_squares(residuals, x0=[-19.5, 38], args=(train,), diff_step=1e-4, bounds=([-24, 0], [-16, 150]))
    x = fit.x
    # 95 % confidence intervals from the Jacobian; they assume independent errors, which is optimistic because
    # the check-ups of one cell are correlated.
    dof = fit.fun.size - x.size
    cov = np.linalg.inv(fit.jac.T @ fit.jac) * 2 * fit.cost / dof
    ci = student_t.ppf(0.975, dof) * np.sqrt(np.diag(cov))
    print(f"SEI solvent diffusivity at 25 degC: {10 ** x[0]:.2e} m2/s (log10 = {x[0]:.2f} ± {ci[0]:.2f})")
    print(f"SEI growth activation energy: {x[1]:.1f} ± {ci[1]:.1f} kJ/mol")
    print(f"RMSE train (25/40 degC cells): {rmse(x, train):.2f} %-points")
    print(f"RMSE test (10 degC cells):     {rmse(x, test):.2f} %-points")
    for c in cells:
        print(f"  cell {c['name']} ({c['chamber']} degC chamber, {c['T']:.1f} degC measured): "
              f"LLI after {c['days'][-1]:.0f} d measured {c['lli'][-1]:.1f} %, model {simulate_lli(x, c['T'], c['days'])[-1]:.1f} %")
    self_check(x)

    T_grid = np.arange(0, 46, 1.0)
    lli_grid = np.array([simulate_lli(x, T, [258])[0] for T in T_grid])
    T_cold = np.mean([c["T"] for c in test])
    gap = np.mean([c["lli"][-1] for c in test]) - np.interp(T_cold, T_grid, lli_grid)
    print(f"10 degC cells: {gap:.1f} %-points more LLI than the SEI model predicts")

    plot_fit(x, cells)
    plot_temperature(cells, T_grid, lli_grid)
    chamber_T = {ch: np.mean([c["T"] for c in cells if c["chamber"] == ch]) for ch in (10, 25, 40)}
    json.dump({"log10_D": x[0], "Ea_kJ_mol": x[1], "chamber_T": chamber_T, "cold_gap_pct": gap,
               "cycles": test[0]["cycles"]},
              open(HERE / "results" / "sei_fit.json", "w"), indent=2)
