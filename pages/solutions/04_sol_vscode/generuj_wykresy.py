"""Generate projectile-motion plots for raport.md.

Run from the repository root with:
    python pages/solutions/04_sol_vscode/generuj_wykresy.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

OUT = Path(__file__).parent / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

G = 9.81
RHO = 1.225
CD = 0.47
DIAMETER = 0.073
MASS = 0.145
AREA = np.pi * (DIAMETER / 2) ** 2
BETA = RHO * CD * AREA / (2 * MASS)
V0 = 20.0
ANGLE_DEG = 45.0
BLUE = "#2463A6"
ORANGE = "#D47722"
GREEN = "#168578"
INK = "#263238"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def vacuum_trajectory(speed, angle_deg, times):
    angle = np.deg2rad(angle_deg)
    vx0 = speed * np.cos(angle)
    vy0 = speed * np.sin(angle)
    x = vx0 * times
    y = vy0 * times - 0.5 * G * times**2
    vx = np.full_like(times, vx0)
    vy = vy0 - G * times
    return x, y, vx, vy


def drag_rhs(_time, state):
    _, _, vx, vy = state
    speed = np.hypot(vx, vy)
    return [vx, vy, -BETA * speed * vx, -G - BETA * speed * vy]


def ground_event(_time, state):
    return state[1]


ground_event.terminal = True
ground_event.direction = -1


def drag_trajectory(speed, angle_deg, sample_count=1200):
    angle = np.deg2rad(angle_deg)
    initial = [0.0, 0.0, speed * np.cos(angle), speed * np.sin(angle)]
    solution = solve_ivp(
        drag_rhs,
        (0.0, 30.0),
        initial,
        events=ground_event,
        dense_output=True,
        rtol=1e-9,
        atol=1e-11,
        max_step=0.02,
    )
    if not solution.success or len(solution.t_events[0]) == 0:
        raise RuntimeError("Nie udało się wyznaczyć chwili lądowania.")
    landing_time = solution.t_events[0][0]
    times = np.linspace(0.0, landing_time, sample_count)
    return times, solution.sol(times), landing_time


def save(fig, filename):
    fig.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


angle = np.deg2rad(ANGLE_DEG)
vacuum_flight_time = 2 * V0 * np.sin(angle) / G
t_vacuum = np.linspace(0.0, vacuum_flight_time, 1200)
x_vacuum, y_vacuum, vx_vacuum, vy_vacuum = vacuum_trajectory(
    V0, ANGLE_DEG, t_vacuum
)
t_drag, state_drag, drag_flight_time = drag_trajectory(V0, ANGLE_DEG)
x_drag, y_drag, vx_drag, vy_drag = state_drag
range_vacuum = x_vacuum[-1]
range_drag = x_drag[-1]

fig, ax = plt.subplots(figsize=(8.5, 5.2), layout="constrained")
ax.plot(x_vacuum, y_vacuum, color=BLUE, label="Bez oporu")
ax.plot(x_drag, y_drag, color=ORANGE, label="Opór kwadratowy")
ax.scatter([0, range_vacuum, range_drag], [0, 0, 0], color=INK, s=24, zorder=3)
ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]",
       title="Tor lotu: model idealny i opór powietrza")
ax.set_ylim(bottom=0)
ax.legend()
save(fig, "01_tory_lotu.png")

fig, axes = plt.subplots(2, 1, figsize=(8.5, 6), sharex=True, layout="constrained")
axes[0].plot(t_vacuum, x_vacuum, color=BLUE, label="Bez oporu")
axes[0].plot(t_drag, x_drag, color=ORANGE, label="Opór kwadratowy")
axes[0].set_ylabel("x [m]")
axes[0].set_title("Współrzędne w funkcji czasu")
axes[0].legend()
axes[1].plot(t_vacuum, y_vacuum, color=BLUE, label="Bez oporu")
axes[1].plot(t_drag, y_drag, color=ORANGE, label="Opór kwadratowy")
axes[1].set(xlabel="Czas t [s]", ylabel="y [m]")
axes[1].axhline(0, color=INK, lw=0.8)
save(fig, "02_polozenie_czas.png")

fig, axes = plt.subplots(2, 1, figsize=(8.5, 6), sharex=True, layout="constrained")
axes[0].plot(t_vacuum, vx_vacuum, color=BLUE, label="Bez oporu")
axes[0].plot(t_drag, vx_drag, color=ORANGE, label="Opór kwadratowy")
axes[0].set_ylabel("v_x [m/s]")
axes[0].set_title("Składowe prędkości")
axes[0].legend()
axes[1].plot(t_vacuum, vy_vacuum, color=BLUE, label="Bez oporu")
axes[1].plot(t_drag, vy_drag, color=ORANGE, label="Opór kwadratowy")
axes[1].axhline(0, color=INK, lw=0.8)
axes[1].set(xlabel="Czas t [s]", ylabel="v_y [m/s]")
save(fig, "03_predkosc_czas.png")

speed_drag = np.hypot(vx_drag, vy_drag)
ax_drag = -BETA * speed_drag * vx_drag
ay_drag = -G - BETA * speed_drag * vy_drag
fig, axes = plt.subplots(2, 1, figsize=(8.5, 6), sharex=True, layout="constrained")
axes[0].plot(t_vacuum, np.zeros_like(t_vacuum), color=BLUE, label="Bez oporu: a_x")
axes[0].plot(t_drag, ax_drag, color=ORANGE, label="Opór: a_x")
axes[0].set_ylabel("a_x [m/s²]")
axes[0].set_title("Przyspieszenie w modelu z oporem")
axes[0].legend()
axes[1].plot(t_vacuum, -G * np.ones_like(t_vacuum), color=BLUE, label="Bez oporu")
axes[1].plot(t_drag, ay_drag, color=ORANGE, label="Opór kwadratowy")
axes[1].set(xlabel="Czas t [s]", ylabel="a_y [m/s²]")
axes[1].legend()
save(fig, "04_przyspieszenie_czas.png")

energy_vacuum = 0.5 * MASS * (vx_vacuum**2 + vy_vacuum**2) + MASS * G * y_vacuum
energy_drag = 0.5 * MASS * (vx_drag**2 + vy_drag**2) + MASS * G * y_drag
fig, ax = plt.subplots(figsize=(8.5, 4.6), layout="constrained")
ax.plot(t_vacuum, energy_vacuum, color=BLUE, label="Bez oporu")
ax.plot(t_drag, energy_drag, color=ORANGE, label="Opór kwadratowy")
ax.set(xlabel="Czas t [s]", ylabel="E_m [J]",
       title="Energia mechaniczna układu piłka–Ziemia")
ax.legend()
save(fig, "05_energia.png")

angles_deg = np.linspace(10.0, 80.0, 36)
angles = np.deg2rad(angles_deg)
ranges_vacuum = V0**2 * np.sin(2 * angles) / G
ranges_drag = np.array([
    drag_trajectory(V0, angle_value, sample_count=300)[1][0, -1]
    for angle_value in angles_deg
])
fig, ax = plt.subplots(figsize=(8.5, 4.8), layout="constrained")
ax.plot(angles, ranges_vacuum, color=BLUE, label="Bez oporu")
ax.plot(angles, ranges_drag, color=ORANGE, label="Opór kwadratowy")
ax.axvline(np.deg2rad(ANGLE_DEG), color=INK, ls=":", lw=1.2,
           label=f"Kąt bazowy: {np.deg2rad(ANGLE_DEG):.3f} rad")
ax.set(xlabel="Kąt wyrzutu α [rad]", ylabel="Zasięg R [m]",
       title="Zasięg a kąt wyrzutu (v₀ = 20 m/s)")
ax.legend()
save(fig, "06_zasieg_kat.png")

print(f"v0 = {V0:.2f} m/s; alpha = {ANGLE_DEG:.2f} deg")
print(f"Model próżniowy: T = {vacuum_flight_time:.4f} s, R = {range_vacuum:.4f} m, H = {np.max(y_vacuum):.4f} m")
print(f"Opór kwadratowy: T = {drag_flight_time:.4f} s, R = {range_drag:.4f} m, H = {np.max(y_drag):.4f} m")
print(f"Zapisano 6 wykresów w: {OUT}")
