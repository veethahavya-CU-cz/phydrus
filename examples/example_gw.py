"""Floating water-table example for HYDRUS-1D via phydrus.

Built in three stages to isolate numerical issues:

    Stage 1:  constant-flux recharge, zero bottom flux, no roots.
                        Verifies basic saturated-unsaturated redistribution.
    Stage 2:  time-varying atmospheric BC (Prec + rSoil).
    Stage 3:  shallow root uptake confined to the vadose zone.

WHY roots in the saturated zone stall the solver
------------------------------------------------
The Feddes stress function has an oxygen-stress cutoff at p0 (default -10 cm
h). Any root node where h > p0 switches off uptake. In our closed-bottom
system the saturated zone is thick, so many root nodes straddle this
sharp discontinuity at every time step. The modified Picard iteration
must resolve this steeply non-linear sink, which drives dt -> dtMin
and the solver hangs for days of wall-clock time.

Root-depth rules (for Stage 3):
 - max root depth must stay >= 50 cm ABOVE (shallower than) the capillary
   fringe top, i.e.  xrmax  <=  |z_wt| - 1/alpha  - 50  [cm]
 - with alpha=0.036 (loam): fringe top ≈ 28 cm above WT
   => for z_wt=-150: safe limit ≈ 150 - 28 - 50 = 72 cm depth
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for terminal runs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import phydrus as ps

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def resolve_hydrus_exe() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    for candidate in [
        repo_root / "hydrus1d" / "bin" / "hydrus",
        Path("/opt/hydrus1d/bin/hydrus"),
    ]:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("HYDRUS-1D executable not found.")


def hydrostatic_h(z: np.ndarray, z_wt: float) -> np.ndarray:
    """Pressure head h [cm] for hydrostatic equilibrium: h(z) = z_wt - z."""
    return z_wt - z


def interpolate_gwl_from_nod_inf(nod: dict) -> pd.DataFrame:
    """Interpolate WT elevation (h=0 crossing) from NOD_INF snapshots.

    Uses read_nod_inf() output (dict of time -> DataFrame with Depth and Head
    columns) so it works regardless of how many obs nodes were registered.
    OBS_NODE.OUT wraps its header line when there are many nodes, which breaks
    pandas; NOD_INF.OUT has a clean one-record-per-node layout at each step.
    """
    rows = []
    for t, df in sorted(nod.items()):
        z = df["Depth"].to_numpy(float)
        h = df["Head"].to_numpy(float)
        mask = ~(np.isnan(z) | np.isnan(h))
        z, h = z[mask], h[mask]
        order = np.argsort(z)
        z, h = z[order], h[order]

        if np.all(h > 0) or np.all(h < 0):
            rows.append((t, np.nan))
            continue

        idx = np.where(np.diff(np.signbit(h)))[0]
        if idx.size == 0:
            rows.append((t, np.nan))
            continue
        i = idx[-1]
        dh = h[i + 1] - h[i]
        z0 = z[i] if abs(dh) < 1e-12 else z[i] + (-h[i]) * (z[i + 1] - z[i]) / dh
        rows.append((t, z0))

    gwl = pd.DataFrame(rows, columns=["Time_d", "z_wt_cm"])
    gwl["depth_wt_cm"] = -gwl["z_wt_cm"]
    return gwl


# ---------------------------------------------------------------------------
# Stage 1 - constant recharge, zero bottom flux, NO roots, NO atmosphere
#
# Profile  : 300 cm deep, dx = 10 cm  →  31 nodes  (fast numerics)
# Soil     : van-Genuchten loam  (1/alpha ≈ 28 cm capillary fringe)
# Top BC   : constant flux  0.30 cm/day  (gentle recharge)
# Bottom BC: constant flux  0.00 cm/day  (sealed base - the "bucket")
# WT init  : z = -150 cm  (hydrostatic initial conditions)
# Run      : 15 days, daily output
# Obs nodes: 15 nodes  z = -80 … -220 cm  (straddles WT + fringe zone)
# ---------------------------------------------------------------------------

Z_TOP = 0.0
Z_BOT = -300.0
Z_WT0 = -150.0  # initial water-table elevation [cm]
DX = 10.0  # grid spacing [cm]
TMAX_STAGE1 = 15
TMAX_STAGE2 = 30
TMAX_STAGE3 = 30
RTOP = -0.30  # constant recharge rate [cm/day] — negative = downward into soil


def _build_base_model(ws: Path, title: str, tmax: int) -> tuple:
    ml = ps.Model(
        exe_name=resolve_hydrus_exe(),
        ws_name=str(ws),
        name="example_gw",
        description=title,
        mass_units="-",
        time_unit="days",
        length_unit="cm",
        print_screen=False,
    )

    ml.add_time_info(
        tinit=0,
        tmax=tmax,
        print_times=True,
        dt=0.05,
        dtmax=0.5,
        dtprint=1,
    )

    return ml


def _add_material_profile_obs(ml) -> pd.DataFrame:
    # Van-Genuchten loam: moderate capillary fringe, numerically well-behaved.
    m = ml.get_empty_material_df(n=1)
    m.loc[1] = [0.078, 0.43, 0.036, 1.56, 24.96, 0.5]
    ml.add_material(m)

    # 31-node profile with hydrostatic initial heads.
    profile = ps.create_profile(top=Z_TOP, bot=Z_BOT, dx=DX, mat=1)
    profile["h"] = hydrostatic_h(profile["x"].to_numpy(float), Z_WT0)
    ml.add_profile(profile)

    # 3 obs depths for head trajectories.
    ml.add_obs_nodes([-80.0, -150.0, -220.0])
    return profile


def build_atmospheric_forcing(tmax: int, include_root: bool = False) -> pd.DataFrame:
    days = np.arange(1, tmax + 1, dtype=float)

    # Positive values are required in ATMOSPH.IN (absolute-value convention).
    # We keep forcing smooth and moderate to avoid timestep collapse.
    prec = 0.05 + 0.10 * (0.5 + 0.5 * np.sin(2.0 * np.pi * days / 10.0))
    storm = (days % 7 == 0).astype(float) * 0.20
    prec = prec + storm

    rsoil = 0.03 + 0.04 * (0.5 + 0.5 * np.sin(2.0 * np.pi * (days - 2.0) / 9.0))
    if include_root:
        rroot = 0.02 + 0.04 * (0.5 + 0.5 * np.sin(2.0 * np.pi * (days + 1.0) / 11.0))
    else:
        rroot = np.zeros_like(days)

    return pd.DataFrame(
        {
            "tAtm": days,
            "Prec": prec,
            "rSoil": rsoil,
            "rRoot": rroot,
            "hCritA": np.full_like(days, 1e5),
            "rB": np.zeros_like(days),
            "hB": np.zeros_like(days),
            "ht": np.zeros_like(days),
        }
    )


def run_stage1(ws: Path) -> tuple:
    ml = _build_base_model(
        ws=ws,
        title="Stage 1: constant recharge, zero bottom flux",
        tmax=TMAX_STAGE1,
    )

    # top_bc=1: constant prescribed flux (rtop); bot_bc=1: zero flux (rbot=0).
    # KodBot=-1 (Neumann) for a zero-flux (sealed) base.
    ml.add_waterflow(
        top_bc=1,
        bot_bc=1,
        rtop=RTOP,
        rbot=0.0,
        rroot=0.0,
        maxit=20,
        tolh=1.0,
        ha=1e-6,
        hb=1e4,
    )

    profile = _add_material_profile_obs(ml)

    ml.write_input()
    result = ml.simulate()
    if result.returncode != 0:
        raise RuntimeError("HYDRUS run failed - see Error.msg in workspace.")

    return ml, profile


def run_stage2(ws: Path) -> tuple:
    ml = _build_base_model(
        ws=ws,
        title="Stage 2: atmospheric forcing, zero bottom flux",
        tmax=TMAX_STAGE2,
    )

    # top_bc=3 activates atmospheric boundary forcing through ATMOSPH.IN.
    ml.add_waterflow(
        top_bc=3,
        bot_bc=1,
        rbot=0.0,
        rroot=0.0,
        maxit=20,
        tolh=1.0,
        ha=1e-6,
        hb=1e4,
    )

    profile = _add_material_profile_obs(ml)

    atm = build_atmospheric_forcing(tmax=TMAX_STAGE2, include_root=False)
    atm.to_csv(ws / "atmospheric_forcing.csv", index=False)
    ml.add_atmospheric_bc(atm)

    ml.write_input()
    result = ml.simulate()
    if result.returncode != 0:
        raise RuntimeError("HYDRUS run failed in Stage 2 - see Error.msg.")
    return ml, profile


def run_stage3(ws: Path) -> tuple:
    ml = _build_base_model(
        ws=ws,
        title="Stage 3: atmospheric forcing + shallow root uptake",
        tmax=TMAX_STAGE3,
    )

    ml.add_waterflow(
        top_bc=3,
        bot_bc=1,
        rbot=0.0,
        rroot=0.0,
        maxit=20,
        tolh=1.0,
        ha=1e-6,
        hb=1e4,
    )

    profile = _add_material_profile_obs(ml)

    atm = build_atmospheric_forcing(tmax=TMAX_STAGE3, include_root=True)
    atm.to_csv(ws / "atmospheric_forcing.csv", index=False)
    ml.add_atmospheric_bc(atm)

    # Feddes uptake with shallow root-depth growth to avoid saturated-zone roots.
    ml.add_root_uptake(model=0, poptm=[-25], p0=-10, p2h=-200, p2l=-800, p3=-8000)
    ml.add_root_growth(
        irootin=2,
        irfak=1,
        trmin=0,
        trmed=0,
        trmax=TMAX_STAGE3,
        xrmin=5,
        xrmed=0,
        xrmax=60,
        trperiod=365,
    )

    ml.write_input()
    result = ml.simulate()
    if result.returncode != 0:
        raise RuntimeError("HYDRUS run failed in Stage 3 - see Error.msg.")
    return ml, profile


# ---------------------------------------------------------------------------
# Post-processing  (shared across stages)
# ---------------------------------------------------------------------------


def postprocess(ml, profile: pd.DataFrame, ws: Path, stage: int) -> None:
    tlevel = ml.read_tlevel()
    tlevel.to_csv(ws / "tlevel.csv")

    # GWL from full nodal snapshots (NOD_INF.OUT) - robust regardless of
    # how many obs nodes exist; avoids the OBS_NODE wrapped-header bug.
    nod = ml.read_nod_inf()
    # read_nod_inf returns a bare DataFrame when exactly one snapshot exists.
    if isinstance(nod, pd.DataFrame):
        nod = {0.0: nod}
    gwl = interpolate_gwl_from_nod_inf(nod)
    gwl.to_csv(ws / "gwl_timeseries.csv", index=False)

    # Save one representative NOD_INF snapshot (last printed time).
    last_t = max(nod.keys())
    nod[last_t].to_csv(ws / "nod_inf_last.csv", index=False)

    # Obs-node head trajectories extracted from NOD_INF (avoid read_obs_node
    # whose pandas c-engine parser breaks when OBS_NODE.OUT has a wide/wrapped
    # header line, which HYDRUS writes when many obs nodes are registered).
    OBS_Z = [-80.0, -150.0, -220.0]  # must match depths in run_stageN
    heads_list = []
    for z_target in OBS_Z:
        series = pd.Series(
            {t: float(df.loc[(df["Depth"] - z_target).abs().idxmin(), "Head"]) for t, df in nod.items()},
            name=f"h_z{int(abs(z_target))}cm",
        )
        heads_list.append(series)
    head_df = pd.concat(heads_list, axis=1)
    head_df.index.name = "Time_d"
    head_df.to_csv(ws / "obs_heads.csv")

    # Saturated base check: deepest node in every NOD_INF snapshot.
    bot_h = pd.Series(
        {t: float(df.loc[df["Depth"].abs().idxmax(), "Head"]) for t, df in nod.items()},
        name="h_bottom_cm",
    )
    bot_h.to_csv(ws / "bottom_head_check.csv", index_label="Time_d")
    bottom_head = bot_h

    # ── Plots ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"Stage {stage}: floating water-table (zero bottom flux)")

    # Panel 1: water-table trajectory
    ax = axes[0]
    ax.plot(gwl["Time_d"], gwl["depth_wt_cm"], lw=2, color="tab:blue")
    ax.invert_yaxis()
    ax.set(xlabel="Time [days]", ylabel="WT depth below surface [cm]", title="Water-table trajectory")
    ax.grid(alpha=0.3)

    # Panel 2: pressure-head traces at three depths (from NOD_INF via head_df)
    ax = axes[1]
    for col in head_df.columns:
        label = col.replace("h_z", "z=−").replace("cm", " cm")
        ax.plot(head_df.index, head_df[col], label=label)
    ax.axhline(0, color="k", ls="--", lw=1, label="h=0 (WT)")
    ax.set(xlabel="Time [days]", ylabel="h [cm]", title="Pressure heads at selected depths")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    for col, label, color in [
        ("vTop", "vTop (actual)", "tab:blue"),
        ("vBot", "vBot (≈0 sealed base)", "tab:orange"),
    ]:
        if col in tlevel.columns:
            ax.plot(tlevel.index, tlevel[col], label=label, color=color)
    ax.axhline(0, color="k", lw=0.8)
    ax.set(xlabel="Time [days]", ylabel="Flux [cm/day]", title="Boundary fluxes")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(ws / f"stage{stage}_overview.png", dpi=150)
    plt.close(fig)

    # ── Analysis summary ────────────────────────────────────────────────────
    summary = pd.Series(
        {
            "stage": stage,
            "wt_initial_cm": Z_WT0,
            "wt_depth_start_cm": float(gwl["depth_wt_cm"].iloc[0]),
            "wt_depth_end_cm": float(gwl["depth_wt_cm"].iloc[-1]),
            "wt_depth_min_cm": float(gwl["depth_wt_cm"].min(skipna=True)),
            "wt_depth_max_cm": float(gwl["depth_wt_cm"].max(skipna=True)),
            "wt_excursion_cm": float(gwl["depth_wt_cm"].max(skipna=True) - gwl["depth_wt_cm"].min(skipna=True)),
            "wt_rise_from_init_cm": float(abs(Z_WT0) - gwl["depth_wt_cm"].min()),
            "bottom_h_min_cm": float(bottom_head.min(skipna=True)),
            "bottom_h_max_cm": float(bottom_head.max(skipna=True)),
            "bottom_stayed_sat": bool(bottom_head.min(skipna=True) > 0),
            "vBot_mean_cm_d": float(tlevel["vBot"].mean()) if "vBot" in tlevel.columns else float("nan"),
            "vBot_maxabs_cm_d": float(tlevel["vBot"].abs().max()) if "vBot" in tlevel.columns else float("nan"),
        }
    )
    summary.to_csv(ws / "analysis_summary.csv", header=["value"])
    print(summary.to_string())


def compare_stages(base_ws: Path) -> None:
    rows = []
    for stage in (1, 2, 3):
        ws = base_ws / f"stage{stage}"
        s = pd.read_csv(ws / "analysis_summary.csv", index_col=0).iloc[:, 0]
        s.name = f"stage{stage}"
        rows.append(s)

    compare = pd.DataFrame(rows)
    compare.to_csv(base_ws / "stage_comparison.csv")

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    for stage in (1, 2, 3):
        ws = base_ws / f"stage{stage}"
        gwl = pd.read_csv(ws / "gwl_timeseries.csv")
        ax.plot(gwl["Time_d"], gwl["depth_wt_cm"], lw=2, label=f"Stage {stage}")
    ax.invert_yaxis()
    ax.set(xlabel="Time [days]", ylabel="WT depth below surface [cm]", title="WT trajectory by stage")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(base_ws / "stage_comparison_wt.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    base_ws = Path(__file__).resolve().parent / "example_gw"
    base_ws.mkdir(parents=True, exist_ok=True)

    print("=== Stage 1: constant recharge, zero bottom flux, no roots ===")
    ws1 = base_ws / "stage1"
    ws1.mkdir(parents=True, exist_ok=True)
    ml1, profile1 = run_stage1(ws1)
    postprocess(ml1, profile1, ws1, stage=1)

    print("=== Stage 2: atmospheric forcing, zero bottom flux, no roots ===")
    ws2 = base_ws / "stage2"
    ws2.mkdir(parents=True, exist_ok=True)
    ml2, profile2 = run_stage2(ws2)
    postprocess(ml2, profile2, ws2, stage=2)

    print("=== Stage 3: atmospheric forcing + shallow roots ===")
    ws3 = base_ws / "stage3"
    ws3.mkdir(parents=True, exist_ok=True)
    ml3, profile3 = run_stage3(ws3)
    postprocess(ml3, profile3, ws3, stage=3)

    compare_stages(base_ws)
    print(f"\nOutputs written to: {base_ws}")


if __name__ == "__main__":
    main()
