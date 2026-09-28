"""
Exploratory analysis of the relative size target z = log(peak + eps) - log(M_t + eps) used by the peak_* models.

Training rows are built with idmodels.peak.series.build_replay_rows, exactly as in the models, from final ILINet
(x percent positive) and FluSurv-NET data (the cached parquet written by validate_ilinet.py). NHSN is deliberately
not loaded: its seasons are the held-out test data. The ILINet series for Puerto Rico (72) and the US Virgin Islands
(78) are dropped before building rows (they are identically zero). Figures (PNG), tables (markdown) and a json of the numbers quoted in relative-size-eda.qmd are
written to analysis/peak-models/eda/.

Usage (from the repository root, with an environment that has idmodels, pandas and matplotlib):
    DYLD_FALLBACK_LIBRARY_PATH=<venv>/lib/python3.12/site-packages/sklearn/.dylibs \
        python analysis/peak-models/relative_size_eda.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from idmodels.peak.base import season_week_to_date  # noqa: E402
from idmodels.peak.series import LOG_EPS, build_replay_rows, build_season_arrays, running_max  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "eda"
ILI_PATH = HERE / "ilinet-validation" / "data-ilinet-flusurvnet.parquet"
DROP_ILINET_LOCATIONS = ["72", "78"]  # Puerto Rico, US Virgin Islands: all-zero series

W0, W1, REPLAY_START, MIN_OBS = 10, 43, 5, 25
ZERO_TOL = 1e-9

# source groups, in plotting order, with colors from a colorblind-validated categorical palette
GROUPS = ["ILINet states", "ILINet national + regions", "FluSurv-NET"]
COLORS = {"ILINet states": "#2a78d6", "ILINet national + regions": "#eb6834", "FluSurv-NET": "#1baf7a"}
ACCENT = "#eb6834"  # non-source accent (running max, medians) in single-source figures
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
Z_YLIM = (-0.3, 9.5)

plt.rcParams.update({
    "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"], "font.size": 10, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "axes.labelcolor": INK2, "axes.edgecolor": "#b7b6b0", "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 160,
    "savefig.bbox": "tight",
})


# ---------------------------------------------------------------------------------------------------------------
# data

def group_of(source: pd.Series, agg_level: pd.Series) -> pd.Series:
    g = np.select(
        [(source == "ilinet") & (agg_level == "state"), source == "ilinet", source == "flusurvnet"],
        ["ILINet states", "ILINet national + regions", "FluSurv-NET"], "other")
    return pd.Series(g, index=source.index)


def load_rows():
    data = pd.read_parquet(ILI_PATH)
    data = data.loc[~((data["source"] == "ilinet") & data["location"].isin(DROP_ILINET_LOCATIONS))]
    arrays = build_season_arrays(data)
    rows = build_replay_rows(arrays, W0, W1, REPLAY_START, MIN_OBS)
    rows["group"] = group_of(rows["source"], rows["agg_level"])
    rows["eps"] = rows["source"].map(LOG_EPS)
    rows["M"] = np.exp(rows["lm"]) - rows["eps"]
    rows["at_zero"] = rows["z"].abs() < ZERO_TOL
    rows["log_peak"] = rows["z"] + rows["lm"]
    return arrays, rows


# ---------------------------------------------------------------------------------------------------------------
# helpers

def week_axis(ax, label=True):
    """Season-week x axis with approximate month labels (dates shift by a few days between seasons)."""
    ticks = [5, 10, 15, 20, 25, 30, 35, 40]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{w}\n{season_week_to_date('2017/18', w):%b}" for w in ticks])
    ax.set_xlim(REPLAY_START - 0.5, W1 + 0.5)
    ax.axvspan(REPLAY_START - 0.5, W0 - 0.5, color="#f1f0ec", zorder=0, lw=0)
    if label:
        ax.set_xlabel("Season week t (approx. month)")


def qtab(df, by, col="z", qs=(0.1, 0.25, 0.5, 0.75, 0.9)):
    g = df.groupby(by)[col]
    out = pd.concat({f"q{int(q * 100)}": g.quantile(q) for q in qs}, axis=1)
    out["n"] = g.size()
    return out


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)


def md_table(df: pd.DataFrame, path: Path, floatfmt=2):
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "|".join("---:" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = [f"{v:.{floatfmt}f}" if isinstance(v, (float, np.floating)) else str(v) for v in r]
        lines.append("| " + " | ".join(cells) + " |")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# figures

def fig_examples(arrays, rows, examples, name, title):
    """For each example series: top the series, running max M_t and the peak; bottom z_t."""
    fig, axes = plt.subplots(2, len(examples), figsize=(3.3 * len(examples), 5.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1]})
    for j, (src, loc, season, lab, unit) in enumerate(examples):
        i = np.flatnonzero((arrays.keys["source"] == src) & (arrays.keys["location"] == loc) &
                           (arrays.keys["season"] == season))[0]
        y = arrays.y[i]
        weeks = np.arange(1, len(y) + 1)
        m = np.array([running_max(arrays.y[[i]], t, W0)[0][0] for t in range(1, W1 + 1)])
        r = rows.query("source == @src and location == @loc and season == @season").sort_values("season_week")
        pk_w = int(r["peak_week"].iloc[0])
        ax = axes[0, j]
        eps = LOG_EPS[src]
        sel = (weeks >= REPLAY_START) & (weeks <= W1)
        pk = y[pk_w - 1]
        ax.plot(weeks[sel], y[sel] + eps, color=INK2, lw=1.2, marker="o", ms=2.5, label="weekly value")
        tt = np.arange(1, W1 + 1)[REPLAY_START - 1:]
        ax.step(tt, m[REPLAY_START - 1:] + eps, where="post", color=ACCENT, lw=2, label="running max $M_t$")
        ax.axhline(pk + eps, color=COLORS["ILINet states"], lw=1, ls="--")
        ax.plot([pk_w], [pk + eps], marker="*", ms=12, color=COLORS["ILINet states"], mec="white", zorder=5)
        ax.annotate(f"peak, week {pk_w}", (pk_w, pk + eps), xytext=(6, 4), textcoords="offset points",
                    fontsize=8.5, color=INK)
        # z at one week is the vertical gap between M_t and the peak on this log scale
        ta = 18
        ax.annotate("", (ta, pk + eps), (ta, m[ta - 1] + eps),
                    arrowprops=dict(arrowstyle="<->", color=COLORS["ILINet states"], lw=1.3))
        za = np.log(pk + eps) - np.log(m[ta - 1] + eps)
        ax.text(ta - 0.6, np.sqrt((pk + eps) * (m[ta - 1] + eps)), f"$z_{{{ta}}}$ = {za:.1f}", ha="right",
                va="center", fontsize=8.5, color=COLORS["ILINet states"], fontweight="bold")
        ax.set_yscale("log")
        ax.set_ylim(eps * 0.6, (pk + eps) * 3)
        ax.set_title(lab, loc="left")
        ax.set_ylabel(f"{unit} + eps (log scale)", fontsize=9)
        week_axis(ax, label=False)
        ax = axes[1, j]
        ax.plot(r["season_week"], r["z"], color=COLORS["ILINet states"], lw=2, marker="o", ms=3)
        ax.axvline(pk_w, color=COLORS["ILINet states"], lw=1, ls="--")
        ax.set_ylim(*Z_YLIM)
        if j == 0:
            ax.set_ylabel("eventual z = log(peak) − log($M_t$)", fontsize=9)
        week_axis(ax)
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper right", ncol=2, fontsize=9, bbox_to_anchor=(0.99, 1.0))
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, name)


def fig_realtime(arrays, src="ilinet", loc="25", season="2017/18", t=18):
    """What a forecaster sees at week t (left) versus what is known only after the season (right)."""
    i = np.flatnonzero((arrays.keys["source"] == src) & (arrays.keys["location"] == loc) &
                       (arrays.keys["season"] == season))[0]
    eps = LOG_EPS[src]
    y = arrays.y[i] + eps
    weeks = np.arange(1, len(y) + 1)
    m = np.array([running_max(arrays.y[[i]], tt, W0)[0][0] for tt in range(1, W1 + 1)]) + eps
    pk_w = int(W0 + np.nanargmax(arrays.y[i, W0 - 1:W1]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, hindsight in zip(axes, [False, True]):
        sel = (weeks >= REPLAY_START) & (weeks <= t)
        ax.plot(weeks[sel], y[sel], color=INK2, lw=1.4, marker="o", ms=3)
        ax.step(np.arange(REPLAY_START, t + 1), m[REPLAY_START - 1:t], where="post", color=ACCENT, lw=2)
        ax.axvline(t, color=INK2, lw=1, ls=":")
        ax.text(t + 0.4, 1.3e-3, f"today: week {t}", fontsize=8.5, color=INK2)
        ax.text(t - 0.6, m[t - 1] * 1.25, f"$M_{{{t}}}$ (known)", fontsize=8.5, color=ACCENT, ha="right")
        if hindsight:
            fut = (weeks > t) & (weeks <= W1)
            ax.plot(weeks[fut], y[fut], color="#b7b6b0", lw=1.4, marker="o", ms=3)
            ax.axhline(y[pk_w - 1], color=COLORS["ILINet states"], lw=1, ls="--")
            ax.plot([pk_w], [y[pk_w - 1]], marker="*", ms=12, color=COLORS["ILINet states"], mec="white", zorder=5)
            ax.annotate("", (t, y[pk_w - 1]), (t, m[t - 1]),
                        arrowprops=dict(arrowstyle="<->", color=COLORS["ILINet states"], lw=1.4))
            zt = np.log(y[pk_w - 1]) - np.log(m[t - 1])
            ax.text(t - 0.6, np.sqrt(y[pk_w - 1] * m[t - 1]), f"$z_{{{t}}}$ = {zt:.1f}\n(computed after\nthe season ends)",
                    fontsize=8.5, color=COLORS["ILINet states"], va="center", ha="right")
            ax.text(pk_w + 0.6, y[pk_w - 1] * 1.3, f"actual peak, week {pk_w}", fontsize=8.5, color=INK)
            ax.set_title("In hindsight: the rest of the season, the peak, and z", loc="left")
        else:
            ax.text(30, 0.3, "future weeks,\npeak and z:\nunknown", fontsize=10, color=INK2, ha="center",
                    va="center", style="italic")
            ax.set_title("In real time: what a forecaster knows at week t", loc="left")
        ax.set_yscale("log")
        ax.set_ylim(eps * 0.8, np.nanmax(y) * 4)
        week_axis(ax)
    axes[0].set_ylabel("ILI × % positive + eps (log scale)", fontsize=9)
    fig.suptitle("Massachusetts, ILINet 2017/18", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig00_realtime_vs_hindsight")


def fig_fan(rows):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for ax, g in zip(axes, GROUPS):
        q = qtab(rows[rows["group"] == g], "season_week")
        c = COLORS[g]
        ax.fill_between(q.index, q["q10"], q["q90"], color=c, alpha=0.18, lw=0, label="10–90%")
        ax.fill_between(q.index, q["q25"], q["q75"], color=c, alpha=0.38, lw=0, label="25–75%")
        ax.plot(q.index, q["q50"], color=c, lw=2, label="median")
        n_ser = rows[rows["group"] == g].groupby(["location", "season"]).ngroups
        ax.set_title(f"{g}\n({n_ser} series)", loc="left")
        ax.set_ylim(*Z_YLIM)
        week_axis(ax)
    axes[0].set_ylabel("eventual z (known only after the season)")
    axes[0].legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig03_z_fan_by_source")


def fig_atom(rows):
    fig, ax = plt.subplots(figsize=(7.5, 4))
    frac = rows.groupby(["group", "season_week"])["at_zero"].mean().unstack(0)
    for g in GROUPS:
        ax.plot(frac.index, frac[g], color=COLORS[g], lw=2, label=g)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_ylim(0, 1.02)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_yticklabels([f"{v:.0%}" for v in np.linspace(0, 1, 6)])
    ax.set_ylabel("share of rows with z = 0\n(peak already reached)")
    week_axis(ax)
    fig.tight_layout()
    save(fig, "fig04_atom_at_zero")
    return frac


def quantile_bars(ax, df, bycol, order, color):
    """Median + 25–75% + 10–90% ranges per category, with the share at z = 0 printed above."""
    for x, cat in enumerate(order):
        v = df.loc[df[bycol] == cat]
        if len(v) < 10:
            continue
        q = v["z"].quantile([0.1, 0.25, 0.5, 0.75, 0.9]).to_numpy()
        ax.plot([x, x], [q[0], q[4]], color=color, lw=1.5, alpha=0.5, solid_capstyle="round")
        ax.plot([x, x], [q[1], q[3]], color=color, lw=7, alpha=0.8, solid_capstyle="butt")
        ax.plot([x - 0.22, x + 0.22], [q[2], q[2]], color=INK, lw=2)
        ax.text(x, Z_YLIM[1] * 0.93, f"{v['at_zero'].mean():.0%}\nn={len(v)}", ha="center", va="top",
                fontsize=7.5, color=INK2)
    ax.set_xticks(range(len(order)))
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.grid(axis="x", visible=False)


def fig_stage(rows, weeks=(16, 22, 28)):
    rows = rows.copy()
    rows["wsm_bin"] = np.where(rows["wks_since_max"] >= 4, "4+", rows["wks_since_max"].astype(int).astype(str))
    rows["rel_bin"] = pd.cut(rows["rel_max"], [-np.inf, -1.0, -0.4, -0.1, -1e-9, 0.1],
                             labels=["< −1", "−1 to −0.4", "−0.4 to −0.1", "−0.1 to 0", "at max (0)"])
    rows["g1_bin"] = pd.cut(rows["g1"], [-np.inf, -0.15, 0.15, 0.4, np.inf],
                            labels=["falling\n< −0.15", "flat\n±0.15", "rising\n0.15–0.4", "fast\n> 0.4"])
    specs = [("wsm_bin", ["0", "1", "2", "3", "4+"], "weeks since the max was set"),
             ("rel_bin", ["at max (0)", "−0.1 to 0", "−0.4 to −0.1", "−1 to −0.4", "< −1"],
              "rel_max = log(current / $M_t$)"),
             ("g1_bin", ["falling\n< −0.15", "flat\n±0.15", "rising\n0.15–0.4", "fast\n> 0.4"],
              "g1 = log growth over the last week")]
    fig, axes = plt.subplots(3, len(weeks), figsize=(12, 10.5), sharey=True)
    for i, (col, order, lab) in enumerate(specs):
        for j, wk in enumerate(weeks):
            ax = axes[i, j]
            quantile_bars(ax, rows[rows["season_week"] == wk], col, order, COLORS["ILINet states"])
            ax.set_xticklabels(order, fontsize=8)
            ax.set_ylim(*Z_YLIM)
            if i == 0:
                ax.set_title(f"Season week {wk} (~{season_week_to_date('2017/18', wk):%b %d})", loc="left")
            if j == 0:
                ax.set_ylabel("eventual z")
            ax.set_xlabel(lab)
    fig.tight_layout()
    save(fig, "fig05_z_by_stage")


def fig_kz(rows):
    d = rows[rows["k"] >= 1]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharex=True, sharey=True)
    stats = {}
    for ax, g in zip(axes, GROUPS):
        v = d[d["group"] == g]
        ax.hexbin(np.log(v["k"]), v["z"], gridsize=(22, 18), extent=(0, np.log(38), 0, Z_YLIM[1]), mincnt=1,
                  bins="log", cmap="Blues", linewidths=0.2)
        med = v.groupby("k")["z"].median()
        med = med[v.groupby("k").size() >= 20]
        ax.plot(np.log(med.index), med.values, color=INK, lw=2, label="median z given k")
        rho = v[["k", "z"]].corr(method="spearman").iloc[0, 1]
        stats[g] = rho
        ax.set_title(f"{g}\nSpearman ρ(k, z) = {rho:.2f}", loc="left")
        kt = [1, 2, 4, 8, 16, 32]
        ax.set_xticks(np.log(kt))
        ax.set_xticklabels(kt)
        ax.set_xlabel("k = weeks until the peak (log scale)")
    axes[0].set_ylabel("eventual z")
    axes[0].set_ylim(0, Z_YLIM[1])
    axes[0].legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig06_z_vs_k")
    return stats


def fig_compare(rows):
    """ILINet vs FluSurv-NET for the same (location, season) pairs: FluSurv-NET sites that are states, and the US."""
    fsn = rows.loc[rows["source"] == "flusurvnet", ["location", "season"]].drop_duplicates()
    ili = rows.loc[rows["source"] == "ilinet", ["location", "season"]].drop_duplicates()
    both = fsn.merge(ili)
    d = rows.merge(both)
    panels = [("FluSurv-NET sites (states) and the matching ILINet states", d[d["location"] != "US"]),
              ("United States", d[d["location"] == "US"])]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    out = {}
    for ax, (title, dd) in zip(axes, panels):
        for src, c, lab in [("ilinet", COLORS["ILINet states"], "ILINet"), ("flusurvnet", COLORS["FluSurv-NET"], "FluSurv-NET")]:
            q = qtab(dd[dd["source"] == src], "season_week")
            ax.fill_between(q.index, q["q25"], q["q75"], color=c, alpha=0.25, lw=0)
            ax.plot(q.index, q["q50"], color=c, lw=2)
            ax.plot(q.index, q["q90"], color=c, lw=1, ls=":")
            ax.plot([], [], color=c, lw=2, label=f"{lab} ({dd.loc[dd['source'] == src].groupby(['location', 'season']).ngroups} series)")
            out[(title, src)] = q
        ax.legend(loc="upper right", fontsize=8.5)
        ax.set_title(title, loc="left")
        ax.set_ylim(*Z_YLIM)
        week_axis(ax)
    axes[0].set_ylabel("eventual z: median, 25–75% (band), 90% (dotted)")
    fig.tight_layout()
    save(fig, "fig07_ilinet_vs_flusurvnet")
    # per-pair agreement: correlation of z between the two sources at the same location, season and week
    w = d.pivot_table(index=["location", "season", "season_week"], columns="source", values="z").dropna()
    w = w[w.index.get_level_values("season_week") >= W0]
    stats = {"n_pairs": int(both.shape[0]), "n_locations": int(both["location"].nunique()),
             "spearman_same_week": round(float(w.corr(method="spearman").iloc[0, 1]), 2),
             "median_abs_diff": round(float((w["ilinet"] - w["flusurvnet"]).abs().median()), 2),
             "median_z": {wk: {s: round(float(d.loc[(d["location"] != "US") & (d["source"] == s) &
                                                    (d["season_week"] == wk), "z"].median()), 2)
                               for s in ["ilinet", "flusurvnet"]} for wk in [10, 15, 20, 25]},
             "peak_week_diff_median": round(float(
                 (d[d["source"] == "ilinet"].drop_duplicates(["location", "season"]).set_index(["location", "season"])["peak_week"]
                  - d[d["source"] == "flusurvnet"].drop_duplicates(["location", "season"]).set_index(["location", "season"])["peak_week"]
                  ).median()), 1)}
    return stats


def fig_scale(rows, week=18):
    peaks = rows.drop_duplicates(["source", "agg_level", "location", "season"])
    zw = rows[rows["season_week"] == week]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    rng = np.random.default_rng(1)
    for ax, (d, col, lab) in zip(axes, [(peaks, "log_peak", "log(peak + eps), in each source's own units"),
                                        (zw, "z", f"z at season week {week}")]):
        for i, g in enumerate(GROUPS):
            v = d.loc[d["group"] == g, col].to_numpy()
            ax.scatter(v, i + rng.uniform(-0.18, 0.18, len(v)), s=6, color=COLORS[g], alpha=0.35, lw=0)
            q = np.quantile(v, [0.25, 0.5, 0.75])
            ax.plot([q[0], q[2]], [i + 0.32, i + 0.32], color=INK, lw=3, solid_capstyle="butt")
            ax.plot([q[1]], [i + 0.32], marker="|", ms=10, mew=2, color="white")
        ax.set_yticks(range(len(GROUPS)))
        ax.set_yticklabels(GROUPS if ax is axes[0] else [""] * len(GROUPS))
        ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
        ax.set_xlabel(lab)
    axes[0].set_title("Absolute peak size: not comparable across sources", loc="left")
    axes[1].set_title("Relative size z: on a common scale", loc="left")
    fig.tight_layout()
    save(fig, "fig08_absolute_vs_relative")
    return (peaks.groupby("group")["log_peak"].median(), zw.groupby("group")["z"].median())


def fig_eps(arrays, rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    near = rows.assign(near=rows["M"] < 10 * rows["eps"]).groupby(["group", "season_week"])["near"].mean().unstack(0)
    for g in GROUPS:
        ax.plot(near.index, near[g], color=COLORS[g], lw=2, label=g)
    ax.legend(loc="upper right", fontsize=8.5)
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_ylabel("share of rows with $M_t$ < 10·eps")
    ax.set_title("How often is the running max tiny?", loc="left")
    week_axis(ax)

    # ILINet (states + national/regions): recompute z with eps scaled by 1/10, 1, 10
    ax = axes[1]
    ili = rows[rows["source"] == "ilinet"]
    peak = np.exp(ili["log_peak"]) - ili["eps"]
    res = {}
    for f, ls in [(0.1, ":"), (1, "-"), (10, "--")]:
        e = ili["eps"] * f
        z = np.log(peak + e) - np.log(ili["M"] + e)
        q = z.groupby(ili["season_week"]).quantile([0.5, 0.9]).unstack()
        res[f] = q
        ax.plot(q.index, q[0.5], color=INK, lw=2, ls=ls)
        ax.plot(q.index, q[0.9], color="#9a9993", lw=2, ls=ls)
        ax.plot([], [], color=INK2, lw=1.5, ls=ls, label=f"eps × {f:g}")
    ax.set_ylim(0, 12)
    ax.set_ylabel("z (ILINet, all levels)")
    ax.set_title("Sensitivity of ILINet z to eps", loc="left")
    ax.plot([], [], color=INK, lw=2, label="median")
    ax.plot([], [], color="#9a9993", lw=2, label="90th percentile")
    ax.legend(fontsize=8.5, ncol=2, loc="upper right")
    week_axis(ax)
    fig.tight_layout()
    save(fig, "fig09_eps")
    return near, res


def fig_seasons(rows):
    d = rows[rows["group"] == "ILINet states"]
    med = d.groupby(["season", "season_week"])["z"].median().unstack(0)
    frac = d.groupby(["season", "season_week"])["at_zero"].mean().unstack(0)
    pooled = d.groupby("season_week")["z"].median()
    pk = d.drop_duplicates(["location", "season"]).groupby("season")["peak_week"].median()
    early, late = pk.idxmin(), pk.idxmax()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, tab, ylab in [(axes[0], med, "median z across states"),
                          (axes[1], frac, "share of states already peaked")]:
        for s in tab.columns:
            hi = s in (early, late)
            c = COLORS["ILINet states"] if s == early else ACCENT if s == late else "#b7b6b0"
            ax.plot(tab.index, tab[s], color=c, lw=2.2 if hi else 1, zorder=3 if hi else 2)
            if hi:
                wk = (14 if s == early else 25) if ax is axes[0] else (17 if s == early else 33)
                ax.annotate(f"{s} (median peak wk {pk[s]:.0f})", (wk, tab.loc[wk, s]), xytext=(6, 2),
                            textcoords="offset points", fontsize=8, color=c, fontweight="bold",
                            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))
        ax.set_ylabel(ylab)
        week_axis(ax)
    axes[0].plot(pooled.index, pooled.values, color=INK, lw=1.5, ls="--", zorder=4)
    axes[0].annotate("all seasons pooled", (11, pooled.loc[11]), xytext=(10, 14), textcoords="offset points",
                     fontsize=8, color=INK)
    axes[0].set_ylim(*Z_YLIM)
    axes[1].set_ylim(0, 1.02)
    axes[1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    axes[0].set_title("Each line is one season (ILINet, states)", loc="left")
    axes[1].set_title("Timing differs by season, and z follows", loc="left")
    fig.tight_layout()
    save(fig, "fig10_season_spaghetti")
    # between-season spread of median z at a few weeks
    return med, pk


# ---------------------------------------------------------------------------------------------------------------

def main():
    OUT.mkdir(exist_ok=True)
    arrays, rows = load_rows()
    rows.to_parquet(OUT / "replay_rows.parquet")
    nums = {"n_rows": len(rows), "n_series": int(rows.groupby(["source", "agg_level", "location", "season"]).ngroups),
            "n_negative_z": int((rows["z"] < -ZERO_TOL).sum())}

    # example series: Massachusetts, a small state, and the US
    small = next(loc for loc in ["50", "56", "30", "44", "33"]
                 if all(((rows["source"] == s) & (rows["location"] == loc) & (rows["season"] == ss)).any()
                        for s, ss in [("ilinet", "2017/18"), ("ilinet", "2019/20"), ("ilinet", "2012/13")]))
    small_name = {"50": "Vermont", "56": "Wyoming", "30": "Montana", "44": "Rhode Island", "33": "New Hampshire"}[small]
    nums["small_state"] = small_name
    ili_u, fsn_u = "ILI × % positive", "hospitalizations per 100k"
    for loc, nm, fname in [("25", "Massachusetts", "fig01_examples_ma"), (small, small_name, "fig02_examples_small")]:
        fig_examples(arrays, rows, [("ilinet", loc, "2017/18", f"{nm}, ILINet 2017/18", ili_u),
                                    ("ilinet", loc, "2019/20", f"{nm}, ILINet 2019/20", ili_u),
                                    ("ilinet", loc, "2012/13", f"{nm}, ILINet 2012/13", ili_u)],
                     fname, f"{nm}: the series, its running maximum, and z")
    fig_examples(arrays, rows, [("ilinet", "US", "2017/18", "US, ILINet 2017/18", ili_u),
                                ("ilinet", "US", "2019/20", "US, ILINet 2019/20", ili_u),
                                ("flusurvnet", "US", "2017/18", "US, FluSurv-NET 2017/18", fsn_u)],
                 "fig02b_examples_us", "United States: the series, its running maximum, and z")

    fig_realtime(arrays)
    fig_fan(rows)
    frac = fig_atom(rows)
    fig_stage(rows)
    rho = fig_kz(rows)
    nums["ilinet_vs_flusurvnet"] = fig_compare(rows)
    lp_med, zw_med = fig_scale(rows)
    near, eps_res = fig_eps(arrays, rows)
    season_med, season_pk = fig_seasons(rows)

    # tables
    weeks = [5, 10, 15, 20, 25, 30, 35, 40]
    q = qtab(rows, ["group", "season_week"])
    q["p0"] = rows.groupby(["group", "season_week"])["at_zero"].mean()
    tab = pd.DataFrame({"week": weeks})
    for g in GROUPS:
        cell = []
        for w in weeks:
            if (g, w) not in q.index:
                cell.append("–")
                continue
            r = q.loc[(g, w)]
            cell.append(f"{r['q50']:.2f} [{r['q25']:.2f}, {r['q75']:.2f}]; {r['p0']:.0%}")
        tab[g] = cell
    md_table(tab, OUT / "table_z_by_week.md")

    cov = rows.groupby("group").agg(series=("z", lambda s: 0), seasons=("season", "nunique"),
                                    locations=("location", "nunique"), rows=("z", "size"))
    cov["series"] = rows.groupby("group").apply(lambda d: d.groupby(["location", "season"]).ngroups,
                                                include_groups=False)
    cov["seasons_span"] = rows.groupby("group")["season"].agg(lambda s: f"{s.min()}–{s.max()}")
    cov = cov.reindex(GROUPS).reset_index().rename(columns={"group": "source"})
    md_table(cov[["source", "seasons_span", "seasons", "locations", "series", "rows"]], OUT / "table_coverage.md")

    # numbers quoted in the text
    med = rows.groupby(["group", "season_week"])["z"].median()
    nums["median_z"] = {g: {w: round(float(med.get((g, w), np.nan)), 2) for w in weeks} for g in GROUPS}
    nums["frac_zero"] = {g: {w: round(float(frac.loc[w, g]), 3) for w in [15, 20, 25, 30, 35, 40]} for g in GROUPS}
    nums["pooled_median_z"] = {w: round(float(rows.loc[rows["season_week"] == w, "z"].median()), 2) for w in weeks}
    nums["pooled_frac_zero"] = {w: round(float(rows.loc[rows["season_week"] == w, "at_zero"].mean()), 3)
                                for w in [15, 20, 25, 30, 35]}
    nums["spearman_k_z"] = {g: round(float(v), 2) for g, v in rho.items()}
    nums["median_log_peak"] = {g: round(float(v), 2) for g, v in lp_med.items()}
    nums["median_log_peak_as_value"] = {g: round(float(np.exp(v)), 3) for g, v in lp_med.items()}
    nums["median_z_week18"] = {g: round(float(v), 2) for g, v in zw_med.items()}
    nums["near_eps"] = {g: {w: round(float(near[g].get(w, np.nan)), 3) for w in [5, 10, 15, 20]} for g in GROUPS}
    nums["eps_sensitivity_median"] = {f: {w: round(float(q.loc[w, 0.5]), 2) for w in [5, 10, 15, 20, 25]}
                                      for f, q in eps_res.items()}
    nums["eps_sensitivity_q90"] = {f: {w: round(float(q.loc[w, 0.9]), 2) for w in [5, 10, 15, 20, 25]}
                                   for f, q in eps_res.items()}
    nums["season_median_z_range"] = {w: [round(float(season_med.loc[w].min()), 2),
                                         round(float(season_med.loc[w].max()), 2)] for w in [15, 20, 25]}
    nums["season_median_peak_week"] = {s: float(v) for s, v in season_pk.items()}
    # stage conditioning numbers at week 22
    w22 = rows[rows["season_week"] == 22]
    nums["week22_by_wsm"] = {str(b): {"frac_zero": round(float(v["at_zero"].mean()), 2),
                                      "median": round(float(v["z"].median()), 2), "q90": round(float(v["z"].quantile(0.9)), 2),
                                      "n": len(v)}
                             for b, v in w22.groupby(np.minimum(w22["wks_since_max"], 4))}
    nums["week22_by_g1"] = {str(b): {"frac_zero": round(float(v["at_zero"].mean()), 2),
                                     "median": round(float(v["z"].median()), 2), "n": len(v)}
                            for b, v in w22.groupby(pd.cut(w22["g1"], [-np.inf, -0.15, 0.15, 0.4, np.inf]), observed=True)}
    w16 = rows[rows["season_week"] == 16]
    nums["week16_by_wsm"] = {str(b): {"frac_zero": round(float(v["at_zero"].mean()), 2),
                                      "median": round(float(v["z"].median()), 2), "n": len(v)}
                             for b, v in w16.groupby(np.minimum(w16["wks_since_max"], 4))}
    kk = rows[(rows["k"] >= 5) & (rows["k"] <= 8)]["z"]
    nums["z_k5to8_q10_q90"] = [round(float(kk.quantile(0.1)), 2), round(float(kk.quantile(0.9)), 2)]
    kz = rows[rows["k"] >= 1]
    nums["median_z_by_k"] = {int(k): round(float(kz.loc[kz["k"] == k, "z"].median()), 2) for k in [1, 2, 4, 8, 16]}
    (OUT / "key_numbers.json").write_text(json.dumps(nums, indent=1, default=str))
    print(json.dumps(nums, indent=1, default=str))


if __name__ == "__main__":
    main()
