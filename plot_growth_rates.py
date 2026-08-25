"""
Plot the per-segment growth rates produced by calculate_growth_rates.py.

Reads growth_rates.csv and a group definition CSV (reactor -> group) and writes:

  - growth_rate_scatter/<unit>.png   growth rate vs time, one plot per reactor
  - growth_rate_scatter/all_reactors.png
                                     the same scatters as small multiples
  - growth_rates_by_group.png        every reactor side by side, ordered and
                                     coloured by group, plus the per-group summary

--reference-line draws a dashed horizontal rule at a given growth rate on every
plot, for comparing against a target or a published value.

Points are shaded by the R^2 of their fit, so a poor regression is visible
rather than silently averaged in. Use --min-r2 / --between-dilutions-only to
drop segments that are not exponential growth (see the README), and
--omit-reactors to leave a reactor out entirely.
"""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from scipy import stats as sstats

RATE_COL = "growth_rate_h-1"
UNIT_COL = "pioreactor_unit"

# --- palette -----------------------------------------------------------------
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#83827d"
GRID = "#e3e2de"

# Categorical slots, assigned to groups in fixed order (never cycled).
GROUP_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# Secondary encoding, so groups are never distinguished by colour alone.
GROUP_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
OTHER_COLOR = "#83827d"
REFERENCE_COLOR = "#4a3aa7"

# Sequential blue ramp (light -> dark) used for R^2.
R2_CMAP = LinearSegmentedColormap.from_list(
    "r2_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#0d366b"]
)


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=3, color=GRID)


def load_growth_rates(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = [UNIT_COL, RATE_COL, "start_time_h", "end_time_h", "r_squared", "segment_position"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}.\n"
            f"Found columns: {', '.join(map(str, df.columns))}"
        )
    # Each segment spans an interval; plot it at its midpoint.
    df["time_h"] = (df["start_time_h"] + df["end_time_h"]) / 2.0
    return df


def _normalise_unit(name: str) -> str:
    # Reactor ids are letter+digits; a typed capital O for a zero is a common slip.
    return str(name).strip().upper().replace("O", "0")


def load_groups(path: Path, units: list[str], omitted: set[str] | None = None) -> pd.DataFrame:
    """Map reactor -> group from a two-column CSV.

    The header of the supplied file has its labels the other way round from its
    contents, so the reactor column is identified by which one actually holds
    reactor ids; the header is only a fallback.
    """
    raw = pd.read_csv(path)
    if raw.shape[1] < 2:
        raise ValueError(f"{path} needs at least two columns (reactor and group).")

    known = {_normalise_unit(u): u for u in units}
    cols = list(raw.columns)
    hits = [sum(_normalise_unit(v) in known for v in raw[c]) for c in cols]
    best = int(np.argmax(hits))
    if hits[best] == 0:
        # Nothing matched, so fall back to the column names.
        lowered = [str(c).strip().lower() for c in cols]
        best = lowered.index("reactor") if "reactor" in lowered else 0
        print(f"WARNING: no reactor id in {path} matches the growth-rate data.")
    reactor_col = cols[best]
    group_col = cols[1] if best == 0 else cols[0]

    out = pd.DataFrame(
        {
            "raw_reactor": raw[reactor_col].astype(str).str.strip(),
            "group": raw[group_col].astype(str).str.strip(),
        }
    )
    out[UNIT_COL] = [known.get(_normalise_unit(r)) for r in out["raw_reactor"]]

    fixed = out[out[UNIT_COL].notna() & (out["raw_reactor"] != out[UNIT_COL])]
    for _, row in fixed.iterrows():
        print(f"NOTE: group file lists '{row['raw_reactor']}', matched to reactor {row[UNIT_COL]}")
    # A reactor the user omitted on purpose is not worth warning about.
    dropped = {_normalise_unit(u) for u in (omitted or ())}
    unmatched = [
        r for r in out.loc[out[UNIT_COL].isna(), "raw_reactor"] if _normalise_unit(r) not in dropped
    ]
    if unmatched:
        print(f"WARNING: group file reactors with no growth-rate data: {', '.join(unmatched)}")

    return out.dropna(subset=[UNIT_COL])[[UNIT_COL, "group"]].drop_duplicates(subset=[UNIT_COL])


def group_styles(groups: list[str]) -> dict:
    # Fixed-order assignment: a group keeps its colour however many are plotted.
    styles = {}
    for i, g in enumerate(groups):
        if i < len(GROUP_COLORS):
            styles[g] = {"color": GROUP_COLORS[i], "marker": GROUP_MARKERS[i]}
        else:
            styles[g] = {"color": OTHER_COLOR, "marker": "."}
    return styles


def rule_handles(references=None, show_average: bool = True) -> list[Line2D]:
    """Legend entries for the dashed rules drawn on the scatter plots."""
    handles = []
    if show_average:
        handles.append(
            Line2D([], [], color="#eb6834", linestyle="--", linewidth=2, label="reactor average")
        )
    if references:
        values = ", ".join(f"{v:g}" for v in references)
        handles.append(
            Line2D(
                [],
                [],
                color=REFERENCE_COLOR,
                linestyle=(0, (6, 3)),
                linewidth=1.6,
                label=f"reference ({values} h$^{{-1}}$)",
            )
        )
    return handles


def add_reference_lines(ax, values, label: bool = True) -> None:
    """Draw the user's --reference-line values as dashed horizontal rules.

    Each line is labelled with its own value, so it is identifiable without a
    legend and cannot be confused with the (orange) reactor-mean line.
    """
    for value in values or ():
        ax.axhline(value, color=REFERENCE_COLOR, linewidth=1.6, linestyle=(0, (6, 3)), zorder=4)
        if label:
            # Labelled inside the axes: outside would collide with the colorbar.
            ax.annotate(
                f"{value:g}",
                xy=(0.995, value),
                xycoords=ax.get_yaxis_transform(),
                xytext=(0, 3),
                textcoords="offset points",
                ha="right",
                va="bottom",
                fontsize=8,
                color=REFERENCE_COLOR,
                zorder=5,
            )


def scatter_one_reactor(ax, unit_df: pd.DataFrame, norm, references=None, show_average=True):
    style_axes(ax)
    sc = ax.scatter(
        unit_df["time_h"],
        unit_df[RATE_COL],
        c=unit_df["r_squared"],
        cmap=R2_CMAP,
        norm=norm,
        s=46,
        linewidths=1.2,
        edgecolors=SURFACE,  # 2px-equivalent surface ring on overlapping marks
        zorder=3,
    )
    ax.axhline(0, color=TEXT_MUTED, linewidth=1, zorder=2)
    mean = unit_df[RATE_COL].mean()
    if show_average:
        ax.axhline(mean, color="#eb6834", linewidth=2, linestyle="--", zorder=2)
    add_reference_lines(ax, references)
    return sc, mean


def plot_per_reactor(df: pd.DataFrame, out_dir: Path, norm, references=None, show_average=True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for unit, unit_df in df.groupby(UNIT_COL, sort=True):
        fig, ax = plt.subplots(figsize=(8, 4.6), facecolor=SURFACE)
        sc, mean = scatter_one_reactor(ax, unit_df, norm, references, show_average)
        ax.set_xlabel("Time since experiment start (h)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_ylabel("Growth rate $\\mu$ (h$^{-1}$)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_title(
            f"{unit} — growth rate per segment  "
            f"(n = {len(unit_df)}, average {mean:.3f} h$^{{-1}}$)",
            color=TEXT_PRIMARY,
            fontsize=12,
            loc="left",
        )
        cbar = fig.colorbar(sc, ax=ax, pad=0.02)
        cbar.set_label("$R^2$ of fit", color=TEXT_SECONDARY, fontsize=9)
        cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
        cbar.outline.set_visible(False)
        legend_handles = rule_handles(references, show_average)
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                loc="best",
                fontsize=9,
                frameon=False,
                labelcolor=TEXT_SECONDARY,
            )
        fig.tight_layout()
        fig.savefig(out_dir / f"{unit}.png", dpi=150, facecolor=SURFACE)
        plt.close(fig)


def _grid_layout(units: list[str], groups: pd.DataFrame | None):
    """Rows of units for the small-multiples grid.

    With a group file, each group gets its own row (reactors with no group fall
    into a trailing row of their own); without one, units simply wrap at three
    per row.
    """
    if groups is None or groups.empty:
        ncols = min(3, len(units))
        return [units[i : i + ncols] for i in range(0, len(units), ncols)], None

    unit_group = groups.set_index(UNIT_COL)["group"].to_dict()
    rows, labels = [], []
    for g in sorted({unit_group[u] for u in units if u in unit_group}):
        rows.append(sorted(u for u in units if unit_group.get(u) == g))
        labels.append(g)
    ungrouped = sorted(u for u in units if u not in unit_group)
    if ungrouped:
        rows.append(ungrouped)
        labels.append("(no group)")
    return rows, labels


def plot_small_multiples(df, out_path: Path, norm, references=None, groups=None, show_average=True) -> None:
    units = sorted(df[UNIT_COL].unique())
    rows, row_labels = _grid_layout(units, groups)
    nrows = len(rows)
    ncols = max(len(r) for r in rows)
    styles = group_styles(row_labels) if row_labels else {}

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.6 * ncols, 3.6 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
        facecolor=SURFACE,
    )
    sc = None
    for r, row_units in enumerate(rows):
        for c in range(ncols):
            ax = axes[r][c]
            if c >= len(row_units):
                ax.set_visible(False)
                continue
            unit = row_units[c]
            sc, _ = scatter_one_reactor(ax, df[df[UNIT_COL] == unit], norm, references, show_average)
            ax.set_title(unit, color=TEXT_PRIMARY, fontsize=11, loc="left")
        if row_labels:
            # One group heading per row, above that row's first panel.
            label = row_labels[r]
            axes[r][0].annotate(
                f"{label}  ({len(row_units)} reactor{'s' if len(row_units) != 1 else ''})",
                xy=(0.0, 1.14),
                xycoords="axes fraction",
                fontsize=12,
                color=styles.get(label, {}).get("color", TEXT_SECONDARY),
                va="bottom",
            )

    # A partly-filled row would otherwise leave the column above it with no x
    # labels, so label the lowest visible panel of every column.
    for col in range(ncols):
        visible = [axes[row][col] for row in range(nrows) if axes[row][col].get_visible()]
        if visible:
            visible[-1].set_xlabel("Time (h)", color=TEXT_SECONDARY, fontsize=9)
            visible[-1].tick_params(labelbottom=True)
    for row in axes:
        row[0].set_ylabel("$\\mu$ (h$^{-1}$)", color=TEXT_SECONDARY, fontsize=9)

    fig.suptitle(
        "Growth rate per segment over time, by reactor"
        + (", grouped" if row_labels else ""),
        color=TEXT_PRIMARY,
        fontsize=13,
        x=0.01,
        y=0.99,
        ha="left",
    )
    # Row headings sit above each row's panels, so the rows need room between them.
    fig.subplots_adjust(
        top=0.90 - 0.02 * nrows if row_labels else 0.93,
        hspace=0.42 if row_labels else 0.22,
        wspace=0.16,
    )
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes, pad=0.015, fraction=0.02)
        cbar.set_label("$R^2$ of fit", color=TEXT_SECONDARY, fontsize=9)
        cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
        cbar.outline.set_visible(False)
    # The panels carry no legend of their own, so name the dashed rules here.
    legend_handles = rule_handles(references, show_average)
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            loc="upper right",
            bbox_to_anchor=(0.995, 0.995),
            fontsize=9,
            frameon=False,
            labelcolor=TEXT_SECONDARY,
        )
    fig.savefig(out_path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def _box(ax, values_by_x: list[np.ndarray], positions, colors, width=0.62):
    bp = ax.boxplot(
        values_by_x,
        positions=positions,
        widths=width,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color=SURFACE, linewidth=2),
        whiskerprops=dict(color=TEXT_MUTED, linewidth=1),
        capprops=dict(color=TEXT_MUTED, linewidth=1),
        boxprops=dict(linewidth=1.2, edgecolor=SURFACE),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    return bp


def compare_groups(merged: pd.DataFrame, level: str = "reactor") -> dict:
    """Test whether the groups differ in growth rate.

    At level="reactor" (the default) each reactor contributes one value, its
    average growth rate. That is the honest unit of replication: segments from
    the same reactor are repeated measures of one culture, so pooling them
    (level="segment") treats n as far larger than the experiment supports and
    makes any difference look better resolved than it is.

    Two groups are compared with Welch's t-test; more than two get a one-way
    ANOVA across all groups plus Holm-adjusted pairwise Welch tests.
    """
    if level == "reactor":
        per_group = {
            g: sub.groupby(UNIT_COL)[RATE_COL].mean().to_numpy()
            for g, sub in merged.groupby("group")
        }
        unit_name = "reactors"
    else:
        per_group = {g: sub[RATE_COL].to_numpy() for g, sub in merged.groupby("group")}
        unit_name = "segments"

    groups = sorted(per_group)
    testable = [g for g in groups if len(per_group[g]) >= 2]
    result = {
        "level": level,
        "unit_name": unit_name,
        "n": {g: len(per_group[g]) for g in groups},
        "mean": {g: float(np.mean(per_group[g])) for g in groups},
        "omnibus": None,
        "pairs": [],
        "skipped": [g for g in groups if g not in testable],
    }
    if len(testable) < 2:
        return result

    if len(testable) > 2:
        f_stat, p = sstats.f_oneway(*(per_group[g] for g in testable))
        result["omnibus"] = {"test": "one-way ANOVA", "statistic": float(f_stat), "p": float(p)}

    pairs = []
    for i, g1 in enumerate(testable):
        for g2 in testable[i + 1 :]:
            t_stat, p = sstats.ttest_ind(per_group[g1], per_group[g2], equal_var=False)
            pairs.append(
                {
                    "groups": (g1, g2),
                    "test": "Welch t-test",
                    "statistic": float(t_stat),
                    "p": float(p),
                    "difference": float(np.mean(per_group[g2]) - np.mean(per_group[g1])),
                }
            )
    result["pairs"] = _holm(pairs)
    return result


def _holm(pairs: list[dict]) -> list[dict]:
    # Holm-Bonferroni: control the family-wise error rate over the pairwise tests.
    order = sorted(range(len(pairs)), key=lambda i: pairs[i]["p"])
    running = 0.0
    for rank, i in enumerate(order):
        adjusted = min(1.0, pairs[i]["p"] * (len(pairs) - rank))
        running = max(running, adjusted)  # keep the adjusted values monotonic
        pairs[i]["p_adjusted"] = running
    return pairs


def stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def print_comparison(result: dict) -> None:
    print(f"\nGroup comparison ({result['unit_name']} as replicates):")
    for g in sorted(result["n"]):
        print(f"  {g}: {result['mean'][g]:.4f} h^-1  (n = {result['n'][g]} {result['unit_name']})")
    for g in result["skipped"]:
        print(f"  WARNING: {g} has fewer than 2 {result['unit_name']}, so it cannot be tested.")
    if result["omnibus"]:
        o = result["omnibus"]
        print(f"  {o['test']}: F = {o['statistic']:.3f}, p = {o['p']:.4g}")
    for pair in result["pairs"]:
        g1, g2 = pair["groups"]
        print(
            f"  {g1} vs {g2}: {pair['test']} t = {pair['statistic']:.3f}, "
            f"p = {pair['p']:.4g}, Holm-adjusted p = {pair['p_adjusted']:.4g} "
            f"{stars(pair['p_adjusted'])}  (difference {pair['difference']:+.4f} h^-1)"
        )


def brackets_to_show(result: dict | None, group_order: list[str]) -> list[dict]:
    """The pairwise tests worth drawing on the figure."""
    if not result:
        return []
    positions = set(group_order)
    shown = [
        p for p in result["pairs"] if p["groups"][0] in positions and p["groups"][1] in positions
    ]
    if len(shown) > 3:
        # Too many brackets to read; keep the ones that survived correction.
        shown = [p for p in shown if p["p_adjusted"] < 0.05]
    return shown


def add_significance_brackets(ax, group_order: list[str], result: dict, y_start: float, step: float):
    """Draw one bracket per tested pair, lowest first. Returns the top y used."""
    positions = {g: i for i, g in enumerate(group_order)}
    # Narrow brackets first, so a wide one never crosses a short one.
    shown = sorted(
        brackets_to_show(result, group_order),
        key=lambda p: abs(positions[p["groups"][1]] - positions[p["groups"][0]]),
    )
    y = y_start
    for pair in shown:
        x1, x2 = positions[pair["groups"][0]], positions[pair["groups"][1]]
        drop = step * 0.22
        ax.plot(
            [x1, x1, x2, x2],
            [y - drop, y, y, y - drop],
            color=TEXT_SECONDARY,
            linewidth=1,
            zorder=4,
            clip_on=False,
        )
        ax.text(
            (x1 + x2) / 2,
            y,
            f"{stars(pair['p_adjusted'])}  p = {pair['p_adjusted']:.3g}",
            ha="center",
            va="bottom",
            fontsize=9,
            color=TEXT_SECONDARY,
        )
        y += step
    return y


def plot_by_group(
    df: pd.DataFrame, groups: pd.DataFrame, out_path: Path, references=None, comparison=None
) -> None:
    merged = df.merge(groups, on=UNIT_COL, how="inner")
    if merged.empty:
        raise ValueError("No reactor in the group file matches the growth-rate data.")

    group_order = sorted(merged["group"].unique())
    styles = group_styles(group_order)
    # Reactors ordered by group, then by name, so each group forms one band.
    unit_order = (
        merged[[UNIT_COL, "group"]]
        .drop_duplicates()
        .sort_values(["group", UNIT_COL])[UNIT_COL]
        .tolist()
    )

    fig, (ax_r, ax_g) = plt.subplots(
        1,
        2,
        figsize=(6 + 0.9 * len(unit_order), 5.4),
        gridspec_kw={"width_ratios": [max(len(unit_order), 3), max(len(group_order) + 1, 2)]},
        sharey=True,
        facecolor=SURFACE,
    )
    style_axes(ax_r)
    style_axes(ax_g)

    unit_group = merged.drop_duplicates(UNIT_COL).set_index(UNIT_COL)["group"]
    values = [merged.loc[merged[UNIT_COL] == u, RATE_COL].to_numpy() for u in unit_order]
    colors = [styles[unit_group[u]]["color"] for u in unit_order]
    _box(ax_r, values, positions=np.arange(len(unit_order)), colors=colors)

    rng = np.random.default_rng(0)  # jitter only, fixed so plots are reproducible
    for i, unit in enumerate(unit_order):
        style = styles[unit_group[unit]]
        y = values[i]
        ax_r.scatter(
            i + rng.uniform(-0.16, 0.16, size=len(y)),
            y,
            s=26,
            color=style["color"],
            marker=style["marker"],
            linewidths=0.9,
            edgecolors=SURFACE,
            alpha=0.95,
            zorder=3,
        )

    # Group bands: a label above each run of reactors, so identity is never colour-alone.
    # Reference lines share the y range, so they must not fall outside it.
    y_top = max([merged[RATE_COL].max(), *(references or ())])
    y_bottom = min([merged[RATE_COL].min(), *(references or ())])
    span = (y_top - y_bottom) or 1.0
    # Significance brackets sit between the data and the per-group averages.
    shown_pairs = brackets_to_show(comparison, group_order)
    bracket_start = y_top + span * 0.06
    bracket_step = span * 0.11
    label_y = bracket_start + bracket_step * len(shown_pairs) + span * 0.03
    for g in group_order:
        members = [i for i, u in enumerate(unit_order) if unit_group[u] == g]
        lo, hi = min(members), max(members)
        ax_r.axvspan(lo - 0.5, hi + 0.5, color=styles[g]["color"], alpha=0.06, zorder=0)
        ax_r.text(
            (lo + hi) / 2,
            label_y,
            f"{g}  (n={len(members)})",
            ha="center",
            va="bottom",
            fontsize=10,
            color=TEXT_SECONDARY,
        )
    ax_r.set_ylim(y_bottom - span * 0.08, label_y + span * 0.1)

    ax_r.set_xticks(np.arange(len(unit_order)))
    ax_r.set_xticklabels(unit_order, fontsize=9, color=TEXT_SECONDARY)
    ax_r.set_xlim(-0.7, len(unit_order) - 0.3)
    ax_r.set_ylabel("Growth rate $\\mu$ (h$^{-1}$)", color=TEXT_SECONDARY, fontsize=10)
    ax_r.set_title("Per reactor", color=TEXT_PRIMARY, fontsize=11, loc="left")
    ax_r.axhline(0, color=TEXT_MUTED, linewidth=1, zorder=1)
    add_reference_lines(ax_r, references, label=False)

    # Right panel: the same rates pooled per group.
    g_values = [merged.loc[merged["group"] == g, RATE_COL].to_numpy() for g in group_order]
    g_colors = [styles[g]["color"] for g in group_order]
    _box(ax_g, g_values, positions=np.arange(len(group_order)), colors=g_colors, width=0.5)
    for i, g in enumerate(group_order):
        y = g_values[i]
        ax_g.scatter(
            i + rng.uniform(-0.13, 0.13, size=len(y)),
            y,
            s=22,
            color=styles[g]["color"],
            marker=styles[g]["marker"],
            linewidths=0.9,
            edgecolors=SURFACE,
            alpha=0.9,
            zorder=3,
        )
        ax_g.text(
            i,
            label_y,
            f"{np.mean(y):.3f} h$^{{-1}}$",
            ha="center",
            va="bottom",
            fontsize=10,
            color=TEXT_SECONDARY,
        )
    ax_g.set_xticks(np.arange(len(group_order)))
    ax_g.set_xticklabels(group_order, fontsize=10, color=TEXT_SECONDARY)
    ax_g.set_xlim(-0.7, len(group_order) - 0.3)
    ax_g.set_title("Pooled per group", color=TEXT_PRIMARY, fontsize=11, loc="left")
    ax_g.axhline(0, color=TEXT_MUTED, linewidth=1, zorder=1)
    add_reference_lines(ax_g, references)
    if shown_pairs:
        add_significance_brackets(ax_g, group_order, comparison, bracket_start, bracket_step)

    if len(group_order) > 1:
        handles = [
            Line2D(
                [],
                [],
                color=styles[g]["color"],
                marker=styles[g]["marker"],
                linestyle="none",
                markersize=8,
                markeredgecolor=SURFACE,
                label=g,
            )
            for g in group_order
        ]
        ax_g.legend(handles=handles, loc="lower left", fontsize=9, frameon=False, labelcolor=TEXT_SECONDARY)

    fig.suptitle(
        "Growth rate by reactor, grouped", color=TEXT_PRIMARY, fontsize=13, x=0.01, ha="left"
    )
    if comparison:
        note = (
            f"Welch t-test on per-{'reactor average' if comparison['level'] == 'reactor' else 'segment'}"
            f" growth rates (n = {', '.join(f'{g}: {n}' for g, n in sorted(comparison['n'].items()))}"
            f" {comparison['unit_name']})"
        )
        if comparison["omnibus"]:
            o = comparison["omnibus"]
            note += f" · {o['test']} p = {o['p']:.3g}"
        if len(comparison["pairs"]) > 1:
            note += " · Holm-adjusted p"
        note += " · ns p≥0.05, * <0.05, ** <0.01, *** <0.001"
        fig.text(0.01, 0.005, note, fontsize=8, color=TEXT_MUTED, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 1) if comparison else None)
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)

    summary = (
        merged.groupby("group")[RATE_COL].agg(["count", "mean", "std"]).round(4).reset_index()
    )
    print("\nPer-group growth rate over all segments (h^-1):")
    print(summary.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--growth-rates",
        type=Path,
        default=Path(__file__).parent / "output" / "growth_rates.csv",
        help="growth_rates.csv written by calculate_growth_rates.py.",
    )
    parser.add_argument(
        "--groups",
        type=Path,
        default=Path(__file__).parent / "groups.csv",
        help="CSV mapping each reactor to a group (two columns: reactor, group).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write the plots (default: alongside growth_rates.csv).",
    )
    parser.add_argument(
        "--min-r2",
        type=float,
        default=0.0,
        help="Drop segments whose fit has a lower R^2 (default: 0.0, keep everything).",
    )
    parser.add_argument(
        "--omit-reactors",
        nargs="+",
        default=[],
        metavar="UNIT",
        help="Reactors to leave out of every plot and out of the group summary, "
        "e.g. --omit-reactors P02. Matched case-insensitively, with a capital O "
        "accepted for a zero.",
    )
    parser.add_argument(
        "--stats",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Test the groups against each other and annotate the grouped plot "
        "with the result (default: on).",
    )
    parser.add_argument(
        "--stats-level",
        choices=("reactor", "segment"),
        default="reactor",
        help="Unit of replication for the group comparison: 'reactor' (default) "
        "uses one average per reactor; 'segment' pools every segment, which "
        "overstates n because segments from one reactor are not independent.",
    )
    parser.add_argument(
        "--average-line",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Draw each reactor's average growth rate as a dashed line on its "
        "scatter plot (default: on; --no-average-line turns it off).",
    )
    parser.add_argument(
        "--reference-line",
        type=float,
        nargs="+",
        metavar="MU",
        help="Draw a dashed horizontal line at this growth rate (h^-1) on every "
        "plot, e.g. a target or a published value. Repeatable: --reference-line 0.6 0.9",
    )
    parser.add_argument(
        "--between-dilutions-only",
        action="store_true",
        help="Keep only segments bounded by dilutions on both sides, dropping the "
        "start-of-experiment (lag phase) and end-of-experiment segments.",
    )
    args = parser.parse_args()

    out_dir = args.output_dir or args.growth_rates.parent

    try:
        df = load_growth_rates(args.growth_rates)
    except (FileNotFoundError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    print(f"Growth rates: {args.growth_rates}  ({len(df)} segments)")

    omitted = set()
    if args.omit_reactors:
        known = {_normalise_unit(u): u for u in df[UNIT_COL].unique()}
        omitted = {known[k] for k in (_normalise_unit(u) for u in args.omit_reactors) if k in known}
        drop = omitted
        unknown = [u for u in args.omit_reactors if _normalise_unit(u) not in known]
        if unknown:
            print(f"WARNING: --omit-reactors: no such reactor in the data: {', '.join(unknown)}")
        if drop:
            n = len(df)
            df = df[~df[UNIT_COL].isin(drop)]
            print(f"Omitted {', '.join(sorted(drop))} ({n - len(df)} segments)")

    if args.between_dilutions_only:
        n = len(df)
        df = df[df["segment_position"] == "between_dilutions"]
        print(f"Kept {len(df)} of {n} segments with segment_position == between_dilutions")
    if args.min_r2 > 0:
        n = len(df)
        df = df[df["r_squared"] >= args.min_r2]
        print(f"Kept {len(df)} of {n} segments with R^2 >= {args.min_r2}")
    if df.empty:
        parser.exit(2, "error: no segments left to plot after filtering.\n")

    out_dir.mkdir(parents=True, exist_ok=True)
    norm = Normalize(vmin=max(0.0, df["r_squared"].min()), vmax=1.0)

    scatter_dir = out_dir / "growth_rate_scatter"
    references = args.reference_line
    try:
        groups = load_groups(args.groups, sorted(df[UNIT_COL].unique()), omitted)
    except (FileNotFoundError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    plot_per_reactor(df, scatter_dir, norm, references, args.average_line)
    plot_small_multiples(
        df, scatter_dir / "all_reactors.png", norm, references, groups, args.average_line
    )
    print(f"Per-reactor scatter plots -> {scatter_dir}")

    group_path = out_dir / "growth_rates_by_group.png"
    comparison = None
    if args.stats:
        merged = df.merge(groups, on=UNIT_COL, how="inner")
        if merged["group"].nunique() < 2:
            print("\nOnly one group in the group file, so there is nothing to compare.")
        else:
            comparison = compare_groups(merged, args.stats_level)
            print_comparison(comparison)
    try:
        plot_by_group(df, groups, group_path, references, comparison)
    except ValueError as exc:
        parser.exit(2, f"error: {exc}\n")
    print(f"\nGrouped plot -> {group_path}")


if __name__ == "__main__":
    main()
