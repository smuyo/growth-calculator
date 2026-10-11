"""
Calculate and plot per-cycle growth rates from a Pioreactor turbidostat export.

One command does the whole analysis. Point it at an export folder and a group
file and it discovers the CSVs it needs, fits a growth rate to every growth
period between dilutions, and writes every figure:

    uv run growth_analysis.py --data-dir input/260720_pio_data \
                              --groups   input/groups.csv

Calculation. For each pioreactor unit, OD readings are split into growth
segments bounded by dilution events. A dilution runs from a "DosingStarted"
event to its matching "DosingStopped" event; every OD reading inside such a
window is discarded, since dosing disturbs the OD signal. Manual (UI-triggered)
pump actions from dosing_events/ are treated as additional exclusion windows by
default, because they perturb OD the same way but are not recorded in
dosing_automation_events/. For every remaining segment, ln(OD) is regressed
against time (hours_since_experiment_created) to obtain the growth rate (h^-1)
and the R^2 of the fit.

Plots. Growth rates are shaded by the R^2 of their fit, so a poor regression is
visible rather than silently averaged in. Use --min-r2 / --between-dilutions-only
to drop segments that are not exponential growth (see the README), and
--omit-reactors to leave a reactor out entirely. --reference-line draws a dashed
horizontal rule at a given growth rate on every plot, for comparing against a
target or a published value.

Outputs (written to --output-dir, by default ./YYMMDD_growth_analysis):
  - growth_rates.csv                    one row per growth segment
  - od_calibration.csv                  the OD calibration applied (if any)
  - segment_fit_plots/<unit>/*.png      ln(OD) vs time with the fitted line,
                                        one plot per segment
  - od_plots/<unit>_od_timeline.png     OD vs time per reactor, exclusion
                                        windows shaded
  - od_plots/all_reactors.png           the OD traces as small multiples
  - growth_rate_plots/<unit>.png        growth rate vs time, one per reactor
  - growth_rate_plots/all_reactors.png  the same as small multiples
  - od_growth_rate_plots/<unit>.png     OD with the growth rates overlaid on
                                        a second axis, one per reactor
  - od_growth_rate_plots/all_reactors.png
                                        the same as small multiples
  - growth_rates_by_group.png           every reactor side by side, ordered and
                                        coloured by group, plus the per-group
                                        summary and the group comparison

With temperature readings in the export (or --temperatures) it also writes:
  - temperature_plots/<unit>.png        temperature vs time, one per reactor
  - temperature_plots/all_reactors.png  the same as small multiples
  - temperature_growth_rate_plots/<unit>.png
                                        temperature with the growth rates
                                        overlaid on a second axis
"""

import argparse
import math
import shutil
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from scipy import stats
from scipy import stats as sstats

# Every PNG is written at print resolution, so figures can go straight onto a
# poster or slide without looking soft. The per-segment regression plots are
# the exception: there are hundreds of them and they are for checking fits, so
# they default to screen resolution (--high-res-regression-plots overrides).
FIGURE_DPI = 300
REGRESSION_DPI = 150

# =============================================================================
# Growth-rate calculation
# =============================================================================

TIME_COL = "hours_since_experiment_created"

# Manual UI pump actions are logged as many individual tiny-volume rows. Rows
# closer together than this are collapsed into one exclusion window.
MANUAL_BURST_GAP_H = 0.05

OD_COLUMNS = ["pioreactor_unit", "od_reading", TIME_COL]
AUTOMATION_COLUMNS = ["pioreactor_unit", "event_name", TIME_COL]
MANUAL_COLUMNS = ["pioreactor_unit", "source_of_event", TIME_COL]


def resolve_input_file(
    explicit: Path | None,
    data_dir: Path | None,
    subdir: str,
    flag: str,
    required: list[str] | None = None,
) -> Path:
    # Pick the CSV to read: an explicit --*-file wins, else the one CSV in data_dir/subdir.
    # With several CSVs there, ones lacking the `required` columns are set aside, so
    # an unrelated file kept alongside (e.g. the OD calibration file) does no harm.
    if explicit is not None:
        if explicit.is_dir():
            raise IsADirectoryError(f"{flag} expects a CSV file, but got the directory {explicit}")
        if not explicit.is_file():
            raise FileNotFoundError(f"{flag}: no such file: {explicit}")
        return explicit

    if data_dir is None:
        raise FileNotFoundError(f"Provide either {flag} or --data-dir to locate the {subdir} CSV.")

    directory = data_dir / subdir
    if not directory.is_dir():
        raise FileNotFoundError(
            f"{directory} does not exist. Pass {flag} to point directly at the CSV."
        )
    files = sorted(directory.glob("*.csv"))
    if len(files) > 1 and required:
        matching = [f for f in files if set(required) <= set(pd.read_csv(f, nrows=0).columns)]
        if len(matching) == 1:
            return matching[0]
    if not files:
        raise FileNotFoundError(
            f"No CSV file found in {directory}. Pass {flag} to point directly at the CSV."
        )
    if len(files) > 1:
        names = "\n  ".join(f.name for f in files)
        raise ValueError(
            f"{directory} contains {len(files)} CSV files, so the input is ambiguous. "
            f"Pass {flag} to choose one:\n  {names}"
        )
    return files[0]


def _read_csv(path: Path, required: list[str], label: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, parse_dates=["timestamp"])
    except ValueError:
        # A file without a 'timestamp' column is still usable; TIME_COL is what we fit on.
        df = pd.read_csv(path)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} file {path} is missing required column(s): {', '.join(missing)}.\n"
            f"Found columns: {', '.join(map(str, df.columns))}"
        )
    return df


def load_od_readings(path: Path) -> pd.DataFrame:
    df = _read_csv(path, OD_COLUMNS, "OD readings")
    return df.sort_values(["pioreactor_unit", TIME_COL]).reset_index(drop=True)


def load_dosing_automation_events(path: Path) -> pd.DataFrame:
    df = _read_csv(path, AUTOMATION_COLUMNS, "Dosing automation events")
    return df.sort_values(["pioreactor_unit", TIME_COL]).reset_index(drop=True)


def load_manual_dosing_events(path: Path) -> pd.DataFrame:
    df = _read_csv(path, MANUAL_COLUMNS, "Dosing events")
    manual = df[df["source_of_event"] == "UI"]
    return manual.sort_values(["pioreactor_unit", TIME_COL]).reset_index(drop=True)


def load_calibration(path: Path) -> pd.DataFrame:
    """Read the OD calibration file: two true OD values per reactor.

    The smaller value is the true OD at the start of the run, the larger the
    true OD at its end; their order in the file does not matter. The OD column
    is whichever one is numeric, so the header names are not relied on.
    """
    raw = pd.read_csv(path)
    if raw.shape[1] < 2:
        raise ValueError(f"{path} needs two columns (reactor and OD).")
    numeric = [c for c in raw.columns if pd.to_numeric(raw[c], errors="coerce").notna().all()]
    if not numeric or len(numeric) == raw.shape[1]:
        raise ValueError(f"{path}: cannot tell the reactor column from the OD column.")
    od_col = numeric[0]
    reactor_col = next(c for c in raw.columns if c not in numeric)

    cal = pd.DataFrame(
        {"key": raw[reactor_col].map(_normalise_unit), "od": pd.to_numeric(raw[od_col])}
    )
    counts = cal.groupby("key")["od"].count()
    wrong = counts[counts != 2]
    if not wrong.empty:
        listed = ", ".join(f"{k} ({n})" for k, n in wrong.items())
        raise ValueError(f"{path}: every reactor needs exactly 2 OD values; got {listed}.")
    out = cal.groupby("key")["od"].agg(true_low="min", true_high="max")
    flat = out.index[out["true_low"] == out["true_high"]]
    if len(flat):
        raise ValueError(f"{path}: low and high OD are equal for {', '.join(flat)}.")
    return out


def calibrate_od(
    od: pd.DataFrame, calibration: pd.DataFrame, skip: int, points: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Linearly rescale each reactor's OD readings onto its calibration values.

    The measured start OD is the median of the `points` readings after the first
    `skip` (the first few readings of a run read high); the measured end OD is
    the median of the last `points` readings, so a few outliers at the end do
    not move it. Those two are mapped onto the true low and high OD.
    """
    od = od.copy()
    rows, missing = [], []
    for unit, unit_od in od.groupby("pioreactor_unit"):
        key = _normalise_unit(unit)
        if key not in calibration.index:
            missing.append(unit)
            continue
        values = unit_od.sort_values(TIME_COL)["od_reading"].to_numpy(dtype=float)
        if len(values) < skip + 2 * points:
            raise ValueError(
                f"{unit} has only {len(values)} OD readings, too few to calibrate with "
                f"--calibration-skip {skip} and --calibration-points {points}."
            )
        measured_low = float(np.median(values[skip : skip + points]))
        measured_high = float(np.median(values[-points:]))
        true_low, true_high = calibration.loc[key, ["true_low", "true_high"]]
        if measured_high <= measured_low:
            raise ValueError(
                f"{unit}: the measured end OD ({measured_high:.4g}) is not above the "
                f"measured start OD ({measured_low:.4g}), so it cannot be calibrated."
            )
        slope = (true_high - true_low) / (measured_high - measured_low)
        intercept = true_low - slope * measured_low
        od.loc[unit_od.index, "od_reading"] = slope * unit_od["od_reading"] + intercept
        rows.append(
            {
                "pioreactor_unit": unit,
                "measured_low": measured_low,
                "measured_high": measured_high,
                "true_low": true_low,
                "true_high": true_high,
                "slope": slope,
                "intercept": intercept,
            }
        )
    if missing:
        raise ValueError(
            f"no calibration values for {', '.join(missing)}. Add them to the calibration "
            "file, or leave those reactors out with --exclude-reactors."
        )
    unused = sorted(set(calibration.index) - {_normalise_unit(r["pioreactor_unit"]) for r in rows})
    if unused:
        print(f"NOTE: calibration file reactors with no OD readings: {', '.join(unused)}")
    return od, pd.DataFrame(rows)


def automation_dilution_windows(events: pd.DataFrame, unit: str) -> list[dict]:
    """Pair each DosingStarted with the following DosingStopped for one unit.

    A trailing DosingStarted with no DosingStopped (dilution still running when
    the export was taken) yields a window with end=None, meaning "discard
    everything from here to the end of the record".
    """
    unit_events = events[events["pioreactor_unit"] == unit].sort_values(TIME_COL)
    windows = []
    pending_start = None
    for time_h, event_name in zip(unit_events[TIME_COL], unit_events["event_name"]):
        if event_name == "DosingStarted":
            if pending_start is not None:
                # Two starts with no stop between them: close the first at the
                # second's start so its readings are still excluded.
                windows.append({"start": pending_start, "end": time_h, "kind": "automation"})
            pending_start = time_h
        elif event_name == "DosingStopped" and pending_start is not None:
            windows.append({"start": pending_start, "end": time_h, "kind": "automation"})
            pending_start = None
    if pending_start is not None:
        windows.append({"start": pending_start, "end": None, "kind": "automation"})
    return windows


def manual_dosing_windows(manual: pd.DataFrame, unit: str) -> list[dict]:
    # Collapse bursts of manual UI pump rows into exclusion windows
    times = manual.loc[manual["pioreactor_unit"] == unit, TIME_COL].sort_values().to_numpy()
    if times.size == 0:
        return []
    split_at = np.flatnonzero(np.diff(times) > MANUAL_BURST_GAP_H)
    starts = np.concatenate([times[:1], times[split_at + 1]])
    ends = np.concatenate([times[split_at], times[-1:]])
    return [{"start": s, "end": e, "kind": "manual"} for s, e in zip(starts, ends)]


def merge_windows(windows: list[dict]) -> list[dict]:
    """Sort and merge overlapping exclusion windows.

    An open-ended window (end=None) absorbs everything after its start.
    """
    if not windows:
        return []
    ordered = sorted(windows, key=lambda w: w["start"])
    merged = [dict(ordered[0])]
    for win in ordered[1:]:
        last = merged[-1]
        if last["end"] is None:
            continue
        if win["start"] <= last["end"]:
            if win["end"] is None:
                last["end"] = None
            else:
                last["end"] = max(last["end"], win["end"])
            if win["kind"] != last["kind"]:
                last["kind"] = "mixed"
        else:
            merged.append(dict(win))
    return merged


def excluded_mask(times: np.ndarray, windows: list[dict]) -> np.ndarray:
    """Flag readings taken while a pump was running.

    Both window bounds are inclusive: a reading stamped exactly at DosingStarted
    or DosingStopped is concurrent with the dosing and must be dropped.
    """
    mask = np.zeros(times.shape, dtype=bool)
    for win in windows:
        if win["end"] is None:
            mask |= times >= win["start"]
        else:
            mask |= (times >= win["start"]) & (times <= win["end"])
    return mask


def build_segments(od_unit: pd.DataFrame, windows: list[dict], min_points: int) -> list[dict]:
    # Split one unit's OD readings into the maximal runs that no pump event touches
    times = od_unit[TIME_COL].to_numpy(dtype=float)
    keep = ~excluded_mask(times, windows)

    # Each contiguous run of kept readings is one growth segment.
    run_id = np.cumsum(~keep)
    segments = []
    for _, idx in pd.Series(np.arange(len(keep))[keep]).groupby(run_id[keep]):
        block = od_unit.iloc[idx.to_numpy()]
        if len(block) >= min_points:
            segments.append({"data": block})

    n = len(segments)
    for i, seg in enumerate(segments):
        seg["segment_index"] = i
        if n == 1:
            seg["segment_position"] = "only_segment"
        elif i == 0:
            seg["segment_position"] = "experiment_start"
        elif i == n - 1:
            seg["segment_position"] = "experiment_end"
        else:
            seg["segment_position"] = "between_dilutions"
    return segments


def fit_growth_rate(seg_data: pd.DataFrame):
    t = seg_data[TIME_COL].to_numpy(dtype=float)
    ln_od = np.log(seg_data["od_reading"].to_numpy(dtype=float))
    result = stats.linregress(t, ln_od)
    fit = {
        "start_time_h": t.min(),
        "end_time_h": t.max(),
        "duration_h": t.max() - t.min(),
        "growth_rate_h-1": result.slope,
        "doubling_time_h": np.log(2) / result.slope if result.slope > 0 else np.nan,
        "r_squared": result.rvalue**2,
        "p_value": result.pvalue,
        "std_err_h-1": result.stderr,
        "intercept_ln_od": result.intercept,
        "od_start": seg_data["od_reading"].iloc[0],
        "od_end": seg_data["od_reading"].iloc[-1],
        "n_points": len(t),
    }
    return fit, t, ln_od, result


def plot_regression(unit, fit, t, ln_od, result, out_dir: Path, dpi: int = REGRESSION_DPI):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(t, ln_od, s=18, color="#1f77b4", label="ln(OD)", zorder=3)
    t_fit = np.array([t.min(), t.max()])
    ax.plot(
        t_fit,
        result.intercept + result.slope * t_fit,
        color="#d62728",
        linewidth=2,
        label="Linear fit",
        zorder=4,
    )
    ax.set_xlabel("Time since experiment start (h)")
    ax.set_ylabel("ln(OD)")
    ax.set_title(
        f"{unit} — segment {fit['segment_index']} ({fit['segment_position'].replace('_', ' ')})"
    )
    doubling = fit["doubling_time_h"]
    doubling_txt = f"{doubling:.2f} h" if np.isfinite(doubling) else "n/a"
    textstr = (
        f"$\\mu$ = {fit['growth_rate_h-1']:.4f} $\\pm$ {fit['std_err_h-1']:.4f} h$^{{-1}}$\n"
        f"$R^2$ = {fit['r_squared']:.4f}\n"
        f"$t_d$ = {doubling_txt}\n"
        f"n = {fit['n_points']}, {fit['duration_h']:.2f} h"
    )
    ax.text(
        0.03,
        0.97,
        textstr,
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#bbbbbb"),
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{unit}_segment_{fit['segment_index']:03d}.png", dpi=dpi)
    plt.close(fig)


def plot_reactor_timeline(unit, od_unit, windows, unit_results, out_dir: Path):
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.plot(
        od_unit[TIME_COL],
        od_unit["od_reading"],
        color="#bbbbbb",
        linewidth=0.7,
        zorder=1,
        label="OD (excluded)",
    )

    for i, row in enumerate(unit_results):
        seg = od_unit[
            (od_unit[TIME_COL] >= row["start_time_h"]) & (od_unit[TIME_COL] <= row["end_time_h"])
        ]
        ax.plot(
            seg[TIME_COL],
            seg["od_reading"],
            color="#1f77b4",
            linewidth=1.2,
            zorder=3,
            label="OD (used for fits)" if i == 0 else None,
        )

    label_seen = set()
    colors = {"automation": "#d62728", "manual": "#ff7f0e", "mixed": "#9467bd"}
    labels = {
        "automation": "Automated dilution (removed)",
        "manual": "Manual dosing (removed)",
        "mixed": "Dilution + manual (removed)",
    }
    t_max = od_unit[TIME_COL].max()
    for win in windows:
        end = win["end"] if win["end"] is not None else t_max
        kind = win["kind"]
        ax.axvspan(
            win["start"],
            end,
            color=colors[kind],
            alpha=0.35,
            zorder=2,
            label=labels[kind] if kind not in label_seen else None,
        )
        label_seen.add(kind)

    ax.set_xlabel("Time since experiment start (h)")
    ax.set_ylabel("OD reading")
    ax.set_title(f"{unit} — OD over time, split at dilutions ({len(unit_results)} segments fitted)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{unit}_od_timeline.png", dpi=FIGURE_DPI)
    plt.close(fig)


# =============================================================================
# Plotting
# =============================================================================

RATE_COL = "growth_rate_h-1"
UNIT_COL = "pioreactor_unit"
TEMP_COL = "temperature_c"
TEMP_TIME_COL = "hours_since_experiment_created"

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
# Temperature reads as a physical setting, so it gets its own warm hue and is
# never one of the categorical group slots.
TEMP_COLOR = "#b3261e"

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


def prepare_growth_rates(df: pd.DataFrame, source: str = "the growth-rate table") -> pd.DataFrame:
    """Check the fitted segments and add the column the plots are drawn against."""
    required = [UNIT_COL, RATE_COL, "start_time_h", "end_time_h", "r_squared", "segment_position"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{source} is missing required column(s): {', '.join(missing)}.\n"
            f"Found columns: {', '.join(map(str, df.columns))}"
        )
    df = df.copy()
    # Each segment spans an interval; plot it at its midpoint.
    df["time_h"] = (df["start_time_h"] + df["end_time_h"]) / 2.0
    return df


def load_growth_rates(path: Path) -> pd.DataFrame:
    return prepare_growth_rates(pd.read_csv(path), str(path))


def _normalise_unit(name: str) -> str:
    # Reactor ids are letter+digits; a typed capital O for a zero is a common slip.
    return str(name).strip().upper().replace("O", "0")


def drop_excluded(df: pd.DataFrame, excluded: list[str], found: set[str]) -> pd.DataFrame:
    """Remove the --exclude-reactors units from one input table.

    Every excluded name this table actually held is added to `found`, so names
    matching none of the inputs can be reported once at the end.
    """
    if not excluded:
        return df
    keys = df[UNIT_COL].map(_normalise_unit)
    hit = keys.isin({_normalise_unit(u) for u in excluded})
    found.update(keys[hit])
    return df[~hit]


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
        fig.savefig(out_dir / f"{unit}.png", dpi=FIGURE_DPI, facecolor=SURFACE)
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


def _multiples_figure(units: list[str], groups: pd.DataFrame | None, panel_width: float = 4.6):
    """Blank small-multiples grid: one panel per reactor, one row per group."""
    rows, row_labels = _grid_layout(units, groups)
    nrows = len(rows)
    ncols = max(len(r) for r in rows)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(panel_width * ncols, 3.6 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
        facecolor=SURFACE,
    )
    return fig, axes, rows, row_labels


def _label_multiples(fig, axes, rows, row_labels, xlabel: str, ylabel: str) -> None:
    """Row headings and axis labels for a grid from _multiples_figure()."""
    styles = group_styles(row_labels) if row_labels else {}
    if row_labels:
        for r, row_units in enumerate(rows):
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
    nrows, ncols = len(rows), max(len(r) for r in rows)
    for col in range(ncols):
        visible = [axes[row][col] for row in range(nrows) if axes[row][col].get_visible()]
        if visible:
            visible[-1].set_xlabel(xlabel, color=TEXT_SECONDARY, fontsize=9)
            visible[-1].tick_params(labelbottom=True)
    for row in axes:
        row[0].set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=9)

    # Row headings sit above each row's panels, so the rows need room between them.
    fig.subplots_adjust(
        top=0.90 - 0.02 * nrows if row_labels else 0.93,
        hspace=0.42 if row_labels else 0.22,
        wspace=0.16,
    )


def _multiples_header(fig, axes, title: str, handles: list, row_labels) -> None:
    """Title and legend in one band just above the grid, spanning only the panels.

    Both are pinned to the panels' outer edges rather than the figure's, so the
    tight bounding box at save time trims the side margins away. Call this last,
    after any colorbar has shifted the panels.
    """
    visible = [ax for row in axes for ax in row if ax.get_visible()]
    left = min(ax.get_position().x0 for ax in visible)
    right = max(ax.get_position().x1 for ax in visible)
    top = max(ax.get_position().y1 for ax in visible)
    # Clear the row heading and panel title above the first row (inches -> figure fraction).
    y = top + (0.75 if row_labels else 0.45) / fig.get_figheight()
    fig.text(left, y, title, color=TEXT_PRIMARY, fontsize=13, ha="left", va="bottom")
    if handles:
        fig.legend(
            handles=handles,
            loc="lower right",
            bbox_to_anchor=(right, y),
            ncol=len(handles),
            fontsize=9,
            frameon=False,
            labelcolor=TEXT_SECONDARY,
            borderaxespad=0,
        )


def _r2_colorbar(fig, ax_or_axes, sc, **kwargs):
    cbar = fig.colorbar(sc, ax=ax_or_axes, **kwargs)
    cbar.set_label("$R^2$ of fit", color=TEXT_SECONDARY, fontsize=9)
    cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    cbar.outline.set_visible(False)
    return cbar


def plot_small_multiples(df, out_path: Path, norm, references=None, groups=None, show_average=True) -> None:
    units = sorted(df[UNIT_COL].unique())
    fig, axes, rows, row_labels = _multiples_figure(units, groups)

    sc = None
    for r, row_units in enumerate(rows):
        for c in range(len(axes[r])):
            ax = axes[r][c]
            if c >= len(row_units):
                ax.set_visible(False)
                continue
            unit = row_units[c]
            sc, _ = scatter_one_reactor(ax, df[df[UNIT_COL] == unit], norm, references, show_average)
            ax.set_title(unit, color=TEXT_PRIMARY, fontsize=11, loc="left")

    _label_multiples(
        fig,
        axes,
        rows,
        row_labels,
        "Time (h)",
        "$\\mu$ (h$^{-1}$)",
    )
    if sc is not None:
        _r2_colorbar(fig, axes, sc, pad=0.015, fraction=0.02)
    # The panels carry no legend of their own, so name the dashed rules here.
    _multiples_header(
        fig,
        axes,
        "Growth rate per segment over time, by reactor" + (", grouped" if row_labels else ""),
        rule_handles(references, show_average),
        row_labels,
    )
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


# --- OD ----------------------------------------------------------------------
OD_USED_COLOR = "#2a78d6"
OD_EXCLUDED_COLOR = "#c9c8c3"


def plot_od_small_multiples(
    od_trace: pd.DataFrame, out_path: Path, groups=None, calibrated: bool = False
) -> None:
    """OD against time for every reactor, laid out like the growth-rate grid.

    `od_trace` has one row per OD reading and a boolean `used` column: readings
    inside a fitted segment are drawn in blue, the rest (dilutions, manual
    dosing, unfitted stretches) in grey underneath.
    """
    units = sorted(od_trace[UNIT_COL].unique())
    # Wider panels than the other grids: a run's many dilution cycles need room.
    fig, axes, rows, row_labels = _multiples_figure(units, groups, panel_width=7.0)

    for r, row_units in enumerate(rows):
        for c in range(len(axes[r])):
            ax = axes[r][c]
            if c >= len(row_units):
                ax.set_visible(False)
                continue
            unit = row_units[c]
            unit_od = od_trace[od_trace[UNIT_COL] == unit]
            style_axes(ax)
            ax.plot(
                unit_od[TIME_COL],
                unit_od["od_reading"],
                color=OD_EXCLUDED_COLOR,
                linewidth=0.9,
                zorder=2,
            )
            # Masking the unused readings breaks the line at every dilution.
            ax.plot(
                unit_od[TIME_COL],
                unit_od["od_reading"].where(unit_od["used"]),
                color=OD_USED_COLOR,
                linewidth=1.3,
                zorder=3,
            )
            ax.set_title(unit, color=TEXT_PRIMARY, fontsize=11, loc="left")

    _label_multiples(
        fig,
        axes,
        rows,
        row_labels,
        "Time (h)",
        "OD (calibrated)" if calibrated else "OD reading",
    )
    _multiples_header(
        fig,
        axes,
        "OD over time, by reactor" + (", grouped" if row_labels else ""),
        [
            Line2D([], [], color=OD_USED_COLOR, linewidth=1.8, label="used for fits"),
            Line2D([], [], color=OD_EXCLUDED_COLOR, linewidth=1.8, label="excluded"),
        ],
        row_labels,
    )
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


# --- temperature -------------------------------------------------------------
def load_temperatures(path: Path) -> pd.DataFrame:
    """Read the temperature export: one row per reading per reactor."""
    df = pd.read_csv(path)
    required = [UNIT_COL, TEMP_COL, TEMP_TIME_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}.\n"
            f"Found columns: {', '.join(map(str, df.columns))}"
        )
    df = df[required].dropna().copy()
    df["time_h"] = df[TEMP_TIME_COL].astype(float)
    # Readings are exported newest-first per reactor; a line plot needs them in order.
    return df.sort_values([UNIT_COL, "time_h"]).reset_index(drop=True)


def temperature_legend_handles(show_average: bool = True) -> list[Line2D]:
    handles = [Line2D([], [], color=TEMP_COLOR, linewidth=1.8, label="temperature")]
    if show_average:
        handles.append(
            Line2D([], [], color="#eb6834", linestyle="--", linewidth=2, label="reactor average")
        )
    return handles


def line_one_reactor(ax, unit_df: pd.DataFrame, show_average: bool = True) -> float:
    """Temperature against time for one reactor. Returns the reactor mean."""
    style_axes(ax)
    ax.plot(
        unit_df["time_h"],
        unit_df[TEMP_COL],
        color=TEMP_COLOR,
        linewidth=1.4,
        solid_capstyle="round",
        zorder=3,
    )
    mean = unit_df[TEMP_COL].mean()
    if show_average:
        ax.axhline(mean, color="#eb6834", linewidth=2, linestyle="--", zorder=2)
    return mean


def plot_temperature_per_reactor(temps: pd.DataFrame, out_dir: Path, show_average=True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for unit, unit_df in temps.groupby(UNIT_COL, sort=True):
        fig, ax = plt.subplots(figsize=(8, 4.6), facecolor=SURFACE)
        mean = line_one_reactor(ax, unit_df, show_average)
        ax.set_xlabel("Time since experiment start (h)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_ylabel("Temperature (\u00b0C)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_title(
            f"{unit} \u2014 temperature over time  "
            f"(n = {len(unit_df)}, average {mean:.2f} \u00b0C)",
            color=TEXT_PRIMARY,
            fontsize=12,
            loc="left",
        )
        ax.legend(
            handles=temperature_legend_handles(show_average),
            loc="best",
            fontsize=9,
            frameon=False,
            labelcolor=TEXT_SECONDARY,
        )
        fig.tight_layout()
        fig.savefig(out_dir / f"{unit}.png", dpi=FIGURE_DPI, facecolor=SURFACE)
        plt.close(fig)


def plot_temperature_small_multiples(
    temps: pd.DataFrame, out_path: Path, groups=None, show_average=True
) -> None:
    units = sorted(temps[UNIT_COL].unique())
    fig, axes, rows, row_labels = _multiples_figure(units, groups)

    for r, row_units in enumerate(rows):
        for c in range(len(axes[r])):
            ax = axes[r][c]
            if c >= len(row_units):
                ax.set_visible(False)
                continue
            unit = row_units[c]
            line_one_reactor(ax, temps[temps[UNIT_COL] == unit], show_average)
            ax.set_title(unit, color=TEXT_PRIMARY, fontsize=11, loc="left")

    _label_multiples(
        fig,
        axes,
        rows,
        row_labels,
        "Time (h)",
        "Temperature (\u00b0C)",
    )
    _multiples_header(
        fig,
        axes,
        "Temperature over time, by reactor" + (", grouped" if row_labels else ""),
        temperature_legend_handles(show_average),
        row_labels,
    )
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def _rate_twin_axis(ax, unit_rates: pd.DataFrame, norm, references=None, show_average=False):
    """Right-hand axis on `ax` with one reactor's growth rates, shaded by R^2.

    Returns the new axis and the scatter (for the colorbar).
    """
    ax2 = ax.twinx()
    ax2.set_facecolor("none")
    for side in ("top", "left"):
        ax2.spines[side].set_visible(False)
    ax2.spines["right"].set_color(GRID)
    ax2.spines["bottom"].set_color(GRID)
    ax2.tick_params(colors=TEXT_SECONDARY, labelsize=9, length=3, color=GRID)
    sc = ax2.scatter(
        unit_rates["time_h"],
        unit_rates[RATE_COL],
        c=unit_rates["r_squared"],
        cmap=R2_CMAP,
        norm=norm,
        s=46,
        linewidths=1.2,
        edgecolors=SURFACE,
        zorder=4,
    )
    if show_average:
        ax2.axhline(
            unit_rates[RATE_COL].mean(), color="#eb6834", linewidth=2, linestyle="--", zorder=3
        )
    add_reference_lines(ax2, references)
    ax2.set_ylabel("Growth rate $\\mu$ (h$^{-1}$)", color=TEXT_SECONDARY, fontsize=10)
    return ax2, sc


def _growth_rate_handle() -> Line2D:
    return Line2D(
        [],
        [],
        color=R2_CMAP(0.9),
        marker="o",
        linestyle="none",
        markersize=6,
        label="growth rate (right)",
    )


def plot_temperature_with_rates(
    temps: pd.DataFrame,
    rates: pd.DataFrame,
    out_dir: Path,
    norm,
    references=None,
    show_average=True,
) -> list[str]:
    """One plot per reactor: temperature as a line, growth rates as points.

    The two quantities share the time axis but nothing else, so \u00b5 goes on a
    right-hand axis. Only reactors present in both inputs are plotted; the
    returned list names them.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    plotted = []
    for unit in sorted(set(temps[UNIT_COL]) & set(rates[UNIT_COL])):
        unit_temps = temps[temps[UNIT_COL] == unit]
        unit_rates = rates[rates[UNIT_COL] == unit]
        # Wider than the other per-reactor figures: two y axes plus a colorbar.
        fig, ax = plt.subplots(figsize=(9.6, 4.8), facecolor=SURFACE)
        mean_t = line_one_reactor(ax, unit_temps, show_average)
        ax.set_xlabel("Time since experiment start (h)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_ylabel("Temperature (\u00b0C)", color=TEMP_COLOR, fontsize=10)
        ax.tick_params(axis="y", colors=TEMP_COLOR)

        ax2, sc = _rate_twin_axis(ax, unit_rates, norm, references)

        mean_mu = unit_rates[RATE_COL].mean()
        ax.set_title(
            f"{unit} \u2014 temperature and growth rate over time\n"
            f"average {mean_t:.2f} \u00b0C, {mean_mu:.3f} h$^{{-1}}$ over {len(unit_rates)} segments",
            color=TEXT_PRIMARY,
            fontsize=11,
            loc="left",
        )
        _r2_colorbar(fig, ax2, sc, pad=0.12)
        handles = [
            Line2D([], [], color=TEMP_COLOR, linewidth=1.8, label="temperature (left)"),
            _growth_rate_handle(),
        ]
        if show_average:
            handles.append(
                Line2D(
                    [], [], color="#eb6834", linestyle="--", linewidth=2, label="average temperature"
                )
            )
        if references:
            handles.extend(h for h in rule_handles(references, show_average=False))
        ax.legend(
            handles=handles, loc="best", fontsize=9, frameon=False, labelcolor=TEXT_SECONDARY
        )
        fig.tight_layout()
        fig.savefig(out_dir / f"{unit}.png", dpi=FIGURE_DPI, facecolor=SURFACE)
        plt.close(fig)
        plotted.append(unit)
    return plotted


OD_TRACE_COLOR = "#52514e"  # neutral, so it never competes with the blue R^2 ramp


def plot_od_with_rates(
    od_trace: pd.DataFrame,
    rates: pd.DataFrame,
    out_dir: Path,
    norm,
    references=None,
    show_average=True,
    calibrated: bool = False,
) -> list[str]:
    """One plot per reactor: the OD trace as a line, growth rates as points.

    Laid out like plot_temperature_with_rates(): OD on the left axis, µ on
    the right. Readings used in fits are drawn dark, the rest light grey.
    Returns the reactors plotted.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    od_label = "OD (calibrated)" if calibrated else "OD reading"
    plotted = []
    for unit in sorted(set(od_trace[UNIT_COL]) & set(rates[UNIT_COL])):
        unit_od = od_trace[od_trace[UNIT_COL] == unit]
        unit_rates = rates[rates[UNIT_COL] == unit]
        fig, ax = plt.subplots(figsize=(9.6, 4.8), facecolor=SURFACE)
        style_axes(ax)
        ax.plot(
            unit_od[TIME_COL],
            unit_od["od_reading"],
            color=OD_EXCLUDED_COLOR,
            linewidth=0.9,
            zorder=2,
        )
        # Masking the unused readings breaks the line at every dilution.
        ax.plot(
            unit_od[TIME_COL],
            unit_od["od_reading"].where(unit_od["used"]),
            color=OD_TRACE_COLOR,
            linewidth=1.1,
            zorder=3,
        )
        ax.set_xlabel("Time since experiment start (h)", color=TEXT_SECONDARY, fontsize=10)
        ax.set_ylabel(od_label, color=OD_TRACE_COLOR, fontsize=10)

        ax2, sc = _rate_twin_axis(ax, unit_rates, norm, references, show_average)
        mean_mu = unit_rates[RATE_COL].mean()
        ax.set_title(
            f"{unit} — OD and growth rate over time\n"
            f"average {mean_mu:.3f} h$^{{-1}}$ over {len(unit_rates)} segments",
            color=TEXT_PRIMARY,
            fontsize=11,
            loc="left",
        )
        _r2_colorbar(fig, ax2, sc, pad=0.12)
        handles = [
            Line2D([], [], color=OD_TRACE_COLOR, linewidth=1.8, label="OD used for fits (left)"),
            Line2D([], [], color=OD_EXCLUDED_COLOR, linewidth=1.8, label="OD excluded (left)"),
            _growth_rate_handle(),
        ]
        handles.extend(rule_handles(references, show_average))
        ax.legend(
            handles=handles, loc="best", fontsize=9, frameon=False, labelcolor=TEXT_SECONDARY
        )
        fig.tight_layout()
        fig.savefig(out_dir / f"{unit}.png", dpi=FIGURE_DPI, facecolor=SURFACE)
        plt.close(fig)
        plotted.append(unit)
    return plotted


def plot_od_rates_small_multiples(
    od_trace: pd.DataFrame,
    rates: pd.DataFrame,
    out_path: Path,
    norm,
    groups=None,
    references=None,
    show_average=True,
    calibrated: bool = False,
) -> None:
    """The plot_od_with_rates() panels as one grid, laid out like the other grids.

    Both y axes are shared across panels: OD on the left via the grid itself,
    growth rate on the right by sharing every twin axis with the first one.
    """
    units = sorted(set(od_trace[UNIT_COL]) & set(rates[UNIT_COL]))
    fig, axes, rows, row_labels = _multiples_figure(units, groups, panel_width=7.0)

    sc = None
    first_twin = None
    for r, row_units in enumerate(rows):
        for c in range(len(axes[r])):
            ax = axes[r][c]
            if c >= len(row_units):
                ax.set_visible(False)
                continue
            unit = row_units[c]
            unit_od = od_trace[od_trace[UNIT_COL] == unit]
            style_axes(ax)
            ax.plot(
                unit_od[TIME_COL],
                unit_od["od_reading"],
                color=OD_EXCLUDED_COLOR,
                linewidth=0.9,
                zorder=2,
            )
            ax.plot(
                unit_od[TIME_COL],
                unit_od["od_reading"].where(unit_od["used"]),
                color=OD_TRACE_COLOR,
                linewidth=1.1,
                zorder=3,
            )
            ax.set_title(unit, color=TEXT_PRIMARY, fontsize=11, loc="left")

            ax2, sc = _rate_twin_axis(
                ax, rates[rates[UNIT_COL] == unit], norm, references, show_average
            )
            if first_twin is None:
                first_twin = ax2
            else:
                ax2.sharey(first_twin)
            # Growth-rate ticks and label only on the last panel of each row.
            if c != len(row_units) - 1:
                ax2.tick_params(labelright=False)
                ax2.set_ylabel("")
            else:
                ax2.set_ylabel("$\\mu$ (h$^{-1}$)", color=TEXT_SECONDARY, fontsize=9)

    if first_twin is not None:
        # Shared autoscaling leaves the extreme points touching the frame, so set
        # the range explicitly, with room for the points and any reference line.
        values = np.concatenate(
            [rates.loc[rates[UNIT_COL].isin(units), RATE_COL], references or []]
        )
        lo, hi = float(values.min()), float(values.max())
        pad = 0.08 * (hi - lo or 1.0)
        first_twin.set_ylim(lo - pad, hi + pad)

    _label_multiples(
        fig,
        axes,
        rows,
        row_labels,
        "Time (h)",
        "OD (calibrated)" if calibrated else "OD reading",
    )
    if sc is not None:
        # Wider pad than the other grids: the right-hand mu axis sits in between.
        _r2_colorbar(fig, axes, sc, pad=0.07, fraction=0.02)
    _multiples_header(
        fig,
        axes,
        "OD and growth rate over time, by reactor" + (", grouped" if row_labels else ""),
        [
            Line2D([], [], color=OD_TRACE_COLOR, linewidth=1.8, label="OD used for fits (left)"),
            Line2D([], [], color=OD_EXCLUDED_COLOR, linewidth=1.8, label="OD excluded (left)"),
            _growth_rate_handle(),
            *rule_handles(references, show_average),
        ],
        row_labels,
    )
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor=SURFACE, bbox_inches="tight")
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
    fig.savefig(out_path, dpi=FIGURE_DPI, facecolor=SURFACE)
    plt.close(fig)

    summary = (
        merged.groupby("group")[RATE_COL].agg(["count", "mean", "std"]).round(4).reset_index()
    )
    print("\nPer-group growth rate over all segments (h^-1):")
    print(summary.to_string(index=False))


# =============================================================================
# Input discovery and command line
# =============================================================================

# A Pioreactor export folder holds one subfolder per exported table. These two
# are what the analysis cannot run without, so they identify such a folder.
EXPORT_MARKER_SUBDIRS = ("od_readings", "dosing_automation_events")

# How deep below --data-dir to look for the export folder, so that pointing at a
# parent (input/) that holds one export (input/260720_pio_data/) also works.
EXPORT_SEARCH_DEPTH = 2

# Looked for beside the data when --groups is not given.
GROUPS_FILENAME = "groups.csv"


def is_export_dir(path: Path) -> bool:
    return all((path / sub).is_dir() for sub in EXPORT_MARKER_SUBDIRS)


def find_export_dir(root: Path) -> Path:
    """Resolve --data-dir to the export folder holding od_readings/ etc.

    The folder itself is used when it is one; otherwise the subfolders are
    searched, so pointing at a parent that holds a single export works too.
    """
    if not root.is_dir():
        raise FileNotFoundError(f"--data-dir: no such directory: {root}")
    if is_export_dir(root):
        return root

    for depth in range(1, EXPORT_SEARCH_DEPTH + 1):
        found = sorted(
            p for p in root.glob("/".join(["*"] * depth)) if p.is_dir() and is_export_dir(p)
        )
        if len(found) == 1:
            return found[0]
        if len(found) > 1:
            names = "\n  ".join(str(p) for p in found)
            raise ValueError(
                f"{root} holds {len(found)} export folders, so the input is ambiguous. "
                f"Pass --data-dir with one of them:\n  {names}"
            )

    expected = "/, ".join(EXPORT_MARKER_SUBDIRS)
    raise FileNotFoundError(
        f"{root} is not a Pioreactor export folder and holds none: expected a "
        f"folder with the subfolders {expected}/ in it (see input/260720_pio_data)."
    )


def find_groups_file(data_dir: Path | None, root: Path | None) -> Path | None:
    """Look for a group file next to the data, then in the working directory.

    Grouping is optional, so a missing file is not an error: None means "plot
    the reactors ungrouped".
    """
    candidates = []
    for base in (data_dir, data_dir.parent if data_dir else None, root, Path.cwd()):
        if base is not None and base not in candidates:
            candidates.append(base)
    for base in candidates:
        candidate = base / GROUPS_FILENAME
        if candidate.is_file():
            return candidate
    return None


def find_calibration_file(data_dir: Path | None, root: Path | None) -> Path | None:
    """Look for an OD calibration file: a CSV with "calibration" in its name.

    Searched in od_readings/, the export folder, the folder holding it and the
    working directory, in that order; the first folder with one wins. Two in the
    same folder is ambiguous, so that is an error rather than a guess.
    """
    candidates = []
    bases = (
        data_dir / "od_readings" if data_dir else None,
        data_dir,
        data_dir.parent if data_dir else None,
        root,
        Path.cwd(),
    )
    for base in bases:
        if base is not None and base.is_dir() and base.resolve() not in candidates:
            candidates.append(base.resolve())
    for base in candidates:
        found = sorted(f for f in base.glob("*.csv") if "calibration" in f.name.lower())
        if len(found) > 1:
            names = "\n  ".join(str(f) for f in found)
            raise ValueError(
                f"{len(found)} calibration files found, so the choice is ambiguous. "
                f"Pass --calibration to choose one (or --no-calibration):\n  {names}"
            )
        if found:
            return found[0]
    return None


def default_output_dir(data_dir: Path | None) -> Path:
    """Today's results folder, placed beside the export folder it came from.

    Keeping it next to the data means several exports each keep their own
    analysis. With no export folder to sit beside, it lands next to the script.
    """
    base = data_dir.parent if data_dir is not None else Path(__file__).parent
    return base / f"{date.today():%y%m%d}_growth_analysis"


def calculate_growth_rates(
    args, parser, data_dir: Path | None, out_dir: Path, excluded_found: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit every growth segment, write growth_rates.csv and the fit plots.

    Returns the fitted segments and the OD trace they were fitted on (one row
    per reading, with a boolean `used` column marking readings inside a fit).
    """
    try:
        calibration = None
        if args.calibration is not None:
            if not args.calibration.is_file():
                raise FileNotFoundError(f"--calibration: no such file: {args.calibration}")
            calibration = load_calibration(args.calibration)
        od_path = resolve_input_file(
            args.od_file, data_dir, "od_readings", "--od-file", OD_COLUMNS
        )
        events_path = resolve_input_file(
            args.dosing_automation_file,
            data_dir,
            "dosing_automation_events",
            "--dosing-automation-file",
            AUTOMATION_COLUMNS,
        )
        od = load_od_readings(od_path)
        events = load_dosing_automation_events(events_path)
    except (FileNotFoundError, IsADirectoryError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    print(f"OD readings:              {od_path}  ({len(od)} rows)")
    print(f"Dosing automation events: {events_path}  ({len(events)} rows)")

    manual = None
    if not args.ignore_manual_dosing:
        try:
            manual_path = resolve_input_file(
                args.dosing_events_file,
                data_dir,
                "dosing_events",
                "--dosing-events-file",
                MANUAL_COLUMNS,
            )
            manual = load_manual_dosing_events(manual_path)
        except (FileNotFoundError, IsADirectoryError, ValueError) as exc:
            if args.dosing_events_file is not None:
                # The user named this file explicitly, so a failure to read it is an error.
                parser.exit(2, f"error: {exc}\n")
            # Otherwise manual-dosing exclusion is an enhancement; carry on loudly without it.
            print(f"WARNING: manual (UI) dosing events not loaded -> {exc}")
            print("WARNING: OD readings taken during manual pump actions will NOT be excluded.")
        else:
            print(f"Dosing events (manual):   {manual_path}  ({len(manual)} UI rows)")
    print()

    # Excluded reactors are dropped before anything is fitted, so they get no
    # rows in growth_rates.csv and no regression or timeline plots.
    od = drop_excluded(od, args.exclude_reactors, excluded_found)
    events = drop_excluded(events, args.exclude_reactors, excluded_found)
    if manual is not None:
        manual = drop_excluded(manual, args.exclude_reactors, excluded_found)

    calibration_table = None
    if calibration is not None:
        # Excluded reactors need no calibration, and are not worth a note either.
        dropped = {_normalise_unit(u) for u in args.exclude_reactors}
        calibration = calibration[~calibration.index.isin(dropped)]
        try:
            od, calibration_table = calibrate_od(
                od, calibration, args.calibration_skip, args.calibration_points
            )
        except ValueError as exc:
            parser.exit(2, f"error: calibration file {args.calibration}: {exc}\n")
        print(f"OD calibrated with {args.calibration}:")
        for _, row in calibration_table.iterrows():
            print(
                f"  {row['pioreactor_unit']}: measured {row['measured_low']:.3f} -> "
                f"{row['true_low']:.3f}, {row['measured_high']:.3f} -> {row['true_high']:.3f}"
                f"  (OD = {row['slope']:.3f} x reading {row['intercept']:+.3f})"
            )
        print()

    regression_dir = out_dir / "segment_fit_plots"
    reactor_dir = out_dir / "od_plots"
    regression_dir.mkdir(parents=True, exist_ok=True)
    reactor_dir.mkdir(parents=True, exist_ok=True)
    if calibration_table is not None:
        calibration_table.to_csv(out_dir / "od_calibration.csv", index=False)

    n_before = len(od)
    od = od[od["od_reading"] > max(args.min_od, 0.0)]
    if len(od) < n_before:
        print(f"Dropped {n_before - len(od)} OD reading(s) at or below {args.min_od}")

    units = sorted(events["pioreactor_unit"].unique())
    all_results = []
    traces = []
    regression_dpi = FIGURE_DPI if args.high_res_regression_plots else REGRESSION_DPI

    for unit in units:
        od_unit = od[od["pioreactor_unit"] == unit].sort_values(TIME_COL)
        if od_unit.empty:
            print(f"{unit}: no OD readings, skipped")
            continue

        windows = automation_dilution_windows(events, unit)
        n_auto = len(windows)
        n_manual = 0
        if manual is not None:
            manual_wins = manual_dosing_windows(manual, unit)
            n_manual = len(manual_wins)
            windows += manual_wins
        windows = merge_windows(windows)

        segments = build_segments(od_unit, windows, args.min_points)

        # Each unit's regression plots go in their own subfolder.
        unit_regression_dir = regression_dir / str(unit)
        if segments:
            unit_regression_dir.mkdir(parents=True, exist_ok=True)

        unit_results = []
        for seg in segments:
            fit, t, ln_od, result = fit_growth_rate(seg["data"])
            fit["segment_index"] = seg["segment_index"]
            fit["segment_position"] = seg["segment_position"]
            row = {
                "pioreactor_unit": unit,
                "segment_index": seg["segment_index"],
                "segment_position": seg["segment_position"],
                **{k: v for k, v in fit.items() if k not in ("segment_index", "segment_position")},
            }
            unit_results.append(row)
            all_results.append(row)
            plot_regression(unit, fit, t, ln_od, result, unit_regression_dir, regression_dpi)

        # The first readings of a run often spike and squash the rest of the
        # trace, so they can be left out of the OD figures; the fits keep them.
        shown = od_unit.iloc[args.hide_first_readings :]
        plot_reactor_timeline(unit, shown, windows, unit_results, reactor_dir)
        t = shown[TIME_COL]
        used = np.zeros(len(shown), dtype=bool)
        for row in unit_results:
            used |= ((t >= row["start_time_h"]) & (t <= row["end_time_h"])).to_numpy()
        traces.append(shown[[UNIT_COL, TIME_COL, "od_reading"]].assign(used=used))
        print(
            f"{unit}: {len(unit_results)} segments fitted | "
            f"{n_auto} automated dilutions, {n_manual} manual bursts "
            f"-> {len(windows)} merged exclusion windows"
        )

    no_dosing = sorted(set(od["pioreactor_unit"].unique()) - set(units))
    if no_dosing:
        print(
            "\nUnits with OD data but no dosing-automation events "
            f"(not under turbidostat control, skipped): {', '.join(no_dosing)}"
        )

    if not all_results:
        parser.exit(2, "error: no growth segment could be fitted, so there is nothing to plot.\n")

    results_df = pd.DataFrame(all_results)[
        [
            "pioreactor_unit",
            "segment_index",
            "segment_position",
            "start_time_h",
            "end_time_h",
            "duration_h",
            "growth_rate_h-1",
            "std_err_h-1",
            "doubling_time_h",
            "r_squared",
            "p_value",
            "intercept_ln_od",
            "od_start",
            "od_end",
            "n_points",
        ]
    ]
    csv_path = out_dir / "growth_rates.csv"
    results_df.to_csv(csv_path, index=False)

    print(f"\nWrote {len(results_df)} growth-rate rows -> {csv_path}")
    print(f"Segment fit plots (one folder per unit) -> {regression_dir}")
    print(f"OD timeline plots -> {reactor_dir}")
    return results_df, pd.concat(traces, ignore_index=True)


def plot_growth_rates(
    args,
    parser,
    df: pd.DataFrame,
    temps: pd.DataFrame | None,
    out_dir: Path,
    groups_path: Path | None,
    od_trace: pd.DataFrame | None = None,
):
    """Draw every growth-rate figure from the fitted segments."""
    omitted = set()
    if args.omit_reactors:
        # A reactor may appear in only one of the two inputs, so --omit-reactors
        # is resolved against both.
        in_data = set(df[UNIT_COL].unique())
        if temps is not None:
            in_data |= set(temps[UNIT_COL].unique())
        known = {_normalise_unit(u): u for u in in_data}
        omitted = {known[k] for k in (_normalise_unit(u) for u in args.omit_reactors) if k in known}
        unknown = [u for u in args.omit_reactors if _normalise_unit(u) not in known]
        if unknown:
            print(f"WARNING: --omit-reactors: no such reactor in the data: {', '.join(unknown)}")
        if omitted:
            n = len(df)
            df = df[~df[UNIT_COL].isin(omitted)]
            print(f"Omitted {', '.join(sorted(omitted))} ({n - len(df)} segments)")
            if temps is not None:
                temps = temps[~temps[UNIT_COL].isin(omitted)]
            if od_trace is not None:
                od_trace = od_trace[~od_trace[UNIT_COL].isin(omitted)]

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

    scatter_dir = out_dir / "growth_rate_plots"
    references = args.reference_line
    groups = None
    if groups_path is not None:
        try:
            # A reactor with temperature data but no growth rates still needs its
            # group, so both inputs contribute the known reactor list.
            known_units = set(df[UNIT_COL].unique())
            if temps is not None:
                known_units |= set(temps[UNIT_COL].unique())
            groups = load_groups(
                groups_path, sorted(known_units), omitted | set(args.exclude_reactors)
            )
        except (FileNotFoundError, ValueError) as exc:
            parser.exit(2, f"error: {exc}\n")

    plot_per_reactor(df, scatter_dir, norm, references, args.average_line)
    plot_small_multiples(
        df, scatter_dir / "all_reactors.png", norm, references, groups, args.average_line
    )
    print(f"Growth rate plots -> {scatter_dir}")

    if od_trace is not None and not od_trace.empty:
        od_grid_path = out_dir / "od_plots" / "all_reactors.png"
        od_grid_path.parent.mkdir(parents=True, exist_ok=True)
        plot_od_small_multiples(od_trace, od_grid_path, groups, args.calibration is not None)
        print(f"OD grid -> {od_grid_path}")

        od_rate_dir = out_dir / "od_growth_rate_plots"
        if plot_od_with_rates(
            od_trace,
            df,
            od_rate_dir,
            norm,
            references,
            args.average_line,
            args.calibration is not None,
        ):
            plot_od_rates_small_multiples(
                od_trace,
                df,
                od_rate_dir / "all_reactors.png",
                norm,
                groups,
                references,
                args.average_line,
                args.calibration is not None,
            )
            print(f"OD + growth rate plots -> {od_rate_dir}")

    if temps is not None:
        if temps.empty:
            parser.exit(2, "error: no temperature readings left to plot.\n")

        temp_dir = out_dir / "temperature_plots"
        plot_temperature_per_reactor(temps, temp_dir, args.average_line)
        plot_temperature_small_multiples(
            temps, temp_dir / "all_reactors.png", groups, args.average_line
        )
        print(f"Temperature plots -> {temp_dir}")

        combined_dir = out_dir / "temperature_growth_rate_plots"
        plotted = plot_temperature_with_rates(
            temps, df, combined_dir, norm, references, args.average_line
        )
        if plotted:
            print(f"Temperature + growth rate plots -> {combined_dir}")
            missing = sorted(set(temps[UNIT_COL]) ^ set(df[UNIT_COL]))
            if missing:
                print(
                    "NOTE: no combined plot for reactors present in only one of the "
                    f"two inputs: {', '.join(missing)}"
                )
        else:
            print(
                "WARNING: no reactor appears in both the temperature and the "
                "growth-rate data, so no combined plots were written."
            )

    if groups is None:
        # Nothing to group by, so the grouped figure and the comparison are moot.
        print("\nNo group file, so no grouped figure and no group comparison.")
        return

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


HELP_DESCRIPTION = """\
Fit a growth rate (mu, h^-1) to every growth period between dilutions in a
Pioreactor turbidostat export, then write growth_rates.csv and every figure.
Full documentation: README.md"""

HELP_EPILOG = """\
examples:
  growth_analysis.py input/260720_pio_data --groups input/groups.csv
  growth_analysis.py input --exclude-reactors P02        ignore a reactor completely
  growth_analysis.py input --between-dilutions-only --min-r2 0.9
  growth_analysis.py input --growth-rates 260826_growth_analysis/growth_rates.csv"""


class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    # A wider option column keeps each flag's help on one line, without
    # overflowing a narrow terminal.
    def __init__(self, prog):
        width = min(100, shutil.get_terminal_size().columns - 2)
        super().__init__(prog, max_help_position=34, width=width)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [DATA_DIR] [options]",
        description=HELP_DESCRIPTION,
        epilog=HELP_EPILOG,
        formatter_class=_HelpFormatter,
    )

    io = parser.add_argument_group("input and output")
    io.add_argument(
        "data_dir_positional",
        type=str,
        nargs="?",
        metavar="DATA_DIR",
        help="export folder, same as --data-dir (default: current folder)",
    )
    io.add_argument(
        "--data-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="export folder to auto-discover the CSVs in; '' turns discovery off",
    )
    io.add_argument(
        "--groups",
        type=Path,
        default=None,
        metavar="CSV",
        help=f"reactor -> group CSV (default: a {GROUPS_FILENAME} beside the data)",
    )
    io.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="results folder (default: YYMMDD_growth_analysis beside the data)",
    )
    io.add_argument(
        "--growth-rates",
        type=Path,
        default=None,
        metavar="CSV",
        help="redraw figures from an earlier growth_rates.csv, without fitting",
    )

    files = parser.add_argument_group(
        "input files", "Each overrides auto-discovery of that one input in the export folder."
    )
    files.add_argument(
        "--od-file",
        type=Path,
        metavar="CSV",
        help="OD readings (od_readings/)",
    )
    files.add_argument(
        "--dosing-automation-file",
        type=Path,
        metavar="CSV",
        help="dilution start/stop events (dosing_automation_events/)",
    )
    files.add_argument(
        "--dosing-events-file",
        type=Path,
        metavar="CSV",
        help="manual UI pump actions (dosing_events/)",
    )
    files.add_argument(
        "--temperatures",
        type=Path,
        default=None,
        metavar="CSV",
        help="temperature readings (temperature_readings/)",
    )

    reactors = parser.add_argument_group(
        "reactor selection",
        "Names are case-insensitive and a capital O is read as a zero (PO2 = P02).",
    )
    reactors.add_argument(
        "--exclude-reactors",
        nargs="+",
        default=[],
        metavar="UNIT",
        help="ignore completely: not fitted, not in growth_rates.csv, in no plot",
    )
    reactors.add_argument(
        "--omit-reactors",
        nargs="+",
        default=[],
        metavar="UNIT",
        help="fit them, but leave them out of the growth-rate figures",
    )

    calc = parser.add_argument_group("fitting")
    calc.add_argument(
        "--min-points",
        type=int,
        default=4,
        metavar="N",
        help="minimum OD readings needed to fit a segment (default: 4)",
    )
    calc.add_argument(
        "--min-od",
        type=float,
        default=0.0,
        metavar="OD",
        help="drop OD readings at or below this, e.g. blanks (default: 0.0)",
    )
    calc.add_argument(
        "--ignore-manual-dosing",
        action="store_true",
        help="keep OD readings taken during manual UI pump actions",
    )

    cal = parser.add_argument_group(
        "OD calibration",
        "Linearly rescale each reactor's OD so its start and end match measured true ODs.",
    )
    cal.add_argument(
        "--calibration",
        type=Path,
        default=None,
        metavar="CSV",
        help="reactor,od CSV with a low (start) and a high (end) OD per reactor "
        "(default: a *calibration*.csv in od_readings/ or beside the data)",
    )
    cal.add_argument(
        "--no-calibration",
        action="store_true",
        help="use the raw OD even if a calibration file is found",
    )
    cal.add_argument(
        "--calibration-skip",
        type=int,
        default=10,
        metavar="N",
        help="first readings ignored when measuring the start OD (default: 10)",
    )
    cal.add_argument(
        "--calibration-points",
        type=int,
        default=10,
        metavar="N",
        help="readings whose median is the measured start/end OD (default: 10)",
    )

    filters = parser.add_argument_group(
        "plot filters", "These shape the figures only; growth_rates.csv keeps every segment."
    )
    filters.add_argument(
        "--min-r2",
        type=float,
        default=0.0,
        metavar="R2",
        help="drop segments whose fit has a lower R^2 (default: 0.0, keep all)",
    )
    filters.add_argument(
        "--between-dilutions-only",
        action="store_true",
        help="drop the start-of-experiment (lag) and end-of-experiment segments",
    )

    looks = parser.add_argument_group("plot appearance")
    looks.add_argument(
        "--reference-line",
        type=float,
        nargs="+",
        metavar="MU",
        help="dashed line at these growth rates (h^-1) on every plot",
    )
    looks.add_argument(
        "--average-line",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="dashed line at each reactor's average growth rate (default: on)",
    )
    looks.add_argument(
        "--hide-first-readings",
        type=int,
        nargs="?",
        const=10,
        default=0,
        metavar="N",
        help="leave each reactor's first N OD readings (default N: 10) out of the "
        "OD figures; fits are unaffected",
    )
    looks.add_argument(
        "--high-res-regression-plots",
        action="store_true",
        help=f"save the per-segment regression plots at {FIGURE_DPI} dpi too "
        f"(default: {REGRESSION_DPI}; slower)",
    )
    looks.add_argument(
        "--no-temperatures",
        action="store_true",
        help="skip the temperature figures even if the readings are there",
    )

    stats_args = parser.add_argument_group("group statistics")
    stats_args.add_argument(
        "--stats",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="compare the groups and annotate the grouped plot (default: on)",
    )
    stats_args.add_argument(
        "--stats-level",
        choices=("reactor", "segment"),
        default="reactor",
        help="replicate = each reactor's mean (default) or each segment",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if args.calibration_skip < 0 or args.calibration_points < 1:
        parser.error("--calibration-skip must be >= 0 and --calibration-points >= 1")
    if args.hide_first_readings < 0:
        parser.error("--hide-first-readings must be >= 0")

    if args.data_dir is not None and args.data_dir_positional is not None:
        parser.exit(2, "error: the data directory was given both positionally and as --data-dir.\n")
    given = args.data_dir if args.data_dir is not None else args.data_dir_positional
    # `--data-dir ''` disables auto-discovery, so every input must be named explicitly.
    explicit = given is not None
    if not explicit:
        root = Path.cwd()
    elif given.strip():
        root = Path(given.strip())
    else:
        root = None  # auto-discovery switched off with --data-dir ''

    data_dir = None
    if root is not None:
        try:
            data_dir = find_export_dir(root)
        except (FileNotFoundError, ValueError) as exc:
            # Only a directory the user actually named is worth failing over; the
            # default is a guess, and the --*-file flags may cover everything.
            if explicit:
                parser.exit(2, f"error: {exc}\n")
        else:
            print(f"Export folder:            {data_dir}")

    # Created by whichever stage writes first, so a failed run leaves no folder.
    out_dir = args.output_dir or default_output_dir(data_dir)
    print(f"Output folder:            {out_dir}\n")

    excluded_found: set[str] = set()
    od_trace = None  # only a fresh calculation has the OD readings to draw
    if args.growth_rates is not None:
        try:
            df = load_growth_rates(args.growth_rates)
        except (FileNotFoundError, ValueError) as exc:
            parser.exit(2, f"error: {exc}\n")
        print(f"Growth rates: {args.growth_rates}  ({len(df)} segments)")
        if args.calibration is not None:
            print("NOTE: --calibration is ignored with --growth-rates, which does not refit.")
        df = drop_excluded(df, args.exclude_reactors, excluded_found)
        if df.empty:
            parser.exit(2, "error: no growth rates left after --exclude-reactors.\n")
    else:
        if args.no_calibration:
            args.calibration = None
        elif args.calibration is None:
            try:
                args.calibration = find_calibration_file(data_dir, root)
            except ValueError as exc:
                parser.exit(2, f"error: {exc}\n")
            if args.calibration is not None:
                print(f"Calibration:              {args.calibration}  (found beside the data)")
        rates, od_trace = calculate_growth_rates(args, parser, data_dir, out_dir, excluded_found)
        df = prepare_growth_rates(rates)

    temps = None
    if not args.no_temperatures:
        try:
            temp_path = resolve_input_file(
                args.temperatures, data_dir, "temperature_readings", "--temperatures"
            )
            temps = load_temperatures(temp_path)
        except (FileNotFoundError, IsADirectoryError, ValueError) as exc:
            if args.temperatures is not None:
                # Named explicitly, so failing to read it is an error.
                parser.exit(2, f"error: {exc}\n")
            print(f"\nNOTE: no temperature readings found, so no temperature plots -> {exc}")
        else:
            print(
                f"\nTemperatures:             {temp_path}  ({len(temps)} readings, "
                f"{temps[UNIT_COL].nunique()} reactors)"
            )
            temps = drop_excluded(temps, args.exclude_reactors, excluded_found)

    if args.exclude_reactors:
        if excluded_found:
            print(f"\nExcluded from the whole analysis: {', '.join(sorted(excluded_found))}")
        unknown = [u for u in args.exclude_reactors if _normalise_unit(u) not in excluded_found]
        if unknown:
            print(f"WARNING: --exclude-reactors: no such reactor in the data: {', '.join(unknown)}")

    groups_path = args.groups
    if groups_path is None:
        groups_path = find_groups_file(data_dir, root)
        if groups_path is not None:
            print(f"\nGroups:                   {groups_path}  (found beside the data)")
        else:
            print(
                f"\nNOTE: no --groups and no {GROUPS_FILENAME} beside the data, so the "
                "reactors are plotted ungrouped."
            )
    elif not groups_path.is_file():
        parser.exit(2, f"error: --groups: no such file: {groups_path}\n")

    print()
    plot_growth_rates(args, parser, df, temps, out_dir, groups_path, od_trace)
    print(f"\nEverything written to {out_dir}")


if __name__ == "__main__":
    main()
