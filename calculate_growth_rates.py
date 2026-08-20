"""
Calculate per-cycle growth rates from Pioreactor turbidostat export data.

For each pioreactor unit, OD readings are split into growth segments bounded by
dilution events. A dilution runs from a "DosingStarted" event to its matching
"DosingStopped" event; every OD reading inside such a window is discarded,
since dosing disturbs the OD signal. Manual (UI-triggered) pump actions from
dosing_events/ are treated as additional exclusion windows by default, because
they perturb OD the same way but are not recorded in dosing_automation_events/.

For every remaining segment, ln(OD) is regressed against time
(hours_since_experiment_created) to obtain the growth rate (h^-1) and the R^2
of the fit.

Outputs (written to --output-dir):
  - growth_rates.csv          one row per growth segment
  - regression_plots/*.png    ln(OD) vs time with fitted line, per segment
  - reactor_plots/*.png       OD vs time per reactor, exclusion windows shaded
"""

import os
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

TIME_COL = "hours_since_experiment_created"

# Manual UI pump actions are logged as many individual tiny-volume rows. Rows
# closer together than this are collapsed into one exclusion window.
MANUAL_BURST_GAP_H = 0.05

DEFAULT_DATA_DIR = Path(os.getcwd())


OD_COLUMNS = ["pioreactor_unit", "od_reading", TIME_COL]
AUTOMATION_COLUMNS = ["pioreactor_unit", "event_name", TIME_COL]
MANUAL_COLUMNS = ["pioreactor_unit", "source_of_event", TIME_COL]


def resolve_input_file(explicit: Path | None, data_dir: Path | None, subdir: str, flag: str) -> Path:
    # Pick the CSV to read: an explicit --*-file wins, else the one CSV in data_dir/subdir
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


def plot_regression(unit, fit, t, ln_od, result, out_dir: Path):
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
    fig.savefig(out_dir / f"{unit}_segment_{fit['segment_index']:03d}.png", dpi=150)
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
    fig.savefig(out_dir / f"{unit}_od_timeline.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=str(DEFAULT_DATA_DIR),
        help="Pioreactor export directory holding od_readings/, dosing_automation_events/ "
        "and dosing_events/. Each input CSV is auto-discovered inside it unless the "
        "matching --*-file flag is given. Pass --data-dir '' to disable auto-discovery "
        f"and require explicit files. (default: Current working directory)",
    )
    parser.add_argument(
        "--od-file",
        type=Path,
        help="CSV of OD readings. Overrides auto-discovery in <data-dir>/od_readings/.",
    )
    parser.add_argument(
        "--dosing-automation-file",
        type=Path,
        help="CSV of dosing automation events (DosingStarted / DosingStopped). "
        "Overrides auto-discovery in <data-dir>/dosing_automation_events/.",
    )
    parser.add_argument(
        "--dosing-events-file",
        type=Path,
        help="CSV of dosing events, used for manual (UI) pump actions. Overrides "
        "auto-discovery in <data-dir>/dosing_events/. Ignored with --ignore-manual-dosing.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "output",
        help="Where to write growth_rates.csv and the plot folders.",
    )
    parser.add_argument(
        "--min-points",
        type=int,
        default=4,
        help="Minimum OD readings required to fit a segment (default: 4).",
    )
    parser.add_argument(
        "--min-od",
        type=float,
        default=0.0,
        help="Discard OD readings at or below this value before fitting. Useful to drop "
        "near-zero blank readings whose logarithm is dominated by noise (default: 0.0).",
    )
    parser.add_argument(
        "--ignore-manual-dosing",
        action="store_true",
        help="Do not treat manual (UI) pump actions from dosing_events/ as exclusion "
        "windows. By default they ARE excluded, since they perturb OD just like "
        "automated dilutions but are absent from dosing_automation_events/.",
    )
    args = parser.parse_args()

    # `--data-dir ''` disables auto-discovery, so every input must be named explicitly.
    data_dir = Path(args.data_dir) if args.data_dir.strip() else None

    try:
        od_path = resolve_input_file(args.od_file, data_dir, "od_readings", "--od-file")
        events_path = resolve_input_file(
            args.dosing_automation_file,
            data_dir,
            "dosing_automation_events",
            "--dosing-automation-file",
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
                args.dosing_events_file, data_dir, "dosing_events", "--dosing-events-file"
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

    regression_dir = args.output_dir / "regression_plots"
    reactor_dir = args.output_dir / "reactor_plots"
    regression_dir.mkdir(parents=True, exist_ok=True)
    reactor_dir.mkdir(parents=True, exist_ok=True)

    n_before = len(od)
    od = od[od["od_reading"] > max(args.min_od, 0.0)]
    if len(od) < n_before:
        print(f"Dropped {n_before - len(od)} OD reading(s) at or below {args.min_od}")

    units = sorted(events["pioreactor_unit"].unique())
    all_results = []

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
            plot_regression(unit, fit, t, ln_od, result, regression_dir)

        plot_reactor_timeline(unit, od_unit, windows, unit_results, reactor_dir)
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

    results_df = pd.DataFrame(all_results)
    column_order = [
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
    results_df = results_df[column_order]
    csv_path = args.output_dir / "growth_rates.csv"
    results_df.to_csv(csv_path, index=False)

    print(f"\nWrote {len(results_df)} growth-rate rows -> {csv_path}")
    print(f"Regression plots -> {regression_dir}")
    print(f"Reactor timeline plots -> {reactor_dir}")


if __name__ == "__main__":
    main()
