# Turbidostat growth-rate calculator

Computes the specific growth rate (µ, h⁻¹) for every growth period between
dilutions in a Pioreactor turbidostat run, plus the periods at the start and end
of the experiment.

## Usage

Point it at an export directory and let it find the CSVs (each subfolder must
hold exactly one):

```bash
uv run calculate_growth_rates.py --data-dir /path/to/pio_export
```

Or name the input files directly. Any `--*-file` flag overrides auto-discovery
for that input, so you can mix the two — here only the OD file is pinned:

```bash
uv run calculate_growth_rates.py --od-file /path/to/od_readings-run3.csv
```

To rely on nothing but explicit files, pass `--data-dir ''` to switch
auto-discovery off:

```bash
uv run calculate_growth_rates.py --data-dir '' \
    --od-file                 exports/od_readings-run3.csv \
    --dosing-automation-file  exports/dosing_automation_events-run3.csv \
    --dosing-events-file      exports/dosing_events-run3.csv \
    --output-dir              results/run3
```

The run prints which file it used for each input, so it is always on record
which data produced a given `growth_rates.csv`.

### Options

| Flag | Default | Purpose |
| --- | --- | --- |
| `--data-dir` | the `20260710_pio_export` folder | Export directory to auto-discover inputs in. `''` disables discovery. |
| `--od-file` | auto | CSV of OD readings (`<data-dir>/od_readings/`) |
| `--dosing-automation-file` | auto | CSV of `DosingStarted`/`DosingStopped` events (`<data-dir>/dosing_automation_events/`) |
| `--dosing-events-file` | auto | CSV holding manual UI pump actions (`<data-dir>/dosing_events/`) |
| `--output-dir` | `./output` | where the CSV and plots go |
| `--min-points` | `4` | minimum OD readings needed to fit a segment |
| `--min-od` | `0.0` | drop OD readings at or below this value before fitting |
| `--ignore-manual-dosing` | off | keep OD readings taken during manual UI pump actions |

Required columns are checked up front, so a mis-specified file fails
immediately with the column list rather than part-way through the analysis.
If a subfolder contains more than one CSV, the run stops and lists the
candidates so you can pick one with the matching `--*-file` flag.

A missing `dosing_events` input is only a warning — the analysis continues
without manual-dosing exclusion, and says so. But if you named
`--dosing-events-file` explicitly and it cannot be read, that is an error.

## Method

1. **Find dilution windows.** Each `DosingStarted` event in
   `dosing_automation_events/` is paired with the following `DosingStopped`.
   The interval between them is a dilution window. A trailing `DosingStarted`
   with no `DosingStopped` (export taken mid-dilution) is treated as open-ended,
   discarding everything after it.
2. **Add manual dosing windows.** `dosing_events/` also contains
   `source_of_event == "UI"` rows — manual pump actions that perturb OD but are
   *not* recorded in `dosing_automation_events/`. Consecutive rows less than
   0.05 h apart are collapsed into one window. Disable with
   `--ignore-manual-dosing`.
3. **Exclude and split.** Every OD reading falling inside any window
   (bounds inclusive) is discarded. Each remaining maximal contiguous run of
   readings is one growth segment.
4. **Fit.** For each segment with at least `--min-points` readings,
   `ln(OD)` is regressed on `hours_since_experiment_created`. The slope is µ in
   h⁻¹; doubling time is `ln(2)/µ`.

Only units that appear in `dosing_automation_events/` are analysed. Units with
OD data but no dosing events are reported and skipped.

## Outputs

- `output/growth_rates.csv` — one row per segment:
  `pioreactor_unit`, `segment_index`, `segment_position`
  (`experiment_start` / `between_dilutions` / `experiment_end`),
  `start_time_h`, `end_time_h`, `duration_h`, `growth_rate_h-1`,
  `std_err_h-1`, `doubling_time_h`, `r_squared`, `p_value`,
  `intercept_ln_od`, `od_start`, `od_end`, `n_points`.
- `output/regression_plots/<unit>_segment_<nnn>.png` — ln(OD) vs time with the
  fitted line, annotated with µ ± SE, R², doubling time, n and duration.
- `output/reactor_plots/<unit>_od_timeline.png` — full OD trace per reactor.
  Blue = readings used in fits, grey = excluded, red = automated dilution
  windows, orange = manual dosing windows.

## Interpreting the results

Not every segment is exponential growth, so filter before averaging:

- `segment_position == "between_dilutions"` selects true turbidostat cycles.
  The `experiment_start` segment includes the lag phase and near-zero blank
  readings, so its R² is usually poor.
- Long segments (`duration_h` of several hours) are stretches where the
  turbidostat stopped diluting — the culture was at steady state or stalled, not
  growing exponentially. They fit poorly and should normally be dropped.
- `r_squared` is the main quality filter; a low value means the segment was not
  log-linear.

A reasonable summary per reactor:

```python
import pandas as pd
df = pd.read_csv("output/growth_rates.csv")
clean = df[(df.segment_position == "between_dilutions") & (df.r_squared > 0.9)]
print(clean.groupby("pioreactor_unit")["growth_rate_h-1"].agg(["count", "mean", "std"]))
```

## Notes on this dataset (`20260710_pio_export`)

- Reactors under turbidostat control: P01–P06, P08. P07, P09, P10, P11 have OD
  data but no dosing events.
- Many dilutions fire back-to-back (< 3 min apart) when OD stays above target;
  those gaps hold too few readings to fit and are skipped via `--min-points`.
- After roughly 45.5 h OD climbs well above target in several reactors despite
  dilution commands, suggesting media exhaustion or a pump problem. Those data
  fall inside dilution windows and are excluded.
