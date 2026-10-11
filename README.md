# Turbidostat growth-rate calculator

Computes the specific growth rate (µ, h⁻¹) for every growth period between
dilutions in a Pioreactor turbidostat run, plus the periods at the start and end
of the experiment.

## Usage

One command does the whole analysis — fit the growth rates and draw every
figure. Point it at a Pioreactor export folder and a group file:

```bash
uv run growth_analysis.py input/260720_pio_data --groups input/groups.csv
```

Every input CSV is discovered inside the export folder, so nothing else needs
naming. A folder that merely *holds* one export works too (up to two levels
deep), which is why this is equivalent:

```bash
uv run growth_analysis.py input --groups input/groups.csv
```

`--groups` is optional: leave it off and a `groups.csv` sitting beside the data
(or in the working directory) is used. With no group file anywhere, the
reactors are simply plotted ungrouped — the per-reactor and small-multiples
figures are still written, the grouped figure and the group comparison are
skipped.

Results go to a `YYMMDD_growth_analysis/` folder (today's date) created next to
the export folder so each export keeps its own analysis and runs on different 
days never overwrite each other.
Override with `--output-dir`.

Any `--*-file` flag overrides auto-discovery for that one input, so the two can
be mixed — here only the OD file is pinned:

```bash
uv run growth_analysis.py input --groups input/groups.csv \
    --od-file exports/od_readings-run3.csv
```

To rely on nothing but explicit files, pass `--data-dir ''` to switch
auto-discovery off:

```bash
uv run growth_analysis.py --data-dir '' \
    --groups                  input/groups.csv \
    --od-file                 exports/od_readings-run3.csv \
    --dosing-automation-file  exports/dosing_automation_events-run3.csv \
    --dosing-events-file      exports/dosing_events-run3.csv \
    --output-dir              results/run3
```

The run prints which file it used for each input, so it is always on record
which data produced a given `growth_rates.csv`.

Common variations:

```bash
uv run growth_analysis.py input --between-dilutions-only --min-r2 0.9
uv run growth_analysis.py input --reference-line 0.75 0.95  # target / expected µ
uv run growth_analysis.py input --omit-reactors P02         # hide a reactor from the figures
uv run growth_analysis.py input --exclude-reactors P02      # ignore a reactor completely
uv run growth_analysis.py input --no-calibration            # use raw OD despite a calibration file
uv run growth_analysis.py input --no-average-line           # hide the reactor average
uv run growth_analysis.py input --no-temperatures           # skip the temperature figures
uv run growth_analysis.py input --hide-first-readings       # hide the start-up spikes in OD plots
```

To redraw the figures from a `growth_rates.csv` an earlier run produced, pass
`--growth-rates`; the calculation step and its plots are then skipped:

```bash
uv run growth_analysis.py input \
    --growth-rates 260826_growth_analysis/growth_rates.csv --min-r2 0.95
```

### Input folder

`--data-dir` expects a Pioreactor export as exported from the UI: one subfolder
per table, each holding a single CSV. `input/260720_pio_data` is the shipped
example:

```
260720_pio_data/
├── od_readings/                   <- required
├── dosing_automation_events/      <- required
├── dosing_events/                 <- optional: manual UI pump actions
├── temperature_readings/          <- optional: adds the temperature figures
└── ...                            <- anything else in the export is ignored
```

A folder is recognised as an export by having `od_readings/` and
`dosing_automation_events/` in it. The other two are used when present:

- no `dosing_events/` — a warning, and OD readings taken during manual pump
  actions are not excluded;
- no `temperature_readings/` — a note, and the temperature figures are skipped.

If a subfolder contains more than one CSV, the ones lacking that input's
required columns are set aside (so an OD calibration file can live in
`od_readings/`). If that still leaves more than one, the run stops and lists the
candidates so you can pick one with the matching `--*-file` flag. Required
columns are checked up front, so a mis-specified file fails immediately with the
column list rather than part-way through the analysis. If you named an input
explicitly and it cannot be read, that is always an error, even for the optional
ones.

### Options

| Flag | Default | Purpose |
| --- | --- | --- |
| `--data-dir` (or the first positional argument) | current working directory | Export folder to auto-discover inputs in. `''` disables discovery. |
| `--groups` | a `groups.csv` beside the data | two-column CSV mapping each reactor to a group; without one the reactors are plotted ungrouped |
| `--output-dir` | `YYMMDD_growth_analysis` beside the export folder | where the CSV and all plots go |
| `--od-file` | auto | CSV of OD readings (`<data-dir>/od_readings/`) |
| `--dosing-automation-file` | auto | CSV of `DosingStarted`/`DosingStopped` events (`<data-dir>/dosing_automation_events/`) |
| `--dosing-events-file` | auto | CSV holding manual UI pump actions (`<data-dir>/dosing_events/`) |
| `--temperatures` | auto | CSV of temperature readings (`<data-dir>/temperature_readings/`) |
| `--no-temperatures` | off | skip the temperature figures even if the readings are there |
| `--hide-first-readings [N]` | off | leave each reactor's first N OD readings (10 if N is omitted) out of the OD figures; fits are unaffected. Put it after the data folder, or write `--hide-first-readings=10` |
| `--high-res-regression-plots` | off | save the per-segment fit plots (`segment_fit_plots/`) at 300 dpi like every other figure (default 150 dpi; slower) |
| `--growth-rates` | none | plot this existing `growth_rates.csv` instead of recalculating |
| `--min-points` | `4` | minimum OD readings needed to fit a segment |
| `--min-od` | `0.0` | drop OD readings at or below this value before fitting |
| `--ignore-manual-dosing` | off | keep OD readings taken during manual UI pump actions |
| `--calibration` | a `*calibration*.csv` beside the data | CSV of true start/end ODs per reactor; OD is rescaled linearly before fitting (see below) |
| `--no-calibration` | off | use the raw OD even if a calibration file is found |
| `--calibration-skip` | `10` | first readings ignored when measuring a reactor's start OD |
| `--calibration-points` | `10` | readings whose median is taken as the measured start / end OD |
| `--min-r2` | `0.0` | drop segments whose fit is worse than this (plots only) |
| `--between-dilutions-only` | off | keep only true turbidostat cycles (plots only) |
| `--omit-reactors` | none | reactors to leave out of the growth-rate figures, e.g. `--omit-reactors P02` |
| `--exclude-reactors` | none | reactors to ignore completely: not fitted, not in `growth_rates.csv`, in no plot |
| `--reference-line` | none | draw a dashed horizontal line at one or more µ values on every plot |
| `--stats` / `--no-stats` | on | compare the groups statistically and annotate the grouped plot |
| `--stats-level` | `reactor` | unit of replication for the comparison: `reactor` or `segment` |
| `--average-line` / `--no-average-line` | on | dashed line at each reactor's average growth rate |

The filters (`--min-r2`, `--between-dilutions-only`, `--omit-reactors`) shape
the figures only: `growth_rates.csv` always holds every fitted segment.
`--exclude-reactors` is different: the named reactors are dropped from every
input as soon as it is read, so they are never fitted, get no rows in
`growth_rates.csv`, and have no regression, timeline, scatter or temperature
plots. Use it for a reactor that failed or was not part of the experiment.

## OD calibration

The Pioreactor's OD signal is in its own units, which differ from reactor to
reactor. With `--calibration` each reactor's readings are rescaled linearly onto
true OD values measured offline, before anything else is done with them.

The calibration file is found automatically: any CSV with `calibration` in its
name, looked for in `od_readings/`, then the export folder, the folder holding
it, and the working directory. The first folder with one wins, and two in the
same folder is an error. The run prints which file it used. Name one explicitly
with `--calibration`, or turn calibration off with `--no-calibration`:

```bash
uv run growth_analysis.py 1H --calibration exports/261011_calibration.csv
```

The file has two columns, a reactor and an OD, and two rows per reactor in any
order. The smaller value is the true OD at the start of the run, the larger the
true OD at its end:

```csv
reactor,od
P01,0.063
P01,0.246
```

For each reactor, two measured values are matched to these:

- **start** — the median of the readings just after the first
  `--calibration-skip` (10), which usually read high. `--calibration-points` (10)
  readings are used.
- **end** — the median of the last `--calibration-points` readings, so a few
  outliers at the very end do not move it.

The line through those two points (`OD = slope × reading + intercept`) is then
applied to every reading of that reactor. The fitted lines are printed and saved
to `od_calibration.csv`, and `od_start` / `od_end` in `growth_rates.csv` and the
timeline plots are in calibrated units. Calibrated values at or below
`--min-od` (at least 0) are dropped, as uncalibrated ones would be.

Because of the intercept, calibration changes the fitted growth rates and not
just the OD scale: ln(OD) is no longer the raw ln(reading) plus a constant.
Every analysed reactor needs calibration values; leave out any that lack them
with `--exclude-reactors`. With `--growth-rates` nothing is refitted, so
`--calibration` has no effect there.

## Method

0. **Calibrate (optional).** With `--calibration`, rescale each reactor's OD
   linearly onto its true start and end OD (see above).
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

Every figure is saved at 300 dpi, ready for posters and slides, except the
per-segment regression plots: there are hundreds of them, so they default to
150 dpi unless `--high-res-regression-plots` is given.

Everything lands in one dated folder beside the export folder,
`YYMMDD_growth_analysis/`:

- `od_calibration.csv` — with `--calibration` only: per reactor, the measured
  and true start/end OD and the slope and intercept applied.
- `growth_rates.csv` — one row per segment:
  `pioreactor_unit`, `segment_index`, `segment_position`
  (`experiment_start` / `between_dilutions` / `experiment_end`),
  `start_time_h`, `end_time_h`, `duration_h`, `growth_rate_h-1`,
  `std_err_h-1`, `doubling_time_h`, `r_squared`, `p_value`,
  `intercept_ln_od`, `od_start`, `od_end`, `n_points`.
- `segment_fit_plots/<unit>/<unit>_segment_<nnn>.png` — ln(OD) vs time with
  the fitted line, annotated with µ ± SE, R², doubling time, n and duration.
  One subfolder per reactor.
- `od_plots/<unit>_od_timeline.png` — full OD trace per reactor.
  Blue = readings used in fits, grey = excluded, red = automated dilution
  windows, orange = manual dosing windows.
- `od_plots/all_reactors.png` — every reactor's OD trace as small
  multiples on shared axes, laid out like `growth_rate_plots/all_reactors.png`
  (one row per group). Blue = readings used in fits, grey = excluded. Not drawn
  with `--growth-rates`, which has no OD readings; `--omit-reactors` applies.
- `od_growth_rate_plots/<unit>.png` — one plot per reactor with the OD trace on the
  left axis (dark = used in fits, light grey = excluded) and that reactor's
  growth rates as R²-shaded points on the right axis, with the reactor-average
  and any `--reference-line` rules. Like the temperature version, but showing
  the OD each µ was fitted from. Also not drawn with `--growth-rates`.
- `od_growth_rate_plots/all_reactors.png` — the same OD + growth-rate panels as small
  multiples, laid out like the other `all_reactors.png` grids (one row per
  group). Both y axes are shared across panels, so OD and µ compare directly
  between reactors; the µ ticks are labelled on the last panel of each row.

## The growth-rate figures

Growth rates are shaded by the R² of their fit, so a poor regression is visible
rather than silently averaged in:

- `growth_rate_plots/<unit>.png` — µ per segment against time, one plot per
  reactor, points shaded by the R² of their fit, dashed line = reactor average.
  A legend names the dashed rules on both the per-reactor plots and the
  small-multiples figure.
- `growth_rate_plots/all_reactors.png` — the same panels as small multiples on
  shared axes, one row per group (labelled with the group name and its reactor
  count), for comparing reactors within and across groups at a glance. Reactors
  missing from the group file get a trailing `(no group)` row.
- `growth_rates_by_group.png` — every reactor side by side (box + individual
  segments), ordered, shaded and labelled by group, with the pooled per-group
  distribution beside it. Per-group mean, SD and count are also printed.

### Temperature

Temperature readings — the columns `pioreactor_unit`, `temperature_c` and
`hours_since_experiment_created` of the Pioreactor temperature export — are
picked up from `<data-dir>/temperature_readings/`, or named with
`--temperatures`. When they are there, these are written too (turn them off
with `--no-temperatures`):

- `temperature_plots/<unit>.png` — temperature against time, one plot per reactor,
  dashed line = reactor average.
- `temperature_plots/all_reactors.png` — the same traces as small multiples on shared
  axes, laid out exactly like `growth_rate_plots/all_reactors.png`: one row per
  group, with a trailing `(no group)` row.
- `temperature_growth_rate_plots/<unit>.png` — one plot per reactor with the
  temperature trace on the left axis and that reactor's growth rates as points
  on the right axis, so a µ can be read against the temperature it was measured
  at. Only reactors present in both inputs get one; any reactor found in just
  one of them is named in a note.

`--omit-reactors` and `--average-line` apply to these figures too, and
`--omit-reactors` accepts a reactor that appears only in the temperature data.
The growth-rate filters (`--min-r2`, `--between-dilutions-only`) affect the
points on the combined plot, not the temperature traces.

### Comparing the groups

With two or more groups, `compare_groups()` tests them against each other and the
result is drawn on `growth_rates_by_group.png` as significance brackets over the
pooled panel, with a footnote naming the test. The same numbers are printed.

- **Two groups** — Welch's t-test (unequal variances).
- **More than two** — a one-way ANOVA across all groups, plus pairwise Welch
  tests with Holm-Bonferroni correction. Only the first three brackets are
  drawn, or the significant ones if there are more; the full table is printed.
- **`--stats-level reactor` (default)** — each reactor contributes one value,
  its average growth rate. This is the honest unit of replication: segments from
  one reactor are repeated measures of the same culture.
- **`--stats-level segment`** — pools every segment instead. It gives far
  smaller p-values, but they answer "do these segments differ?", not "do these
  treatments differ?". Use it only to describe the spread, never to claim an
  effect.

With three reactors per group, the reactor-level test has very little power: a
real effect can easily come back `ns`. Read the effect size and the per-reactor
panel alongside the p-value.

`--omit-reactors` drops the named reactors before anything is plotted, so they
appear in neither the scatters, the grouped figure, nor the per-group summary.
Names are matched case-insensitively and a capital `O` is accepted for a zero; a
name that is not in the data is reported as a warning rather than failing the
run. Note that a per-reactor PNG written by an earlier run is not deleted.

`--reference-line` takes one or more growth rates and draws each as a dashed
violet rule, labelled with its own value, on every figure — a target µ, a
published value, or a previous run's mean. It is distinct from the orange dashed
reactor-average line, and the y range is widened if needed so the line is visible.

The group file has two columns, one row per reactor:

```csv
group,reactor
P01,Heated
P02,Heated
```

Which column holds the reactor is detected from the contents rather than the
header (the shipped `input/groups.csv` has its headers the other way round), and a
capital `O` typed for a zero — `PO8` — is matched to `P08`, with a note printed.
Reactors missing from the group file are left out of the grouped plot. Without
any group file, `growth_rates_by_group.png` and the comparison are not produced
and the small multiples wrap three per row instead of one row per group.

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
df = pd.read_csv("260826_growth_analysis/growth_rates.csv")
clean = df[(df.segment_position == "between_dilutions") & (df.r_squared > 0.9)]
print(clean.groupby("pioreactor_unit")["growth_rate_h-1"].agg(["count", "mean", "std"]))
```