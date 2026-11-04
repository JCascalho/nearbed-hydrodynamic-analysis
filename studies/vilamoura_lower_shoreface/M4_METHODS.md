# M4 — Near-Bed Hydrodynamics and Tracer Mobility

**Vilamoura lower-shoreface study (Portugal)**

This study-specific adaptation of the author's validated `M3_oceanographic_data_processing_T_equiv.py` (October 2026) quantifies combined wave–current forcing near the seabed and compares maximum bed shear stress with critical thresholds relevant to fluorescent tracer and native sediment. **The numerical formulae, sampling intervals, threshold values, and core analysis logic have not been modified.** Only the module name, documentation, and default output filename were adapted.

## What the module calculates

- Current magnitude and azimuth; current-related friction factor, shear velocity (`u_star_current`) and shear stress (`tau_c`).
- Wave orbital parameters, equivalent period (`T_equiv = 1.281 × Tz`), wave friction and wave shear stress (`tau_w`).
- Wave–current angle, mean combined shear stress (`tau_m`), maximum combined shear stress (`tau_max`) and equivalent combined-flow shear velocity (`u_star_max`).
- Per-interval stress statistics, threshold-exceedance record counts and filtered observation sheets.
- Native-sediment mixing-layer estimates (`delta_mix`, `delta_mix_p95_tau`, `delta_mix_mean_tau`) from Harris and Wiberg (1997), with maximum, 95th-percentile, or mean interval `tau_max`, respectively.

## Critical distinction: tracer and native thresholds

The program evaluates configured `tau_max` exceedance thresholds of **0.259 N m⁻² (tracer)** and **0.303 N m⁻² (native sediment)**, together with a zero reference threshold. The `tau_cr` column in the input represents **native sediment** critical stress. Accordingly, `tau_excess`, `is_mobilized_native`, `u_star_cr`, and all `delta_mix` estimates are **native-threshold-based**; they are **not** tracer-only mobility diagnostics.

The `tau_max_thresholds` summaries do identify records for which `tau_max > 0.259 N m⁻²`, but counts are **individual observations**, not independent mobility episodes or verified tracer transport events. Consecutive records must be grouped explicitly to estimate independent events and their durations. A tracer-exceedance percentage can be calculated from the valid record count for an interval, but this would not automatically equal a time-weighted duration if record spacing or coverage is irregular.

## Required input

An Excel workbook with a worksheet containing the columns:

| Field | Description | Units |
| --- | --- | --- |
| `Time` | Observation timestamp | date/time |
| `ve`, `vn` | Eastward and northward currents | m s⁻¹ |
| `D50` | Native-sediment median grain diameter | m |
| `z` | Current sensor height above bed | m |
| `h` | Water depth | m |
| `Tz` | Mean zero-crossing wave period | s |
| `Hs` | Significant wave height | m |
| `dir_w` | Wave direction | degrees |
| `tau_cr` | Native-sediment critical bed shear stress | N m⁻² |

`Tau_cr` and `TAU_CR` are accepted and normalized. The original code applies a basic required-column check; it is not a comprehensive physical-input quality-control implementation.

## Synthetic demonstration (no experimental data)

The provided `examples/INPUT_CURRENTS_WAVES_SYNTHETIC.xlsx` contains invented values with dates falling within the study's configured six monitoring intervals. They demonstrate code execution **only** and must not be used as Vilamoura observations, manuscript results, or event statistics.

Install dependencies and run from this case-study directory:

```bash
python -m pip install numpy pandas openpyxl
python M4_nearbed_hydrodynamics_tracer_mobility.py --input examples/INPUT_CURRENTS_WAVES_SYNTHETIC.xlsx --output M4_DEMO_RESULTS.xlsx
```

For private experimental use, point `--input` to the actual workbook. The script will overwrite an existing output file with the same name.

## Outputs

- `C0_C7`: complete processed time series, with all derived variables.
- `tau_max_summary`: interval summaries, exceedance fractions for the **native** threshold, and mixing-layer statistics.
- `tau_max_thresholds`: interval/threshold summaries and **record counts** for `tau_max > 0`, `> 0.259`, and `> 0.303` N m⁻² when present.
- `C0_C2` through `C6_C7`: sampled interval records.
- Threshold-specific worksheets: records exceeding each configured threshold.

A nominal native threshold of 0.303 N m⁻² and a tracer threshold of 0.259 N m⁻² were used in the original Vilamoura configuration. `delta_mix_p95_tau` is the estimate subsequently used for tracer-derived sediment transport calculations in the manuscript, but is calculated with the **native** threshold.

## Availability policy

The full Vilamoura experimental current–wave database and corresponding numerical results are **not included** in this repository package and are intended for release following publication of the research article, subject to the journal's data-availability requirements. Only synthetic demonstration input is supplied here.

## References and provenance

- Bagnold (1956).
- Soulsby and Smallman (1986).
- Harris and Wiberg (1997).
- Soulsby and Clarke (2005).
- Soulsby (2006), equivalent-period approach.

These author/year citations reproduce the references stated in the original script; consult the corresponding manuscript bibliography for full bibliographic details.

**Author:** João Cascalho, FCUL / Instituto Dom Luiz, University of Lisbon.

**Scientific scope:** Combined wave–current potential sediment mobilization; threshold exceedance is not a direct measurement of sediment transport. No scientific claims should be derived from the synthetic example.
