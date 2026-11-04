# M4 — Supplementary Hydrodynamic Visualization Tools (Vilamoura lower shoreface)

These four **independent** Python plotting scripts complement the M4 near-bed hydrodynamic and tracer-mobility calculation. They reproduce the original scientific classification methods and the graphical conventions of the four validated legacy routines. Harmonization is limited to file names, default input selection, default sheet (`C0_C7`) and consistent output locations. Each script still accepts `--input`, `--sheet`, and `--output` (the rose scripts also accept `--legend` and `--normalize`).

## Files and purpose

| File | Scientific purpose | Variables |
| --- | --- | --- |
| `M4a_current_speed_roses.py` | 32-sector current-speed directional frequency rose | `dir_c`, `u_z` |
| `M4b_wave_height_roses.py` | 32-sector significant wave-height directional frequency rose | `dir_w`, `Hs` |
| `M4c_wave_period_roses.py` | 32-sector wave-period directional frequency rose | `dir_w`, `Tz`; optional `Tp` |
| `M4d_along_cross_shelf_currents.py` | Along-shelf and cross-shelf current velocity time series with running mean | `Time`, `ve`, `vn` |

The rose scripts use **32 compass sectors**, centred on the indicated bearings, with sector width 11.25°. Their original fixed classes are preserved: current speed <0.10, 0.10–0.20, 0.20–0.30, 0.30–0.40, ≥0.40 m/s; Hs <0.2, 0.2–0.4, 0.4–0.6, 0.6–0.8, 0.8–1.0, 1.0–1.2, ≥1.2 m; and Tz <3, 3–6, 6–9, 9–12, 12–15, 15–18, ≥18 s. These are **plotting classes**, not mobility thresholds.

For M4a `dir_c` is a current-flow **towards** bearing, while `dir_w` in M4b and M4c represents wave **origin (from)** bearings. They are deliberately **not** made numerically equivalent: no 180° reversal is applied automatically. Verify the convention of any external input.

M4d retains the original **112° clockwise-from-North positive along-shelf axis**, equivalent to the original 22° transformation, with the positive cross-shelf component 90° counter-clockwise from that axis. The original 96-*observation* rolling mean is retained. It represents 24 h only for continuous 15-minute input; irregular sampling or missing records mean it is **not** necessarily a 24-h time-weighted average. Confirm the orientation and sign convention against all other manuscript velocity and tracer components before interpretation.

## Input data policy

`examples/INPUT_HYDRODYNAMICS_SYNTHETIC.xlsx` contains **invented observations** across six artificial campaign sheets (`C0_C2` ... `C6_C7`) and `C0_C7`. These data are not from Vilamoura and are **not suitable for manuscript results**. No experimental time series or figure computed from the original observations is included. The eventual experimental data archive should be deposited separately after the article is published, subject to journal requirements.

## Installation and execution

From the `plotting` directory:

```bash
python -m pip install numpy pandas matplotlib openpyxl pillow
python M4a_current_speed_roses.py
python M4b_wave_height_roses.py
python M4c_wave_period_roses.py
python M4d_along_cross_shelf_currents.py
```

The four scripts read the **synthetic** workbook by default and save images under `outputs/` (not included in the release). Run a particular interval with a specified **private** workbook, for example:

```bash
python M4a_current_speed_roses.py --input "C:/private/OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx" --sheet C2_C3 --output outputs/M4a_C2_C3.png --legend outputs/M4a_legend.png
python M4b_wave_height_roses.py --input "C:/private/OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx" --sheet C2_C3 --output outputs/M4b_C2_C3.png
python M4c_wave_period_roses.py --input "C:/private/OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx" --sheet C2_C3 --output outputs/M4c_C2_C3.png
python M4d_along_cross_shelf_currents.py --input "C:/private/OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx" --sheet C0_C7 --azimuth 112 --window 96 --output outputs/M4d_C0_C7.png
```

For the rose scripts, `--normalize global` (default) gives percentages of all valid records in the chosen sheet; `--normalize per-direction` instead describes classes within each directional sector, and should **not** be interpreted as sector-occurrence frequency. Specify output names when changing sheets to avoid overwriting previous plots. The input `Tp` option in M4c requires a genuine `Tp` column; do not relabel derived `T_equiv` as measured peak period.

**Scientific provenance:** the original scripts by João Cascalho (FCUL / Instituto Dom Luiz, University of Lisbon), November 2025, revised October 2026; harmonized visualization packaging for M4, November 2026. This package does **not** modify the released M4 hydrodynamic computation or its scientific results.
