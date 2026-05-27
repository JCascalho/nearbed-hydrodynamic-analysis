# Near-Bed Hydrodynamic Analysis

Python routine for processing synchronized near-bed current and wave records and
calculating hydrodynamic forcing parameters relevant to sediment mobility.

The workflow computes:

- current magnitude and direction;
- current-only shear velocity and shear stress;
- wave orbital velocity and wave shear stress;
- wave-current angle;
- mean and maximum combined bed shear stress;
- combined-flow shear velocity;
- critical shear velocity;
- sediment mobilization flags;
- mixing-layer thickness estimates;
- interval and threshold summary statistics.

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Input Data

The input file can be CSV or Excel.

Required columns:

```text
ve, vn, D50, z, h, Tz, Hs, dir_w, Time, tau_cr
```

Where:

- `ve`, `vn`: near-bed current velocity components;
- `D50`: median native sediment grain size;
- `z`: current measurement elevation above bed;
- `h`: water depth;
- `Tz`: zero-crossing wave period;
- `Hs`: significant wave height;
- `dir_w`: wave direction;
- `Time`: timestamp;
- `tau_cr`: critical shear stress for native sediment.

## Usage

Excel input:

```bash
python nearbed_hydrodynamic_analysis.py --input currents_waves.xlsx --output hydrodynamic_results.xlsx
```

CSV input:

```bash
python nearbed_hydrodynamic_analysis.py --input currents_waves.csv --output hydrodynamic_results.xlsx
```

Use a specific Excel sheet:

```bash
python nearbed_hydrodynamic_analysis.py --input currents_waves.xlsx --sheet Sheet1 --output hydrodynamic_results.xlsx
```

Use custom `tau_max` thresholds:

```bash
python nearbed_hydrodynamic_analysis.py --input currents_waves.xlsx --thresholds 0,0.259,0.303
```

Use custom time intervals:

```bash
python nearbed_hydrodynamic_analysis.py --input currents_waves.xlsx --intervals-csv intervals.csv
```

The intervals CSV must contain:

```text
start,end,label
2020-08-18 17:50:08,2020-08-30 23:50:08,C0_C2
```

## Outputs

The output Excel workbook includes:

- `All_Data`: full processed time series;
- `tau_max_summary`: interval statistics;
- `tau_max_thresholds`: threshold exceedance statistics;
- one sheet per time interval;
- thresholded interval subsets;
- `metadata`;
- `SheetName_Map`.

## Notes

The default intervals and thresholds reproduce the original study workflow, but
can be replaced using command-line arguments. Users should verify physical
constants, roughness assumptions, critical shear-stress values, and interval
definitions for their own deployment and sediment conditions.
