# M4 supplementary visualization — M4d_along_cross_shelf_currents.py
# Harmonized file names, default input and output paths only.
# Scientific classification, graphical conventions, and numerical formulas are retained from M_3d_along_cross_shelf_currents_original_convention(3).py.
# The default input is synthetic. Use --input to process local research data privately.

"""
===============================================================================
Along-shelf and cross-shelf current velocity time-series generator
===============================================================================

Purpose
-------
Generate the time series of along-shelf (vas) and cross-shelf (vcs) current
velocities from eastward (ve) and northward (vn) velocity components in the
original oceanographic Excel workbook. Plot the original 15-minute observations
and a 96-observation (24-hour) moving average, in a two-panel figure.

Input
-----
An Excel workbook with a worksheet containing Time, ve, and vn. The default
input is INPUT_CURRENTS_WAVES_DATA_BASE(1).xlsx, worksheet C0_C7.

Coordinate convention
---------------------
The along-shelf axis is specified as a clockwise azimuth from geographic North
(positive direction). Default: 112 degrees (22 degrees south of East), matching the original 2025 script.
The positive cross-shelf direction is 90 degrees counterclockwise from this
along-shelf axis (toward azimuth 22 degrees for the default axis).

    vas = ve * sin(azimuth) + vn * cos(azimuth)
    vcs = -ve * cos(azimuth) + vn * sin(azimuth)

The default 112-degree azimuth reproduces the original 2025 rotation:
vas = ve*cos(22 deg) - vn*sin(22 deg);
vcs = ve*sin(22 deg) + vn*cos(22 deg).

Usage
-----
python M_3d_along_cross_shelf_currents.py
python M_3d_along_cross_shelf_currents.py --input "data.xlsx" --sheet C0_C7 --azimuth 112 --output "currents.png"

Outputs
-------
A high-resolution figure (PNG, JPEG, or PDF, inferred from extension).
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=str(Path(__file__).resolve().parent / 'examples' / 'INPUT_HYDRODYNAMICS_SYNTHETIC.xlsx'), help='Input Excel workbook')
    parser.add_argument('--sheet', default='C0_C7', help='Worksheet containing Time, ve and vn')
    parser.add_argument('--azimuth', type=float, default=112.0, help='Positive along-shelf axis clockwise from North, degrees (default: 112; original 2025 convention)')
    parser.add_argument('--window', type=int, default=96, help='Moving-average window in observations (default: 96)')
    parser.add_argument('--output', default='outputs/M4d_along_cross_shelf_currents.png', help='Output image filename')
    parser.add_argument('--dpi', type=int, default=300, help='Output image resolution (default: 300)')
    args = parser.parse_args()

    if args.window < 1:
        parser.error('--window must be positive')
    df = pd.read_excel(args.input, sheet_name=args.sheet, usecols=['Time', 've', 'vn'])
    df['Time'] = pd.to_datetime(df['Time'], errors='coerce')
    df['ve'] = pd.to_numeric(df['ve'], errors='coerce')
    df['vn'] = pd.to_numeric(df['vn'], errors='coerce')
    df = df.dropna(subset=['Time', 've', 'vn']).sort_values('Time').drop_duplicates('Time')
    if df.empty:
        raise ValueError('No valid Time, ve and vn observations found.')

    theta = np.deg2rad(args.azimuth)
    df['vas'] = df['ve'] * np.sin(theta) + df['vn'] * np.cos(theta)
    df['vcs'] = -df['ve'] * np.cos(theta) + df['vn'] * np.sin(theta)
    # Rolling window uses 96 consecutive 15-min observations; min_periods=1
    # reproduces the original moving average, including initial partial windows.
    df['vas_avg'] = df['vas'].rolling(args.window, min_periods=1).mean()
    df['vcs_avg'] = df['vcs'].rolling(args.window, min_periods=1).mean()

    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    specs = [
        ('vas', 'vas_avg', 'Along-shelf Velocity (vas)', '#e6ae32'),
        ('vcs', 'vcs_avg', 'Cross-shelf Velocity (vcs)', '#9fe2a3'),
    ]
    for ax, (raw, avg, title, color) in zip(axes, specs):
        ax.plot(df['Time'], df[raw], color=color, lw=0.45, alpha=0.75, label='Original', rasterized=True)
        ax.plot(df['Time'], df[avg], color='black', lw=1.1, label=f'Moving Avg ({args.window} data lines)')
        ax.set_title(title, fontsize=19, fontweight='bold', pad=8)
        ax.set_ylabel('Velocity (m/s)', fontsize=16, fontweight='bold')
        ax.set_ylim(-0.3, 0.3)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(loc='upper left' if raw == 'vas' else 'upper right', frameon=False, fontsize=14)
        ax.tick_params(axis='x', labelbottom=True)
    for ax in axes:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.tick_params(axis='x', rotation=45, labelsize=18)
    axes[-1].set_xlabel('Time', fontsize=18, fontweight='bold')
    fig.tight_layout(h_pad=1.2)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'Figure saved: {output} | {len(df)} observations | along-shelf azimuth: {args.azimuth:g} degrees')


if __name__ == '__main__':
    main()
