
"""
===============================================================================
M4 — Near-Bed Hydrodynamics and Tracer Mobility
File: M4_nearbed_hydrodynamics_tracer_mobility.py
===============================================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Original date: November 2025
Revision date: October 2026
M4 study-module documentation adaptation: November 2026

Purpose
-------
This script documents the Python-based computational workflow developed to
derive hydrodynamic parameters, including current-induced and wave-induced
bed shear stresses, current-only and combined-flow shear velocities,
wave–current interaction metrics, active mixing-layer thickness (delta_mix),
and statistical descriptors for each sampling interval.

The analysis is based on synchronized current and wave records collected
during the monitoring period (August 2020 to October 2021).

Scientific Basis
----------------
The workflow follows established formulations for sediment mobility and
wave–current interaction, including those of Bagnold (1956),
Soulsby and Smallman (1986), and Soulsby and Clarke (2005).

These formulations are used to calculate:
    tau_c   : current-induced bed shear stress
    tau_w   : wave-induced bed shear stress
    tau_m   : mean combined bed shear stress
    tau_max : maximum combined wave–current bed shear stress

Wave-Orbital Parameters
-----------------------
Near-bed root-mean-square wave-orbital velocity (urms) is estimated from
significant wave height (Hs), mean zero-crossing period (Tz), and water
depth (h), following Soulsby and Smallman (1986).

The equivalent orbital velocity amplitude is calculated as:

    uw = sqrt(2) * urms

An equivalent wave period is estimated as:

    T_equiv = 1.281 * Tz

This equivalent period, rather than a directly measured spectral peak
period (Tp), is used to calculate the near-bed orbital excursion:

    A = uw * T_equiv / (2 * pi)

The equivalent-period approach follows the formulation described by
Soulsby (2006).

Active Mixing-Layer Thickness
-----------------------------
The active mixing-layer thickness (delta_mix) is estimated using the
parameterization of Harris and Wiberg (1997):

    delta_mix = 0.07 * max(tau_stat - tau_cr, 0) + 6 * D50

where:
    tau_stat : selected interval statistic of tau_max
    tau_cr   : critical bed shear stress for native sediment
    D50      : median grain size of the native sediment

Three estimates are calculated for each sampling interval:

    delta_mix          : based on the absolute interval maximum of tau_max
    delta_mix_p95_tau  : based on the 95th percentile of tau_max
    delta_mix_mean_tau : based on the interval mean of tau_max

The 95th-percentile estimate (delta_mix_p95_tau) is the quantity adopted
in the tracer-derived sediment transport calculations described in the
manuscript. The other estimates are retained for comparison and
sensitivity assessment.

When the selected shear-stress statistic does not exceed tau_cr,
the estimated mixing-layer thickness reduces to the baseline value
of 6 * D50.

Critical Shear-Stress Thresholds
--------------------------------
The native-sediment critical shear stress (tau_cr), based on the native
sediment D50, is used to evaluate sediment mobility and calculate
mixing-layer thickness.

The separate tracer threshold (0.259 N m^-2) is explicitly included in
tau_max_thresholds for exceedance summaries and threshold-specific subsets.
The native-sediment threshold (nominally 0.303 N m^-2) is supplied as
tau_cr in the input workbook. The native threshold, not the tracer
threshold, is used to calculate tau_excess, is_mobilized_native and delta_mix.
Exceedance row counts are not counts of independent entrainment events.

Input Data Assumptions
----------------------
The input Excel file must contain synchronized, burst-averaged current
and wave measurements with the following mandatory fields:

    ve, vn : horizontal current velocity components (m s^-1)
    D50    : median native sediment grain size (m)
    z      : current measurement elevation above the seabed (m)
    h      : water depth (m)
    Tz     : mean zero-crossing wave period (s)
    Hs     : significant wave height (m)
    dir_w  : wave direction (degrees)
    Time   : timestamp (datetime format)
    tau_cr : critical bed shear stress for native sediment (N m^-2)

Records lacking mandatory input values are excluded from processing
when the corresponding calculations cannot be performed.

Implementation Details
----------------------
- Language: Python 3.10+
- Libraries: pandas, numpy, openpyxl, argparse, logging
- Input: Excel workbook containing current and wave measurements
- Output: Excel workbook containing calculated hydrodynamic parameters,
          mixing-layer thickness estimates, and interval statistics
- Structure: Modular functions for physical calculations, data processing,
             statistical aggregation, and export routines
"""



# =============================================================================
# 1. CONFIGURATION
# =============================================================================

import pandas as pd
import numpy as np
import logging
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Dict

# -- Physical constants -------------------------------------------------------
@dataclass
class OceanographicConstants:
    GRAVITY: float = 9.81
    RHO_SEAWATER: float = 1025.0
    ROUGHNESS_FACTOR: float = 2.5
    VON_KARMAN: float = 0.41
    WAVE_PERIOD_FACTOR: float = 1.281  # Equivalent wave period / Tz ratio (not measured spectral Tp)

# -- Interval configuration ---------------------------------------------------
@dataclass
class TimeInterval:
    start_date: str
    end_date: str
    sheet_name: str

# -- Processing configuration -------------------------------------------------
@dataclass
class ProcessingConfig:
    input_file: str = 'INPUT_CURRENTS_WAVES_DATA_BASE.xlsx'
    output_file: str = 'OUTPUT_M4_NEARBED_HYDRODYNAMICS_TRACER_MOBILITY.xlsx'
    tau_max_thresholds: Tuple[float, ...] = (0.000, 0.259, 0.303)  # 0.259 tracer; 0.303 native

# Global definitions
CONSTANTS = OceanographicConstants()
CONFIG = ProcessingConfig()

TIME_INTERVALS = [
    TimeInterval('2020-08-18 17:50:08', '2020-08-30 23:50:08', 'C0_C2'),
    TimeInterval('2020-08-31 00:05:08', '2020-09-28 23:50:08', 'C2_C3'),
    TimeInterval('2020-09-29 00:05:08', '2020-11-14 23:50:08', 'C3_C4'),
    TimeInterval('2020-11-15 00:05:08', '2021-01-14 23:50:08', 'C4_C5'),
    TimeInterval('2021-01-15 00:05:08', '2021-04-22 23:50:08', 'C5_C6'),
    TimeInterval('2021-04-23 00:05:08', '2021-10-27 23:50:08', 'C6_C7')
]

# Note: we use 'tau_cr' (all lower-case) inside the code.
# If the Excel uses 'Tau_cr' or 'TAU_CR', we normalize later in main().
REQUIRED_COLUMNS = ['ve', 'vn', 'D50', 'z', 'h', 'Tz', 'Hs', 'dir_w', 'Time', 'tau_cr']


# =============================================================================
# 2. PHYSICAL AND OCEANOGRAPHIC FORMULATIONS
# =============================================================================

def calculate_current_magnitude(ve, vn):
    return np.hypot(ve, vn)

def calculate_current_azimuth(ve, vn):
    return (np.degrees(np.arctan2(ve, vn)) + 360) % 360

def calculate_wave_current_angle(cur_dir, wave_dir):
    return np.degrees(np.arctan2(np.sin(np.radians(wave_dir - cur_dir)),
                                 np.cos(np.radians(wave_dir - cur_dir))))

def calculate_nikuradse_roughness(D50):
    return CONSTANTS.ROUGHNESS_FACTOR * D50

def calculate_bed_roughness_length(kn):
    return kn / 30.0

def calculate_natural_scaling_period(h):
    return np.sqrt(h / CONSTANTS.GRAVITY)

def calculate_equivalent_wave_period(Tz):
    return CONSTANTS.WAVE_PERIOD_FACTOR * Tz

def calculate_period_ratio(Tn, Tz):
    return Tn / Tz

def calculate_wave_correction_factor(t):
    return (6500 + (0.56 + 15.54 * t) ** 6) ** (1 / 6)

def calculate_orbital_velocity_amplitude(Hs, Tn, Wcf, t):
    return (0.25 * Hs) / (Tn * (1 + Wcf * t ** 2) ** 3)

def calculate_max_orbital_velocity(urms):
    return np.sqrt(2) * urms

def calculate_semi_orbital_excursion(uw, T_equiv):
    return (uw * T_equiv) / (2 * np.pi)

def calculate_relative_roughness(A, kn):
    return A / kn

def calculate_current_friction_factor(z, z0):
    # Soulsby/Prandtl-von Kármán log-law based bulk friction factor
    return 2 * (CONSTANTS.VON_KARMAN / np.log(z / z0)) ** 2

def calculate_shear_velocity(fc, uz):
    # u_star = u_z * sqrt(fc/2)
    return uz * np.sqrt(fc / 2)

def calculate_current_shear_stress(fc, uz):
    return 0.5 * CONSTANTS.RHO_SEAWATER * fc * uz ** 2

def calculate_wave_friction_factor(r):
    return 0.237 * r ** (-0.52)

def calculate_wave_shear_stress(fwr, uw):
    return 0.5 * CONSTANTS.RHO_SEAWATER * fwr * uw ** 2

def calculate_mean_bed_shear_stress(tau_c, tau_w):
    ratio = tau_w / (tau_c + tau_w)
    return tau_c * (1 + 1.2 * ratio ** 3.2)

def calculate_max_bed_shear_stress(tau_m, tau_w, phi_deg):
    phi = np.radians(phi_deg)
    return np.hypot(tau_m + tau_w * np.cos(phi), tau_w * np.sin(phi))


# =============================================================================
# 3. DATA PROCESSING WORKFLOW
# =============================================================================

def process_oceanographic_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df['Time'] = pd.to_datetime(df['Time'], errors='coerce')

    # Currents
    df['u_z'] = calculate_current_magnitude(df['ve'], df['vn'])
    df['dir_c'] = calculate_current_azimuth(df['ve'], df['vn'])

    # Roughness and current friction factor
    df['kn'] = df['D50'].apply(calculate_nikuradse_roughness)
    df['z0'] = df['kn'].apply(calculate_bed_roughness_length)
    df['fc'] = df.apply(lambda r: calculate_current_friction_factor(r['z'], r['z0']), axis=1)

    # Current-only shear velocity and current shear stress
    # IMPORTANT: this is NOT equivalent to tau_max. It represents current-only forcing.
    df['u_star_current'] = df.apply(lambda r: calculate_shear_velocity(r['fc'], r['u_z']), axis=1)

    # Backward-compatible alias for older plots/tables.
    # Use u_star_current explicitly in new figures to avoid confusion.
    df['u_star'] = df['u_star_current']

    # Current shear stress
    df['tau_c'] = df.apply(lambda r: calculate_current_shear_stress(r['fc'], r['u_z']), axis=1)

    # Waves
    df['Tn'] = df['h'].apply(calculate_natural_scaling_period)
    df['T_equiv'] = df['Tz'].apply(calculate_equivalent_wave_period)
    df['t']  = df.apply(lambda r: calculate_period_ratio(r['Tn'], r['Tz']), axis=1)
    df['W_cf'] = df['t'].apply(calculate_wave_correction_factor)
    df['urms'] = df.apply(lambda r: calculate_orbital_velocity_amplitude(r['Hs'], r['Tn'], r['W_cf'], r['t']), axis=1)
    df['uw']   = df['urms'].apply(calculate_max_orbital_velocity)
    df['A']    = df.apply(lambda r: calculate_semi_orbital_excursion(r['uw'], r['T_equiv']), axis=1)
    df['r']    = df.apply(lambda r: calculate_relative_roughness(r['A'], r['kn']), axis=1)
    df['fwr']  = df['r'].apply(calculate_wave_friction_factor)
    df['tau_w'] = df.apply(lambda r: calculate_wave_shear_stress(r['fwr'], r['uw']), axis=1)

    # Coupling and maxima
    df['angle_phi'] = df.apply(lambda r: calculate_wave_current_angle(r['dir_c'], r['dir_w']), axis=1)
    df['tau_m'] = df.apply(lambda r: calculate_mean_bed_shear_stress(r['tau_c'], r['tau_w']), axis=1)
    df['tau_max'] = df.apply(lambda r: calculate_max_bed_shear_stress(r['tau_m'], r['tau_w'], r['angle_phi']), axis=1)

    # Equivalent combined-flow shear velocity from tau_max.
    # This is the correct variable to compare against u_star_cr when assessing
    # threshold exceedance under combined wave-current forcing.
    df['u_star_max'] = np.sqrt(df['tau_max'] / CONSTANTS.RHO_SEAWATER)

    # Critical shear velocity corresponding to the tau_cr used in this file.
    # This guarantees exact consistency with tau_max > tau_cr:
    # tau_max > tau_cr  <=>  u_star_max > u_star_cr
    df['u_star_cr'] = np.sqrt(df['tau_cr'] / CONSTANTS.RHO_SEAWATER)
    
    # -----------------------------------------------------------------
    # Harris & Wiberg (1997) active mixing-layer thickness delta_mix
    # -----------------------------------------------------------------
    # Three mixing-layer thickness estimates are calculated using the
    # maximum, mean, and 95th percentile of tau_max within each
    # sampling interval.
    #
    # The 95th-percentile estimate (delta_mix_p95_tau) is adopted
    # for the tracer-derived sediment transport calculations.
    # This reduces the influence of isolated extreme values.
    #
    # The absolute-maximum estimate (delta_mix) and the mean-based
    # estimate (delta_mix_mean_tau) are retained for comparison
    # and sensitivity assessment.
    # -----------------------------------------------------------------

    df['tau_max_interval_max'] = np.nan
    df['tau_max_interval_mean'] = np.nan
    df['tau_max_interval_p95'] = np.nan

    # Time-resolved excess shear above the native-sediment critical threshold.
    df['tau_excess'] = np.maximum(df['tau_max'] - df['tau_cr'], 0.0)
    df['is_mobilized_native'] = df['tau_max'] > df['tau_cr']

    for interval in TIME_INTERVALS:
        mask = (
            (df['Time'] >= pd.to_datetime(interval.start_date)) &
            (df['Time'] <= pd.to_datetime(interval.end_date))
        )
        if not mask.any():
            continue

        tau_interval = pd.to_numeric(df.loc[mask, 'tau_max'], errors='coerce').dropna()
        if tau_interval.empty:
            continue

        df.loc[mask, 'tau_max_interval_max'] = tau_interval.max()
        df.loc[mask, 'tau_max_interval_mean'] = tau_interval.mean()
        df.loc[mask, 'tau_max_interval_p95'] = tau_interval.quantile(0.95)

    # Primary event-scale delta_mix: maximum tau_max within each sampling interval.
    df['delta_mix'] = (
        0.07 * np.maximum(df['tau_max_interval_max'] - df['tau_cr'], 0.0)
        + 6.0 * df['D50']
    )

    # Diagnostic alternatives for sensitivity analysis.
    df['delta_mix_mean_tau'] = (
        0.07 * np.maximum(df['tau_max_interval_mean'] - df['tau_cr'], 0.0)
        + 6.0 * df['D50']
    )
    df['delta_mix_p95_tau'] = (
        0.07 * np.maximum(df['tau_max_interval_p95'] - df['tau_cr'], 0.0)
        + 6.0 * df['D50']
    )

    # Backward-compatible alias for older downstream routines that expected the
    # column name 'tau_max_interval'. It remains the interval maximum.
    df['tau_max_interval'] = df['tau_max_interval_max']

    return df


# =============================================================================
# 4. STATISTICAL ANALYSIS BY INTERVAL
# =============================================================================

def calculate_interval_statistics(df: pd.DataFrame,
                                  intervals: List[TimeInterval],
                                  thresholds: Tuple[float, ...]) -> Dict:
    """Compute τ_max statistics per time interval and threshold."""
    stats = {'interval_stats': [], 'threshold_stats': []}
    for interval in intervals:
        mask = (
            (df['Time'] >= pd.to_datetime(interval.start_date)) &
            (df['Time'] <= pd.to_datetime(interval.end_date))
        )
        subset = df.loc[mask]
        if subset.empty:
            continue

        row = {
            'Interval': interval.sheet_name,
            'Mean_tau_max': subset['tau_max'].mean(),
            'P95_tau_max': subset['tau_max'].quantile(0.95),
            'Max_tau_max': subset['tau_max'].max(),
            'Min_tau_max': subset['tau_max'].min(),
            'Mean_tau_excess': subset['tau_excess'].mean(),
            'P95_tau_excess': subset['tau_excess'].quantile(0.95),
            'Max_tau_excess': subset['tau_excess'].max(),
            'Mobilized_fraction_percent': 100.0 * subset['is_mobilized_native'].mean(),
            'Mean_delta_mix_max_tau': subset['delta_mix'].mean(),
            'Mean_delta_mix_p95_tau': subset['delta_mix_p95_tau'].mean(),
            'Mean_delta_mix_mean_tau': subset['delta_mix_mean_tau'].mean(),
            'Mean_u_star_current': subset['u_star_current'].mean(),
            'Max_u_star_current': subset['u_star_current'].max(),
            'Min_u_star_current': subset['u_star_current'].min(),
            'Mean_u_star_max': subset['u_star_max'].mean(),
            'P95_u_star_max': subset['u_star_max'].quantile(0.95),
            'Max_u_star_max': subset['u_star_max'].max(),
            'Mean_u_star_cr': subset['u_star_cr'].mean(),
            'Mobilized_fraction_u_star_percent': 100.0 * (subset['u_star_max'] > subset['u_star_cr']).mean(),
            'Count': len(subset)
        }
        stats['interval_stats'].append(row)

        for thr in thresholds:
            sel = subset[subset['tau_max'] > thr]
            if not sel.empty:
                stats['threshold_stats'].append({
                    'Interval': f"{interval.sheet_name}_tau_max>{thr:.3f}",
                    'Threshold': thr,
                    'Count': len(sel),
                    'Max_tau_max': sel['tau_max'].max(),
                    'Mean_tau_max_over_thr': sel['tau_max'].mean(),
                    'Mean_tau_excess_over_thr': sel['tau_excess'].mean(),
                    'Mean_u_star_current_over_thr': sel['u_star_current'].mean(),
                    'Mean_u_star_max_over_thr': sel['u_star_max'].mean(),
                    'Mean_u_star_cr_over_thr': sel['u_star_cr'].mean()
                })

    return {k: pd.DataFrame(v) for k, v in stats.items() if v}


# =============================================================================
# 5. OUTPUT AND REPRODUCIBILITY
# =============================================================================

def save_results(output: str,
                 processed: pd.DataFrame,
                 stats: Dict,
                 intervals: List[TimeInterval],
                 thresholds: Tuple[float, ...]):
    """Export all computed results to a structured Excel workbook."""
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # All data
        processed.to_excel(writer, sheet_name='C0_C7', index=False)

        # Interval statistics
        if 'interval_stats' in stats:
            stats['interval_stats'].to_excel(writer, sheet_name='tau_max_summary', index=False)
        if 'threshold_stats' in stats:
            stats['threshold_stats'].to_excel(writer, sheet_name='tau_max_thresholds', index=False)

        # Per-interval sheets and thresholded subsets
        for interval in intervals:
            mask = (
                (processed['Time'] >= pd.to_datetime(interval.start_date)) &
                (processed['Time'] <= pd.to_datetime(interval.end_date))
            )
            subset = processed.loc[mask]
            if not subset.empty:
                subset.to_excel(writer, sheet_name=interval.sheet_name[:31], index=False)
                for thr in thresholds:
                    sel = subset[subset['tau_max'] > thr]
                    if not sel.empty:
                        sheet_name = f"{interval.sheet_name}_tau_max>{thr:.3f}"[:31]
                        sel.to_excel(writer, sheet_name=sheet_name, index=False)


# =============================================================================
# 6. MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compute shear stresses, wave–current metrics, u_star_current, u_star_max, u_star_cr, and δ_mix."
    )
    parser.add_argument('-i', '--input', default=CONFIG.input_file)
    parser.add_argument('-o', '--output', default=CONFIG.output_file)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    df = pd.read_excel(args.input)

    # Normalize possible tau_cr column names from Excel
    df = df.rename(columns={'Tau_cr': 'tau_cr', 'TAU_CR': 'tau_cr'})

    # Basic input integrity check
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required input columns: {missing}")

    processed = process_oceanographic_data(df)
    stats = calculate_interval_statistics(processed, TIME_INTERVALS, CONFIG.tau_max_thresholds)
    save_results(args.output, processed, stats, TIME_INTERVALS, CONFIG.tau_max_thresholds)

    print(f"Processing complete. Results saved to {args.output}")


if __name__ == "__main__":
    main()
