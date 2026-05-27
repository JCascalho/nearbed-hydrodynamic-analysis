"""Near-bed hydrodynamic forcing analysis.

This command-line routine computes current, wave, and combined-flow
near-bed forcing metrics from synchronized current and wave observations.

Calculated outputs include:

- current magnitude and direction;
- current-only shear velocity and shear stress;
- wave orbital velocity and wave shear stress;
- wave-current angle;
- mean and maximum combined bed shear stress;
- combined-flow shear velocity;
- critical shear velocity;
- native-sediment mobilization flags;
- event-scale and diagnostic mixing-layer thickness estimates;
- interval and threshold summary statistics.

The default physical formulations follow classical sediment-transport and
wave-current interaction approaches used in the original research workflow.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ["ve", "vn", "D50", "z", "h", "Tz", "Hs", "dir_w", "Time", "tau_cr"]


@dataclass
class PhysicalConstants:
    """Physical constants and empirical factors."""

    gravity: float = 9.81
    rho_seawater: float = 1025.0
    roughness_factor: float = 2.5
    von_karman: float = 0.41
    wave_period_factor: float = 1.281


@dataclass
class TimeInterval:
    """Named time interval for summary statistics."""

    start: str
    end: str
    label: str


@dataclass
class Config:
    """Runtime settings."""

    constants: PhysicalConstants
    thresholds: tuple[float, ...] = (0.000, 0.259, 0.303)
    sheet_name: Optional[str] = None


DEFAULT_INTERVALS = [
    TimeInterval("2020-08-18 17:50:08", "2020-08-30 23:50:08", "C0_C2"),
    TimeInterval("2020-08-31 00:05:08", "2020-09-28 23:50:08", "C2_C3"),
    TimeInterval("2020-09-29 00:05:08", "2020-11-14 23:50:08", "C3_C4"),
    TimeInterval("2020-11-15 00:05:08", "2021-01-14 23:50:08", "C4_C5"),
    TimeInterval("2021-01-15 00:05:08", "2021-04-22 23:50:08", "C5_C6"),
    TimeInterval("2021-04-23 00:05:08", "2021-10-27 23:50:08", "C6_C7"),
]


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Strip spaces and normalize known tau_cr variants."""
    out = df.copy()
    out.columns = [str(col).strip() for col in out.columns]
    tau_aliases = {"Tau_cr": "tau_cr", "TAU_CR": "tau_cr", "tau_CR": "tau_cr"}
    out = out.rename(columns={old: new for old, new in tau_aliases.items() if old in out.columns})
    return out


def read_input_table(input_path: Path, sheet_name: Optional[str] = None) -> pd.DataFrame:
    """Read CSV or Excel input."""
    suffix = input_path.suffix.lower()
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(input_path)
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(input_path, sheet_name=sheet_name or 0)
    raise ValueError(f"Unsupported input file type: {input_path.suffix}")


def default_output_path(input_path: Path, user_output: Optional[str] = None) -> Path:
    """Create output path when not provided."""
    if user_output:
        return Path(user_output)
    return input_path.with_name(f"OUTPUT_{input_path.stem}_HYDRODYNAMIC_PARAMETERS.xlsx")


def parse_thresholds(text: str) -> tuple[float, ...]:
    """Parse comma-separated tau_max thresholds."""
    values = []
    for item in str(text).split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    if not values:
        raise argparse.ArgumentTypeError("At least one threshold is required.")
    return tuple(values)


def load_intervals(intervals_csv: Optional[str | Path]) -> list[TimeInterval]:
    """Load intervals from CSV or use built-in defaults."""
    if intervals_csv is None:
        return DEFAULT_INTERVALS

    table = pd.read_csv(intervals_csv)
    required = {"start", "end", "label"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError("Intervals CSV must contain columns: start, end, label")

    return [
        TimeInterval(str(row["start"]), str(row["end"]), str(row["label"]))
        for _, row in table.iterrows()
    ]


def validate_input(df: pd.DataFrame) -> pd.DataFrame:
    """Validate and clean required columns."""
    out = normalize_columns(df)
    missing = [col for col in REQUIRED_COLUMNS if col not in out.columns]
    if missing:
        raise ValueError(f"Missing required input columns: {missing}")

    out = out.copy()
    out["Time"] = pd.to_datetime(out["Time"], errors="coerce")
    numeric_cols = [col for col in REQUIRED_COLUMNS if col != "Time"]
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    before = len(out)
    out = out.dropna(subset=REQUIRED_COLUMNS).reset_index(drop=True)
    dropped = before - len(out)
    if out.empty:
        raise ValueError("No valid rows remain after required-column cleaning.")
    if dropped:
        print(f"Warning: dropped {dropped} rows with missing or invalid required values.")
    return out


def current_magnitude(ve, vn):
    return np.hypot(ve, vn)


def current_azimuth(ve, vn):
    return (np.degrees(np.arctan2(ve, vn)) + 360.0) % 360.0


def wave_current_angle(current_dir, wave_dir):
    return np.degrees(
        np.arctan2(
            np.sin(np.radians(wave_dir - current_dir)),
            np.cos(np.radians(wave_dir - current_dir)),
        )
    )


def nikuradse_roughness(d50, constants: PhysicalConstants):
    return constants.roughness_factor * d50


def bed_roughness_length(kn):
    return kn / 30.0


def natural_scaling_period(h, constants: PhysicalConstants):
    return np.sqrt(h / constants.gravity)


def peak_wave_period(tz, constants: PhysicalConstants):
    return constants.wave_period_factor * tz


def wave_correction_factor(t):
    return (6500.0 + (0.56 + 15.54 * t) ** 6.0) ** (1.0 / 6.0)


def orbital_velocity_amplitude(hs, tn, wcf, t):
    return (0.25 * hs) / (tn * (1.0 + wcf * t**2.0) ** 3.0)


def wave_friction_factor(relative_roughness):
    return 0.237 * relative_roughness ** (-0.52)


def process_hydrodynamic_data(df: pd.DataFrame, intervals: list[TimeInterval], constants: PhysicalConstants) -> pd.DataFrame:
    """Compute near-bed hydrodynamic parameters."""
    out = validate_input(df)

    # Currents
    out["u_z"] = current_magnitude(out["ve"], out["vn"])
    out["dir_c"] = current_azimuth(out["ve"], out["vn"])

    # Roughness and current friction factor
    out["kn"] = nikuradse_roughness(out["D50"], constants)
    out["z0"] = bed_roughness_length(out["kn"])
    out["z_over_z0"] = out["z"] / out["z0"]
    out.loc[out["z_over_z0"] <= 1.0, "z_over_z0"] = np.nan

    out["fc"] = 2.0 * (constants.von_karman / np.log(out["z_over_z0"])) ** 2.0
    out["u_star_current"] = out["u_z"] * np.sqrt(out["fc"] / 2.0)
    out["u_star"] = out["u_star_current"]
    out["tau_c"] = 0.5 * constants.rho_seawater * out["fc"] * out["u_z"] ** 2.0

    # Waves
    out["Tn"] = natural_scaling_period(out["h"], constants)
    out["Tp"] = peak_wave_period(out["Tz"], constants)
    out["t"] = out["Tn"] / out["Tz"]
    out["W_cf"] = wave_correction_factor(out["t"])
    out["urms"] = orbital_velocity_amplitude(out["Hs"], out["Tn"], out["W_cf"], out["t"])
    out["uw"] = np.sqrt(2.0) * out["urms"]
    out["A"] = (out["uw"] * out["Tp"]) / (2.0 * np.pi)
    out["r"] = out["A"] / out["kn"]
    out.loc[out["r"] <= 0.0, "r"] = np.nan
    out["fwr"] = wave_friction_factor(out["r"])
    out["tau_w"] = 0.5 * constants.rho_seawater * out["fwr"] * out["uw"] ** 2.0

    # Combined wave-current forcing
    out["angle_phi"] = wave_current_angle(out["dir_c"], out["dir_w"])
    ratio = out["tau_w"] / (out["tau_c"] + out["tau_w"])
    out["tau_m"] = out["tau_c"] * (1.0 + 1.2 * ratio**3.2)
    phi_rad = np.radians(out["angle_phi"])
    out["tau_max"] = np.hypot(
        out["tau_m"] + out["tau_w"] * np.cos(phi_rad),
        out["tau_w"] * np.sin(phi_rad),
    )

    out["u_star_max"] = np.sqrt(out["tau_max"] / constants.rho_seawater)
    out["u_star_cr"] = np.sqrt(out["tau_cr"] / constants.rho_seawater)
    out["tau_excess"] = np.maximum(out["tau_max"] - out["tau_cr"], 0.0)
    out["is_mobilized_native"] = out["tau_max"] > out["tau_cr"]

    # Interval-level tau_max references for mixing-layer estimates
    out["interval_label"] = pd.NA
    out["tau_max_interval_max"] = np.nan
    out["tau_max_interval_mean"] = np.nan
    out["tau_max_interval_p95"] = np.nan

    for interval in intervals:
        start = pd.to_datetime(interval.start)
        end = pd.to_datetime(interval.end)
        mask = (out["Time"] >= start) & (out["Time"] <= end)
        if not mask.any():
            continue
        tau_interval = pd.to_numeric(out.loc[mask, "tau_max"], errors="coerce").dropna()
        if tau_interval.empty:
            continue
        out.loc[mask, "interval_label"] = interval.label
        out.loc[mask, "tau_max_interval_max"] = tau_interval.max()
        out.loc[mask, "tau_max_interval_mean"] = tau_interval.mean()
        out.loc[mask, "tau_max_interval_p95"] = tau_interval.quantile(0.95)

    out["delta_mix"] = 0.07 * np.maximum(out["tau_max_interval_max"] - out["tau_cr"], 0.0) + 6.0 * out["D50"]
    out["delta_mix_mean_tau"] = 0.07 * np.maximum(out["tau_max_interval_mean"] - out["tau_cr"], 0.0) + 6.0 * out["D50"]
    out["delta_mix_p95_tau"] = 0.07 * np.maximum(out["tau_max_interval_p95"] - out["tau_cr"], 0.0) + 6.0 * out["D50"]
    out["tau_max_interval"] = out["tau_max_interval_max"]

    return out


def calculate_interval_statistics(
    processed: pd.DataFrame,
    intervals: list[TimeInterval],
    thresholds: tuple[float, ...],
) -> dict[str, pd.DataFrame]:
    """Compute interval and threshold statistics."""
    interval_rows = []
    threshold_rows = []

    for interval in intervals:
        start = pd.to_datetime(interval.start)
        end = pd.to_datetime(interval.end)
        subset = processed.loc[(processed["Time"] >= start) & (processed["Time"] <= end)]
        if subset.empty:
            continue

        interval_rows.append({
            "Interval": interval.label,
            "Start": start,
            "End": end,
            "Count": len(subset),
            "Mean_tau_max": subset["tau_max"].mean(),
            "P95_tau_max": subset["tau_max"].quantile(0.95),
            "Max_tau_max": subset["tau_max"].max(),
            "Min_tau_max": subset["tau_max"].min(),
            "Mean_tau_excess": subset["tau_excess"].mean(),
            "P95_tau_excess": subset["tau_excess"].quantile(0.95),
            "Max_tau_excess": subset["tau_excess"].max(),
            "Mobilized_fraction_percent": 100.0 * subset["is_mobilized_native"].mean(),
            "Mean_delta_mix_max_tau": subset["delta_mix"].mean(),
            "Mean_delta_mix_p95_tau": subset["delta_mix_p95_tau"].mean(),
            "Mean_delta_mix_mean_tau": subset["delta_mix_mean_tau"].mean(),
            "Mean_u_star_current": subset["u_star_current"].mean(),
            "Max_u_star_current": subset["u_star_current"].max(),
            "Min_u_star_current": subset["u_star_current"].min(),
            "Mean_u_star_max": subset["u_star_max"].mean(),
            "P95_u_star_max": subset["u_star_max"].quantile(0.95),
            "Max_u_star_max": subset["u_star_max"].max(),
            "Mean_u_star_cr": subset["u_star_cr"].mean(),
            "Mobilized_fraction_u_star_percent": 100.0 * (subset["u_star_max"] > subset["u_star_cr"]).mean(),
        })

        for threshold in thresholds:
            selected = subset.loc[subset["tau_max"] > threshold]
            if selected.empty:
                continue
            threshold_rows.append({
                "Interval": interval.label,
                "Threshold_tau_max": threshold,
                "Count": len(selected),
                "Max_tau_max": selected["tau_max"].max(),
                "Mean_tau_max_over_threshold": selected["tau_max"].mean(),
                "Mean_tau_excess_over_threshold": selected["tau_excess"].mean(),
                "Mean_u_star_current_over_threshold": selected["u_star_current"].mean(),
                "Mean_u_star_max_over_threshold": selected["u_star_max"].mean(),
                "Mean_u_star_cr_over_threshold": selected["u_star_cr"].mean(),
            })

    stats = {}
    if interval_rows:
        stats["interval_stats"] = pd.DataFrame(interval_rows)
    if threshold_rows:
        stats["threshold_stats"] = pd.DataFrame(threshold_rows)
    return stats


def safe_sheet_name(name: str, used: set[str], max_len: int = 31) -> str:
    """Create a valid unique Excel sheet name."""
    invalid = {":": "-", "\\": "-", "/": "-", "?": "", "*": "", "[": "(", "]": ")"}
    out = str(name)
    for bad, good in invalid.items():
        out = out.replace(bad, good)
    out = out.strip().strip("'") or "Sheet"
    base = out[:max_len]
    candidate = base
    index = 1
    while candidate in used:
        suffix = f"_{index}"
        candidate = base[: max_len - len(suffix)] + suffix
        index += 1
    used.add(candidate)
    return candidate


def metadata_table(cfg: Config, intervals: list[TimeInterval]) -> pd.DataFrame:
    """Build metadata table."""
    constants = cfg.constants
    return pd.DataFrame({
        "parameter": [
            "routine",
            "gravity",
            "rho_seawater",
            "roughness_factor",
            "von_karman",
            "wave_period_factor",
            "thresholds_tau_max",
            "n_intervals",
            "required_columns",
            "delta_mix_formula",
        ],
        "value": [
            "nearbed_hydrodynamic_analysis",
            constants.gravity,
            constants.rho_seawater,
            constants.roughness_factor,
            constants.von_karman,
            constants.wave_period_factor,
            ",".join(f"{x:g}" for x in cfg.thresholds),
            len(intervals),
            ",".join(REQUIRED_COLUMNS),
            "0.07 * max(tau_max_interval_reference - tau_cr, 0) + 6 * D50",
        ],
    })


def save_results(
    output_path: Path,
    processed: pd.DataFrame,
    stats: dict[str, pd.DataFrame],
    intervals: list[TimeInterval],
    cfg: Config,
) -> None:
    """Export computed results to Excel."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    used = set()
    sheet_map = []

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        sheet = safe_sheet_name("All_Data", used)
        processed.to_excel(writer, sheet_name=sheet, index=False)
        sheet_map.append({"requested_sheet_name": "All_Data", "excel_sheet_name": sheet})

        if "interval_stats" in stats:
            sheet = safe_sheet_name("tau_max_summary", used)
            stats["interval_stats"].to_excel(writer, sheet_name=sheet, index=False)
            sheet_map.append({"requested_sheet_name": "tau_max_summary", "excel_sheet_name": sheet})

        if "threshold_stats" in stats:
            sheet = safe_sheet_name("tau_max_thresholds", used)
            stats["threshold_stats"].to_excel(writer, sheet_name=sheet, index=False)
            sheet_map.append({"requested_sheet_name": "tau_max_thresholds", "excel_sheet_name": sheet})

        for interval in intervals:
            start = pd.to_datetime(interval.start)
            end = pd.to_datetime(interval.end)
            subset = processed.loc[(processed["Time"] >= start) & (processed["Time"] <= end)]
            if subset.empty:
                continue

            sheet = safe_sheet_name(interval.label, used)
            subset.to_excel(writer, sheet_name=sheet, index=False)
            sheet_map.append({"requested_sheet_name": interval.label, "excel_sheet_name": sheet})

            for threshold in cfg.thresholds:
                selected = subset.loc[subset["tau_max"] > threshold]
                if selected.empty:
                    continue
                requested = f"{interval.label}_tau_max_gt_{threshold:.3f}"
                sheet = safe_sheet_name(requested, used)
                selected.to_excel(writer, sheet_name=sheet, index=False)
                sheet_map.append({"requested_sheet_name": requested, "excel_sheet_name": sheet})

        sheet = safe_sheet_name("metadata", used)
        metadata_table(cfg, intervals).to_excel(writer, sheet_name=sheet, index=False)
        sheet_map.append({"requested_sheet_name": "metadata", "excel_sheet_name": sheet})

        sheet = safe_sheet_name("SheetName_Map", used)
        pd.DataFrame(sheet_map).to_excel(writer, sheet_name=sheet, index=False)


def process_file(input_path: Path, output_path: Path, cfg: Config, intervals: list[TimeInterval]) -> Path:
    """Run the full workflow."""
    print("=== Near-bed hydrodynamic forcing analysis ===")
    print(f"Input: {input_path}")

    raw = read_input_table(input_path, cfg.sheet_name)
    processed = process_hydrodynamic_data(raw, intervals, cfg.constants)
    stats = calculate_interval_statistics(processed, intervals, cfg.thresholds)
    save_results(output_path, processed, stats, intervals, cfg)

    print(f"Output saved: {output_path}")
    print("Processing completed successfully.")
    return output_path


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description="Compute near-bed current, wave, and combined-flow hydrodynamic forcing metrics."
    )
    parser.add_argument("--input", "-i", required=True, type=Path, help="Input CSV or Excel file.")
    parser.add_argument("--output", "-o", default=None, help="Output Excel workbook.")
    parser.add_argument("--sheet", default=None, help="Excel sheet name. Ignored for CSV input.")
    parser.add_argument("--intervals-csv", default=None, help="Optional CSV with columns start,end,label.")
    parser.add_argument(
        "--thresholds",
        type=parse_thresholds,
        default=(0.000, 0.259, 0.303),
        help="Comma-separated tau_max thresholds, e.g. 0,0.259,0.303.",
    )
    parser.add_argument("--gravity", type=float, default=9.81)
    parser.add_argument("--rho-seawater", type=float, default=1025.0)
    parser.add_argument("--roughness-factor", type=float, default=2.5)
    parser.add_argument("--von-karman", type=float, default=0.41)
    parser.add_argument("--wave-period-factor", type=float, default=1.281)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    constants = PhysicalConstants(
        gravity=args.gravity,
        rho_seawater=args.rho_seawater,
        roughness_factor=args.roughness_factor,
        von_karman=args.von_karman,
        wave_period_factor=args.wave_period_factor,
    )
    cfg = Config(constants=constants, thresholds=args.thresholds, sheet_name=args.sheet)
    intervals = load_intervals(args.intervals_csv)
    output_path = default_output_path(args.input, args.output)
    process_file(args.input, output_path, cfg, intervals)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
