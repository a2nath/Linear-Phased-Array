#!/usr/bin/env python3
"""Run the actual PhaseArray executable with every JSON scenario in tests/."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


# The fixtures use legacy names and nested arrays; the CLI accepts CSV rows.
OPTIONS = {
    "timeslot": ("timeslot",),
    "timeslots": ("timeslots",),
    "frequency": ("frequency",),
    "bandwidth": ("bandwidth",),
    "symbolrate": ("symbolrate",),
    "blockpersymbol": ("blockpersymbol", "blockspersym"),
    "height": ("height", "antenna_height"),
    "ms_grx": ("ms_grx", "ms_grx_dB"),
    "system_noise": ("system_noise", "system_noise_dB"),
    "mobile_stations": ("mobile_stations",),
    "base_stations": ("base_stations",),
    "slimit": ("slimit", "sinr_limit_dB"),
    "area": ("area",),
    "base_station_theta_c": ("base_station_theta_c",),
    "base_station_location": ("base_station_location",),
    "mobile_station_location": ("mobile_station_location",),
    "base_station_antenna_counts": ("base_station_antenna_counts",),
    "base_station_power_range_dBm": ("base_station_power_range_dBm",),
    "base_station_scan_angle_range_deg": ("base_station_scan_angle_range_deg",),
    "antenna_spacing": ("antenna_spacing",),
    "antenna_dims": ("antenna_dims",),
    "antenna_txpower": ("antenna_txpower", "base_station_power_dBm_lut"),
    "scan_angle": ("scan_angle", "base_station_scan_alpha_deg_lut"),
    "ms_selection": ("ms_selection", "base_to_mobile_station_id_selection_lut"),
}


def cli_value(value):
    if isinstance(value, list):
        if value and isinstance(value[0], list):
            return " ".join(",".join(str(item) for item in row) for row in value)
        return ",".join(str(item) for item in value)
    return str(value)


def scenario_arguments(scenario, output):
    arguments = []
    timeslots = 1
    for option, names in OPTIONS.items():
        for name in names:
            if name in scenario:
                value = scenario[name]
                arguments.extend([f"--{option}", cli_value(value)])
                if option in ("antenna_txpower", "scan_angle", "ms_selection"):
                    if isinstance(value, list) and value and isinstance(value[0], list):
                        timeslots = max(timeslots, len(value))
                break
    if "timeslots" not in scenario:
        arguments.extend(["--timeslots", str(timeslots)])
    arguments.extend(["--plot", "--output_dir", str(output)])
    return arguments


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--tests-dir", type=Path, default=root / "tests")
    parser.add_argument("--output-dir", type=Path, default=root / "build/test-results")
    parser.add_argument("--timeout", type=int, default=300, help="Seconds per JSON case")
    args = parser.parse_args()
    executable = args.executable.resolve()
    if not executable.is_file():
        parser.error(f"Executable not found: {executable}")
    fixtures = sorted(args.tests_dir.rglob("*.json"))
    if not fixtures:
        parser.error(f"No JSON scenarios found in {args.tests_dir}")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for fixture in fixtures:
        relative = fixture.relative_to(args.tests_dir)
        case_output = output / relative.with_suffix("")
        case_output.mkdir(parents=True, exist_ok=True)
        result = {"scenario": relative.as_posix(), "status": "error", "returncode": None}
        started = time.monotonic()
        try:
            scenario = json.loads(fixture.read_text(encoding="utf-8-sig"))
            shutil.copyfile(fixture, case_output / "input.json")
            if (root / "font").is_dir():
                shutil.copytree(root / "font", case_output / "font", dirs_exist_ok=True)
            command = [str(executable), *scenario_arguments(scenario, case_output)]
            result["command"] = command
            with (case_output / "stdout.log").open("w") as stdout, \
                    (case_output / "stderr.log").open("w") as stderr:
                process = subprocess.Popen(command, cwd=case_output, stdout=stdout, stderr=stderr)
                try:
                    result["returncode"] = process.wait(timeout=args.timeout)
                    result["status"] = "passed" if result["returncode"] == 0 else "failed"
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                    result["status"] = "timeout"
        except (OSError, ValueError, TypeError, AttributeError) as error:
            result["error"] = str(error)
        result["duration_seconds"] = round(time.monotonic() - started, 3)
        (case_output / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        results.append(result)
        print(f"{result['status'].upper()}: {relative} (exit {result['returncode']})", flush=True)

    failed = sum(result["status"] != "passed" for result in results)
    (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    summary = ["## PhaseArray JSON scenarios", "", f"Cases: {len(results)}; failures/errors/timeouts: {failed}.",
               "", "| Scenario | Result | Exit code |", "|---|---|---|"]
    summary.extend(f"| {r['scenario']} | {r['status']} | {r['returncode']} |" for r in results)
    summary_text = "\n".join(summary) + "\n"
    (output / "summary.md").write_text(summary_text, encoding="utf-8")

    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(summary_text)

    if failed:
        print(f"::warning::{failed}/{len(results)} JSON scenarios failed; see test artifacts.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
