"""Automated tests for the NESTbench25 benchmark codes.

Quick tests (default) verify that every C++ engine builds (standalone and
library) and runs, that every Python engine runs, that results land in
report_answer.dat, and that default-protocol results are consistent with the
values reported in the paper (loose tolerances at reduced statistics).

Set NESTBENCH_FULL_TESTS=1 to also run the slower cross-language agreement
checks at higher statistics.

Usage:  pip install pytest
        pytest tests/ -v
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PROBLEMS = ["trap_overdamped", "trap_underdamped", "erasure", "ising", "abp"]
FULL = os.environ.get("NESTBENCH_FULL_TESTS", "0") == "1"

# quick-run arguments and expected default-protocol order parameter
# (value, absolute tolerance) at the reduced statistics used here
QUICK = {
    "trap_overdamped":  (["-n", "500", "-s", "4"],  8.33,  0.6),
    "trap_underdamped": (["-n", "500", "-s", "4"],  2.905, 0.4),
    "erasure":          (["-n", "500", "-s", "4"],  0.478, 0.06),
    "ising":            (["-n", "25",  "-s", "2"],  86.5,  8.0),
    "abp":              (["-n", "4000"],            0.009, 0.006),  # Delta at reduced n (noise-biased upward)
}

PY_QUICK = {
    "trap_overdamped":  (["-n", "2000"], 8.33,  0.6),
    "trap_underdamped": (["-n", "2000"], 2.905, 0.4),
    "erasure":          (["-n", "2000"], 0.478, 0.06),
    "ising":            (["-n", "100"],  86.5,  8.0),
    "abp":              (["-n", "4000"], 0.009, 0.006),
}


def has_torch():
    try:
        import torch  # noqa: F401
        return True
    except Exception:
        return False


def parse_order_parameter(report: Path) -> float:
    """Return the order parameter from the last line of report_answer.dat."""
    text = report.read_text().strip()
    assert text, f"{report} is empty"
    last = text.splitlines()[-1]
    m = re.search(r"order parameter\s*=\s*([-+0-9.eE]+)", last)
    assert m, f"could not parse order parameter from: {last!r}"
    return float(m.group(1))


def parse_report(report: Path):
    """Return (value, uncertainty, n_samples) from the last line of report_answer.dat."""
    last = report.read_text().strip().splitlines()[-1]
    m = re.search(r"order parameter\s*=\s*(\S+)\s*±\s*(\S+),\s*n_samples\s*=\s*(\d+)", last)
    assert m, f"could not parse value ± uncertainty from: {last!r}"
    return float(m.group(1)), float(m.group(2)), int(m.group(3))


def build_sim(cwd):
    if not (cwd / "sim").exists():
        r = run(["make", "standalone"], cwd)
        assert r.returncode == 0, r.stderr


def run(cmd, cwd, timeout=1200):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                          timeout=timeout, env={**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE"})


# ----------------------------------------------------------------------
# C++ builds
# ----------------------------------------------------------------------
@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_builds_standalone(problem):
    cwd = ROOT / problem / "cpp_code"
    r = run(["make", "standalone"], cwd)
    assert r.returncode == 0, r.stderr
    assert (cwd / "sim").exists()


@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_builds_library(problem):
    cwd = ROOT / problem / "cpp_code"
    r = run(["make", "library"], cwd)
    assert r.returncode == 0, r.stderr
    assert list(cwd.glob("libengine_*.a")), "static library not produced"
    # rebuild the standalone binary for the tests that follow
    r = run(["make", "standalone"], cwd)
    assert r.returncode == 0, r.stderr


# ----------------------------------------------------------------------
# C++ default-protocol runs
# ----------------------------------------------------------------------
@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_default_protocol(problem):
    cwd = ROOT / problem / "cpp_code"
    if not (cwd / "sim").exists():
        r = run(["make", "standalone"], cwd)
        assert r.returncode == 0, r.stderr
    args, expected, tol = QUICK[problem]
    r = run(["./sim"] + args, cwd)
    assert r.returncode == 0, r.stderr
    op = parse_order_parameter(cwd / "report_answer.dat")
    assert abs(op - expected) < tol, (
        f"{problem}: order parameter {op} differs from expected {expected} by more than {tol}")


@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_learned_protocol_loads(problem):
    """-p <file> must work without recompilation."""
    cwd = ROOT / problem / "cpp_code"
    learned = cwd / "input_control_parameters_learned.dat"
    if not learned.exists():
        pytest.skip("no learned protocol provided for this problem")
    if not (cwd / "sim").exists():
        r = run(["make", "standalone"], cwd)
        assert r.returncode == 0, r.stderr
    args = QUICK[problem][0]
    r = run(["./sim", "-p", learned.name] + args, cwd)
    assert r.returncode == 0, r.stderr
    op = parse_order_parameter(cwd / "report_answer.dat")
    assert op == op, "order parameter is NaN"  # NaN check


def test_cpp_erasure_learned_protocol_beats_default():
    cwd = ROOT / "erasure" / "cpp_code"
    if not (cwd / "sim").exists():
        r = run(["make", "standalone"], cwd)
        assert r.returncode == 0, r.stderr
    r = run(["./sim", "-p", "input_control_parameters_learned.dat", "-n", "500", "-s", "4"], cwd)
    assert r.returncode == 0, r.stderr
    err = parse_order_parameter(cwd / "report_answer.dat")
    assert err < 0.1, f"learned erasure protocol should give error rate << default (~0.48); got {err}"


def test_cpp_abp_visualization_does_not_crash():
    """Regression test: visualize_protocol() used to segfault with the default
    protocol (degenerate work histogram, W identically zero)."""
    cwd = ROOT / "abp" / "cpp_code"
    if not (cwd / "sim").exists():
        r = run(["make", "standalone"], cwd)
        assert r.returncode == 0, r.stderr
    # run only the C++ part; the python movie step may be missing deps, which is fine
    r = run(["./sim", "-v"], cwd, timeout=3600)
    assert r.returncode in (0, 1, 2), f"visualization crashed (rc={r.returncode}): {r.stderr[-500:]}"
    assert "Segmentation fault" not in r.stderr


def test_cpp_underdamped_velocity_histogram_scale(tmp_path):
    """Regression test: the velocity histogram was written in nm/microsecond
    (thermal width ~0.007) but plotted on an axis spanning -25..25, so the movie
    showed a single spike at zero. It is now written in units of the thermal
    velocity, where the initial (equilibrium) distribution has unit variance."""
    cwd = ROOT / "trap_underdamped" / "cpp_code"
    build_sim(cwd)
    # stub out the movie step so that the histogram files are kept and no
    # plotting dependencies are needed
    stub = tmp_path / "python"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    env = {**os.environ, "PATH": f"{tmp_path}{os.pathsep}{os.environ['PATH']}"}
    r = subprocess.run(["./sim", "-v", "-p", "input_control_parameters_learned.dat"], cwd=cwd,
                       capture_output=True, text=True, timeout=3600, env=env)
    assert r.returncode == 0, r.stderr[-500:]
    rows = [tuple(map(float, line.split()))
            for line in (cwd / "report_velocity_time_0.dat").read_text().splitlines()]
    for f in cwd.glob("report_*_time_*.dat"):
        f.unlink()
    v = [row[0] for row in rows]
    p = [row[1] for row in rows]
    norm = sum(p)
    mean = sum(vi * pi for vi, pi in zip(v, p)) / norm
    std = (sum((vi - mean) ** 2 * pi for vi, pi in zip(v, p)) / norm) ** 0.5
    assert 0.8 < std < 1.2, f"initial velocity distribution should have unit width; got {std}"
    assert 0.3 < max(p) < 0.5, f"peak height should be ~0.4; got {max(p)}"
    movie = (cwd / "movie.py").read_text()
    assert "ax_v.set_xlim(-10, 10)" in movie


# ----------------------------------------------------------------------
# C++ uncertainties and seeding
# ----------------------------------------------------------------------
# small runs with several samples: (arguments, minimum plausible standard error)
UNCERTAINTY = {
    "trap_overdamped":  (["-n", "200", "-s", "8"], 1e-3),
    "trap_underdamped": (["-n", "200", "-s", "8"], 1e-3),
    "erasure":          (["-n", "200", "-s", "8"], 1e-4),
    "ising":            (["-n", "10",  "-s", "6"], 1e-2),
    "abp":              (["-n", "2000", "-s", "8"], 1e-6),
}


@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_uncertainty_is_meaningful(problem):
    """Regression test: the engines used to re-seed from the wall clock for every
    sample, so samples taken within the same second were bit-identical and the
    reported uncertainty was ~0 or nan."""
    cwd = ROOT / problem / "cpp_code"
    build_sim(cwd)
    args, min_err = UNCERTAINTY[problem]
    r = run(["./sim"] + args, cwd)
    assert r.returncode == 0, r.stderr
    value, err, n_samples = parse_report(cwd / "report_answer.dat")
    assert n_samples == int(args[-1])
    assert err == err, "uncertainty is NaN"
    assert err > min_err, (
        f"{problem}: uncertainty {err} is implausibly small; samples are not independent")
    assert err < abs(value) + 1.0, f"{problem}: uncertainty {err} is implausibly large"


@pytest.mark.parametrize("problem", PROBLEMS)
def test_cpp_seed_is_reproducible(problem):
    """-r <seed> gives identical output; different seeds give different output."""
    cwd = ROOT / problem / "cpp_code"
    build_sim(cwd)
    args = UNCERTAINTY[problem][0]
    outputs = []
    for seed in ("11", "11", "12"):
        r = run(["./sim", "-r", seed] + args, cwd)
        assert r.returncode == 0, r.stderr
        outputs.append((cwd / "report_answer.dat").read_text())
    assert outputs[0] == outputs[1], "same seed gave different results"
    assert outputs[0] != outputs[2], "different seeds gave identical results"


# ----------------------------------------------------------------------
# Data files shared by the C++ and Python implementations
# ----------------------------------------------------------------------
@pytest.mark.parametrize("problem", PROBLEMS)
def test_shared_data_files_identical(problem):
    cpp, py = ROOT / problem / "cpp_code", ROOT / problem / "python_code"
    for f in sorted(py.glob("*.dat")):
        if f.name.startswith("report_") or not (cpp / f.name).exists():
            continue
        assert f.read_bytes() == (cpp / f.name).read_bytes(), (
            f"{f.relative_to(ROOT)} differs from the copy in cpp_code")


# ----------------------------------------------------------------------
# Python engines
# ----------------------------------------------------------------------
@pytest.mark.skipif(not has_torch(), reason="torch not installed")
@pytest.mark.parametrize("problem", PROBLEMS)
def test_python_default_protocol(problem):
    cwd = ROOT / problem / "python_code"
    args, expected, tol = PY_QUICK[problem]
    r = run([sys.executable, f"engine_{problem}.py"] + args, cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    op = parse_order_parameter(cwd / "report_answer.dat")
    assert abs(op - expected) < tol, (
        f"{problem} (python): order parameter {op} differs from expected {expected} by more than {tol}")


@pytest.mark.skipif(not has_torch(), reason="torch not installed")
def test_python_abp_learned_protocol():
    cwd = ROOT / "abp" / "python_code"
    r = run([sys.executable, "engine_abp.py", "-p", "input_control_parameters_learned.dat", "-n", "4000"], cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    op = parse_order_parameter(cwd / "report_answer.dat")
    # work extraction (<W> ~ -2.96) pulls the order parameter well below the
    # default-protocol value of ~0.009
    assert op < 0.005, f"learned-protocol order parameter should be ~0 or negative: {op}"


@pytest.mark.skipif(not has_torch(), reason="torch not installed")
@pytest.mark.parametrize("problem", PROBLEMS)
def test_python_seed_is_reproducible(problem):
    """-r <seed> gives identical output; different seeds give different output."""
    cwd = ROOT / problem / "python_code"
    args = PY_QUICK[problem][0]
    outputs = []
    for seed in ("11", "11", "12"):
        r = run([sys.executable, f"engine_{problem}.py", "-r", seed] + args, cwd)
        assert r.returncode == 0, r.stderr[-2000:]
        outputs.append((cwd / "report_answer.dat").read_text())
    assert outputs[0] == outputs[1], "same seed gave different results"
    assert outputs[0] != outputs[2], "different seeds gave identical results"


# ----------------------------------------------------------------------
# Cross-language agreement (slow; enable with NESTBENCH_FULL_TESTS=1)
# ----------------------------------------------------------------------
@pytest.mark.skipif(not FULL, reason="set NESTBENCH_FULL_TESTS=1 to run")
@pytest.mark.skipif(not has_torch(), reason="torch not installed")
@pytest.mark.parametrize("problem,args_c,args_p,tol", [
    ("abp",              ["-n", "200000"], ["-n", "200000"], 0.001),
    ("trap_overdamped",  ["-n", "10000", "-s", "20"], ["-n", "200000"], 0.05),
    ("trap_underdamped", ["-n", "10000", "-s", "20"], ["-n", "200000"], 0.05),
    ("erasure",          ["-n", "10000", "-s", "20"], ["-n", "200000"], 0.01),
])
def test_cpp_python_agreement(problem, args_c, args_p, tol):
    cwd_c = ROOT / problem / "cpp_code"
    cwd_p = ROOT / problem / "python_code"
    if not (cwd_c / "sim").exists():
        r = run(["make", "standalone"], cwd_c)
        assert r.returncode == 0, r.stderr
    r = run(["./sim"] + args_c, cwd_c, timeout=3600)
    assert r.returncode == 0, r.stderr
    op_c = parse_order_parameter(cwd_c / "report_answer.dat")
    r = run([sys.executable, f"engine_{problem}.py"] + args_p, cwd_p, timeout=3600)
    assert r.returncode == 0, r.stderr[-2000:]
    op_p = parse_order_parameter(cwd_p / "report_answer.dat")
    assert abs(op_c - op_p) < tol, (
        f"{problem}: C++ ({op_c}) and Python ({op_p}) disagree by more than {tol}")
