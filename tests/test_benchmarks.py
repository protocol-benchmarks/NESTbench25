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
    text = (cwd / "report_answer.dat").read_text()
    m = re.search(r"mean_work\s*=\s*([-+0-9.eE]+)", text)
    assert m, text
    assert abs(float(m.group(1)) - (-2.96)) < 0.15, f"learned-protocol work should be ~ -2.96: {text}"


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
