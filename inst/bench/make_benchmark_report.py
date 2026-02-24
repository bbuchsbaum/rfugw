#!/usr/bin/env python3
"""Generate a Markdown report from rfugw/POT/thread benchmark CSV files."""

from __future__ import annotations

import csv
import datetime as dt
import math
import sys
from pathlib import Path


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def fmt(x: float, nd: int = 3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    return f"{x:.{nd}f}"


def to_float(d: dict, key: str):
    try:
        return float(d[key])
    except Exception:
        return float("nan")


def md_table(headers, rows):
    out = []
    out.append("| " + " | ".join(headers) + " |")
    out.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


def section_rfugw(rows):
    if not rows:
        return "## rfugw\n\nNo rfugw benchmark CSV found.\n"

    fgw = [r for r in rows if r["suite"] == "fgw"]
    fugw = [r for r in rows if r["suite"] == "fugw"]

    lines = ["## rfugw\n"]
    if fgw:
      lines.append("### FGW\n")
      trows = []
      for r in sorted(fgw, key=lambda x: (int(x["n"]), x["method"])):
          trows.append(
              [
                  r["n"],
                  r["method"],
                  fmt(to_float(r, "median_ms")),
                  fmt(to_float(r, "iter_per_sec")),
                  str(int(float(r["mem_bytes"]))),
              ]
          )
      lines.append(md_table(["n", "method", "median_ms", "iter_per_sec", "mem_bytes"], trows))
      lines.append("")

    if fugw:
      lines.append("### FUGW\n")
      trows = []
      for r in sorted(fugw, key=lambda x: (int(x["n"]), x["method"])):
          trows.append(
              [
                  r["n"],
                  r["method"],
                  fmt(to_float(r, "median_ms")),
                  fmt(to_float(r, "iter_per_sec")),
                  str(int(float(r["mem_bytes"]))),
                  fmt(to_float(r, "warm_inner_total"), 1),
                  fmt(to_float(r, "warm_inner_feat_mean")),
                  fmt(to_float(r, "warm_inner_samp_mean")),
                  fmt(to_float(r, "warm_inner_warm_frac_feat")),
                  fmt(to_float(r, "warm_inner_warm_frac_samp")),
                  fmt(to_float(r, "warm_inner_fallback_frac_feat")),
                  fmt(to_float(r, "warm_inner_fallback_frac_samp")),
              ]
          )
      lines.append(
          md_table(
              [
                  "n",
                  "method",
                  "median_ms",
                  "iter_per_sec",
                  "mem_bytes",
                  "warm_inner_total",
                  "warm_inner_feat_mean",
                  "warm_inner_samp_mean",
                  "warm_frac_feat",
                  "warm_frac_samp",
                  "fallback_frac_feat",
                  "fallback_frac_samp",
              ],
              trows,
          )
      )
      lines.append("")

    return "\n".join(lines)


def section_pot_compare(rfugw_rows, pot_rows):
    if not rfugw_rows or not pot_rows:
        return "## POT vs rfugw\n\nNo comparable POT/rfugw data available.\n"

    lines = ["## POT vs rfugw\n"]
    comp_rows = []

    rf_idx = {(r["suite"], int(r["n"]), r["method"]): r for r in rfugw_rows}
    pot_idx = {(r["suite"], int(r["n"]), r["method"]): r for r in pot_rows}

    # FGW compare: entropic and exact
    fgw_pairs = [
        ("fgw_pgd", "pot_fgw_pgd"),
        ("fgw_pgd_double", "pot_fgw_pgd"),
        ("fgw_pgd_mixed", "pot_fgw_pgd"),
        ("fgw_pgd_mixed_lr64", "pot_fgw_pgd"),
        ("fgw_exact_cg", "pot_fgw_exact"),
    ]
    for rf_method, pot_method in fgw_pairs:
        for n in sorted({int(r["n"]) for r in rfugw_rows if r["suite"] == "fgw" and r["method"] == rf_method}):
            rk = ("fgw", n, rf_method)
            pk = ("fgw", n, pot_method)
            if rk in rf_idx and pk in pot_idx:
                r = rf_idx[rk]
                p = pot_idx[pk]
                rf_ms = to_float(r, "median_ms")
                pot_ms = to_float(p, "median_ms")
                speedup = pot_ms / rf_ms if rf_ms > 0 else float("nan")
                comp_rows.append([f"fgw:{rf_method}", str(n), fmt(rf_ms), fmt(pot_ms), fmt(speedup)])

    # FUGW compare: fugw_kl / fugw_kl_mixed vs pot_fugw_kl
    fugw_methods = ("fugw_kl", "fugw_kl_mixed")
    for rf_method in fugw_methods:
        for n in sorted({int(r["n"]) for r in rfugw_rows if r["suite"] == "fugw" and r["method"] == rf_method}):
            rk = ("fugw", n, rf_method)
            pk = ("fugw", n, "pot_fugw_kl")
            if rk in rf_idx and pk in pot_idx:
                r = rf_idx[rk]
                p = pot_idx[pk]
                rf_ms = to_float(r, "median_ms")
                pot_ms = to_float(p, "median_ms")
                speedup = pot_ms / rf_ms if rf_ms > 0 else float("nan")
                comp_rows.append([f"fugw:{rf_method}", str(n), fmt(rf_ms), fmt(pot_ms), fmt(speedup)])

    if not comp_rows:
        lines.append("No overlapping rows to compare.\n")
        return "\n".join(lines)

    lines.append(md_table(["suite", "n", "rfugw_median_ms", "pot_median_ms", "pot/rfugw"], comp_rows))
    lines.append("")
    return "\n".join(lines)


def section_threads(rows):
    if not rows:
        return "## Thread Scaling\n\nNo thread scaling CSV found.\n"

    lines = ["## Thread Scaling\n"]
    baseline = {}
    for r in rows:
        if int(r["thread_count"]) == 1:
            baseline[(r["suite"], r["method"])] = to_float(r, "median_ms")

    trows = []
    for r in sorted(rows, key=lambda x: (x["suite"], x["method"], int(x["thread_count"]))):
        key = (r["suite"], r["method"])
        base = baseline.get(key, float("nan"))
        cur = to_float(r, "median_ms")
        speedup = (base / cur) if (base > 0 and cur > 0) else float("nan")
        trows.append(
            [
                r["suite"],
                r["method"],
                r["thread_count"],
                fmt(cur),
                fmt(speedup),
                fmt(to_float(r, "iter_per_sec")),
            ]
        )

    lines.append(md_table(["suite", "method", "threads", "median_ms", "speedup_vs_1", "iter_per_sec"], trows))
    lines.append("")
    return "\n".join(lines)


def main():
    rfugw_csv = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path("rfugw/inst/bench/results/benchmark_latest.csv")
    pot_csv = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("rfugw/inst/bench/results/pot_benchmark_latest.csv")
    thread_csv = Path(sys.argv[3]) if len(sys.argv) >= 4 else Path("rfugw/inst/bench/results/thread_scaling_latest.csv")
    out_md = Path(sys.argv[4]) if len(sys.argv) >= 5 else Path("rfugw/inst/bench/results/benchmark_report.md")

    rf_rows = read_csv(rfugw_csv)
    pot_rows = read_csv(pot_csv)
    th_rows = read_csv(thread_csv)

    out_md.parent.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now().isoformat(timespec="seconds")
    chunks = [
        "# Benchmark Report",
        "",
        f"Generated: `{now}`",
        "",
        f"Inputs: `{rfugw_csv}`, `{pot_csv}`, `{thread_csv}`",
        "",
        section_rfugw(rf_rows),
        section_pot_compare(rf_rows, pot_rows),
        section_threads(th_rows),
    ]
    out_md.write_text("\n".join(chunks))
    print(f"Wrote benchmark report: {out_md}")


if __name__ == "__main__":
    main()
