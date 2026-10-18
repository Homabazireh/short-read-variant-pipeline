#!/usr/bin/env python3
"""Report local prerequisites without installing or changing the environment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("apptainer", "docker"), default="apptainer")
    parser.add_argument("--workdir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    checks = {"python_310_or_newer": sys.version_info >= (3, 10)}
    tools = {}
    for name in ("java", "nextflow", "samtools", "bcftools", args.profile):
        tools[name] = shutil.which(name)
        checks[name] = tools[name] is not None
    if args.profile == "docker" and checks["docker"]:
        proc = subprocess.run(["docker", "info", "--format", "{{.ServerVersion}}"],
                              capture_output=True, text=True, timeout=20)
        checks["docker_daemon_accessible"] = proc.returncode == 0
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        checks["avx_for_deepvariant"] = " avx " in cpuinfo.read_text()
    free_gib = round(shutil.disk_usage(args.workdir).free / 1024**3, 1)
    print(json.dumps({"checks": checks, "tools": tools, "free_disk_gib": free_gib,
                      "note": "Available disk is informational; full datasets, caches and temporary work need a run-specific storage budget."}, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
