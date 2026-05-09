"""Run a list of MD production configs sequentially.

Usage:
  python scripts/md_batch.py configs/md/production/*_rep1.toml --resume
  python scripts/md_batch.py configs/md/production/*.toml --dry-run
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--allow-non-production", action="store_true")
    args = parser.parse_args()

    runner = Path(__file__).with_name("md_production.py")
    for config in args.configs:
        cmd = [sys.executable, str(runner), "--config", str(config)]
        if args.resume:
            cmd.append("--resume")
        if args.dry_run:
            cmd.append("--dry-run")
        if args.allow_non_production:
            cmd.append("--allow-non-production")
        print("\n$ " + " ".join(cmd), flush=True)
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
