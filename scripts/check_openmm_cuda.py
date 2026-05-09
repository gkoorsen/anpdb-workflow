"""Preflight check for OpenMM CUDA on a GPU workstation.

Usage:
  python scripts/check_openmm_cuda.py
  python scripts/check_openmm_cuda.py --device-index 0 --precision mixed
"""

from __future__ import annotations

import argparse
import platform as py_platform
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-index", default="0")
    parser.add_argument("--precision", default="mixed", choices=["single", "mixed", "double"])
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()

    try:
        import openmm
        from openmm import Context, LangevinMiddleIntegrator, Platform, System, Vec3, unit
    except Exception as exc:  # pragma: no cover - environment-specific
        print(f"OpenMM import failed: {exc}", file=sys.stderr)
        return 2

    print(f"Python: {sys.version.split()[0]} ({py_platform.platform()})")
    print(f"OpenMM: {openmm.version.version}")

    platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print("Available OpenMM platforms: " + ", ".join(platforms))

    try:
        cuda = Platform.getPlatformByName("CUDA")
    except Exception as exc:
        print(f"CUDA platform unavailable: {exc}", file=sys.stderr)
        return 0 if args.allow_cpu else 1

    properties = {
        "DeviceIndex": str(args.device_index),
        "Precision": args.precision,
    }
    system = System()
    system.addParticle(39.948 * unit.amu)
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)

    try:
        context = Context(system, integrator, cuda, properties)
        context.setPositions([Vec3(0, 0, 0)] * unit.nanometer)
        context.getState(getEnergy=True)
    except Exception as exc:
        print(f"CUDA context creation failed: {exc}", file=sys.stderr)
        return 1

    print("CUDA context: OK")
    for prop in ("DeviceIndex", "DeviceName", "Precision"):
        try:
            print(f"{prop}: {cuda.getPropertyValue(context, prop)}")
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
