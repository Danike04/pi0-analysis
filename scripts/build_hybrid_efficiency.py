#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_table(path):
    comments = []
    rows = []

    with Path(path).open() as f:
        for raw in f:
            line = raw.rstrip("\n")

            if not line.strip():
                continue

            if line.lstrip().startswith("#"):
                comments.append(line)
                continue

            cols = line.split()
            if len(cols) < 12:
                raise RuntimeError(
                    f"{path}: expected at least 12 columns, got {len(cols)}"
                )

            values = [float(x) for x in cols[:12]]

            key = (
                values[0],  # E_low
                values[1],  # E_high
                values[3],  # theta_low
                values[4],  # theta_high
            )

            rows.append((key, line))

    return comments, rows


def same_bin(value, target):
    return abs(value - target) < 1.0e-6


def main():
    parser = argparse.ArgumentParser(
        description="Replace selected broad-MC paper bins with dedicated high-stat MC tables."
    )

    parser.add_argument("--broad", required=True)
    parser.add_argument(
        "--replace",
        nargs=3,
        action="append",
        metavar=("ELOW", "EHIGH", "TABLE"),
        default=[],
    )
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    comments, broad_rows = read_table(args.broad)

    replacement_rows = {}
    provenance = []

    for elo_s, ehi_s, table in args.replace:
        elo = float(elo_s)
        ehi = float(ehi_s)

        _, dedicated_rows = read_table(table)

        broad_keys = {
            key
            for key, _ in broad_rows
            if same_bin(key[0], elo) and same_bin(key[1], ehi)
        }

        dedicated = {
            key: line
            for key, line in dedicated_rows
            if same_bin(key[0], elo) and same_bin(key[1], ehi)
        }

        if not broad_keys:
            raise RuntimeError(f"Bin {elo:g}-{ehi:g} not found in broad table")

        if set(dedicated) != broad_keys:
            raise RuntimeError(
                f"Angular-bin mismatch for {elo:g}-{ehi:g}: "
                f"broad has {len(broad_keys)}, dedicated has {len(dedicated)}"
            )

        replacement_rows.update(dedicated)
        provenance.append(
            f"# hybrid replacement: {elo:g}-{ehi:g} MeV <- {table}"
        )

    output = Path(args.output)

    with output.open("w") as f:
        f.write("# Hybrid paper-bin detection-efficiency table\n")
        f.write(f"# broad source: {args.broad}\n")

        for line in provenance:
            f.write(line + "\n")

        for line in comments:
            f.write(line + "\n")

        for key, line in broad_rows:
            f.write(replacement_rows.get(key, line) + "\n")

    print(f"Wrote {output}")
    print(f"Replaced rows: {len(replacement_rows)}")


if __name__ == "__main__":
    main()
