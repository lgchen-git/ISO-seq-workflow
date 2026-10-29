#!/usr/bin/env python3
"""Create a deterministic, isolated small-read fixture under project tmp/."""
import argparse
import csv
import gzip
import json
import re
import shutil
from pathlib import Path


def set_yaml_scalar(text, key, value):
    pattern = re.compile(rf"^{re.escape(key)}:\s*.*$", re.MULTILINE)
    encoded = json.dumps(value, ensure_ascii=False) if isinstance(value, str) else str(value)
    updated, count = pattern.subn(f"{key}: {encoded}", text)
    if count != 1:
        raise ValueError(f"Expected one {key!r} field in run config, found {count}")
    return updated


def is_within(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=int, default=20000)
    parser.add_argument("--destination", default="tmp/fixture")
    parser.add_argument("--samples", nargs="+", help="Optional sample IDs to include; defaults to all samples")
    args = parser.parse_args()
    if args.reads < 1:
        parser.error("--reads must be positive")

    root = Path.cwd().resolve()
    destination = (root / args.destination).resolve()
    tmp_root = (root / "tmp").resolve()
    if not is_within(destination, tmp_root):
        raise SystemExit("Fixture destination must be under project tmp/")
    if destination == tmp_root:
        raise SystemExit("Fixture destination must be a child of project tmp/")
    if destination.exists() and any(destination.iterdir()):
        raise SystemExit(f"Fixture destination is not empty; choose a new path: {destination}")

    config_text = (root / "config/run.yaml").read_text(encoding="utf-8")
    sample_path = root / "config/samples.tsv"
    rows = []
    destination.mkdir(parents=True, exist_ok=True)
    for name in ["config", "inputs", "log", "outputs", "reference", "tmp"]:
        (destination / name).mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        root / "pipeline",
        destination / "pipeline",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    for name in [".gitignore", ".dockerignore"]:
        shutil.copy2(root / name, destination / name)

    with sample_path.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not source_rows:
        raise SystemExit(f"No samples in {sample_path}")
    if args.samples:
        available = {row["sample_id"].strip() for row in source_rows}
        unknown = sorted(set(args.samples) - available)
        if unknown:
            raise SystemExit(f"Unknown fixture sample IDs: {', '.join(unknown)}")
        requested = set(args.samples)
        source_rows = [row for row in source_rows if row["sample_id"].strip() in requested]

    for row in source_rows:
        sample_id = row["sample_id"].strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample_id):
            raise SystemExit(f"Invalid sample_id: {sample_id!r}")
        source = (root / row["fastq"].strip()).resolve()
        if not is_within(source, root / "inputs") or not source.is_file():
            raise SystemExit(f"FASTQ must exist under inputs/: {row['fastq']}")
        target = destination / "inputs" / f"{sample_id}.fastq.gz"
        count = 0
        with source.open("rt", encoding="ascii", newline="") as fastq, gzip.open(
            target, "wt", encoding="ascii", newline="", compresslevel=1
        ) as output:
            while count < args.reads:
                record = [fastq.readline() for _ in range(4)]
                if not record[0]:
                    break
                if any(line == "" for line in record):
                    raise SystemExit(f"Truncated FASTQ record in {source}")
                if not record[0].startswith("@") or not record[2].startswith("+"):
                    raise SystemExit(f"Malformed FASTQ record {count + 1} in {source}")
                if len(record[1].rstrip("\r\n")) != len(record[3].rstrip("\r\n")):
                    raise SystemExit(f"Sequence/quality length mismatch in {source}, record {count + 1}")
                output.writelines(record)
                count += 1
        if count == 0:
            raise SystemExit(f"No reads copied from {source}")
        rows.append({"sample_id": sample_id, "fastq": f"inputs/{target.name}"})

    with (destination / "config/samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "fastq", "group", "batch"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows({**row, "group": "", "batch": ""} for row in rows)

    (destination / "config/contrasts.tsv").write_text(
        "contrast_id\ttest_group\treference_group\n", encoding="utf-8"
    )

    changes = {
        "samples_tsv": "config/samples.tsv",
        "outputs_root": "outputs/fixture",
        "scratch_root": "tmp/analysis",
        "threads": 4,
        "catalog_threads": 8,
        "project_id": "ISO-seq-representative-fixture",
    }
    for key, value in changes.items():
        config_text = set_yaml_scalar(config_text, key, value)
    (destination / "config/run.yaml").write_text(config_text, encoding="utf-8")
    shutil.copy2(root / "config/analysis_layout.yaml", destination / "config/analysis_layout.yaml")
    print(destination)


if __name__ == "__main__":
    main()
