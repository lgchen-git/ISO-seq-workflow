#!/usr/bin/env python3
"""Convert exact per-sample workflow summaries into conserved read-flow tables."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(samples_tsv: Path, outputs: Path, destination: Path, mapq: int) -> None:
    with samples_tsv.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if not samples:
        raise ValueError("Sample sheet has no rows")
    nodes: list[dict] = []
    links: list[dict] = []
    views: list[dict] = []
    for order, row in enumerate(samples):
        sample = row["sample_id"].strip()
        if not sample:
            raise ValueError("Empty sample_id")
        view = f"{sample}_reads"
        raw = load(outputs / sample / "01_raw_qc/qc.json")["reads"]
        trim = load(outputs / sample / "02_trim/trim_summary.json")
        alignment = load(outputs / sample / "03_alignment/alignment.json")
        raw_trim = int(trim["input_reads"])
        dual = int(trim["dual_adapter_reads"])
        adapter_loss = int(trim["adapter_unassigned_reads"])
        retained = int(trim["retained_reads"])
        short = int(trim["short_insert_reads"])
        aligned_input = int(alignment["input_reads"])
        q20 = int(alignment["primary_q20_reads"])
        low = int(alignment["low_mapq_reads"])
        unmapped = int(alignment["unmapped_reads"])
        if raw != raw_trim:
            raise ValueError(f"{sample}: raw FASTQ ({raw}) != cutadapt input ({raw_trim})")
        if raw != dual + adapter_loss or dual != retained + short:
            raise ValueError(f"{sample}: trimming fates do not conserve raw reads")
        if retained != aligned_input or aligned_input != q20 + low + unmapped:
            raise ValueError(f"{sample}: alignment fates do not conserve retained reads")
        chain = [
            (0, "raw", "Raw reads", raw, "retained"),
            (1, "dual_adapter", "Reads with both adapters", dual, "retained"),
            (1, "adapter_unassigned", "No complete adapter pair", adapter_loss, "loss"),
            (2, "retained", "Retained inserts", retained, "retained"),
            (2, "short_insert", "Insert below length threshold", short, "loss"),
            (3, "primary_q20", f"Primary mapped, MAPQ ≥ {mapq}", q20, "retained"),
            (3, "low_mapq", f"Mapped below MAPQ {mapq}", low, "loss"),
            (3, "unmapped", "Unmapped reads", unmapped, "loss"),
        ]
        nodes.extend({"view_id": view, "node_id": node, "stage": stage,
                      "label": label, "value": count, "unit": "reads",
                      "scope": "sample_library", "node_kind": kind}
                     for stage, node, label, count, kind in chain)
        edges = [
            ("raw", "dual_adapter", dual, "Cutadapt detected both adapter ends"),
            ("raw", "adapter_unassigned", adapter_loss, "Cutadapt did not detect both adapter ends"),
            ("dual_adapter", "retained", retained, "Insert meets configured minimum length"),
            ("dual_adapter", "short_insert", short, "Insert is shorter than configured minimum"),
            ("retained", "primary_q20", q20, f"Primary alignment meets MAPQ {mapq}"),
            ("retained", "low_mapq", low, f"Primary alignment below MAPQ {mapq}"),
            ("retained", "unmapped", unmapped, "No mapped primary alignment"),
        ]
        links.extend({"view_id": view, "link_id": f"{source}_to_{target}",
                      "source": source, "target": target, "value": count, "reason": reason}
                     for source, target, count, reason in edges)
        views.append({"view_id": view, "display_name": sample, "order": order,
                      "is_default": "true" if order == 0 else "false"})
    write(destination / "read_flow_nodes.input.tsv",
          ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"], nodes)
    write(destination / "read_flow_links.input.tsv",
          ["view_id", "link_id", "source", "target", "value", "reason"], links)
    write(destination / "read_flow_views.tsv", ["view_id", "display_name", "order", "is_default"], views)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--outputs", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--mapq", type=int, default=20)
    args = parser.parse_args()
    build(args.samples, args.outputs, args.destination, args.mapq)


if __name__ == "__main__":
    main()
