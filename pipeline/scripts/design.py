"""Sample/contrast table parsing and preflight validation for ISO-seq."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path


SAMPLE_COLUMNS = {"sample_id", "fastq", "group", "batch"}
CONTRAST_COLUMNS = {"contrast_id", "test_group", "reference_group"}
ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")


def read_table(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Table has no header: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Table has duplicate column names: {path}")
        return [
            {key: (value or "").strip() for key, value in row.items() if key is not None}
            for row in reader
        ]


def _rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    if not matrix:
        return 0
    a = [row[:] for row in matrix]
    nrow, ncol = len(a), len(a[0])
    rank = 0
    for col in range(ncol):
        if rank >= nrow:
            break
        pivot = max(range(rank, nrow), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) <= tolerance:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [v / scale for v in a[rank]]
        for row in range(nrow):
            if row != rank:
                scale = a[row][col]
                a[row] = [v - scale * p for v, p in zip(a[row], a[rank])]
        rank += 1
        if rank == nrow:
            break
    return rank


def _design_rank(rows: list[dict[str, str]]) -> tuple[int, int]:
    groups = sorted({row["group"] for row in rows})
    batches = sorted({row["batch"] for row in rows if row["batch"]})
    columns = [[1.0] for _ in rows]
    for level in groups[1:]:
        for i, row in enumerate(rows):
            columns[i].append(float(row["group"] == level))
    for level in batches[1:]:
        for i, row in enumerate(rows):
            columns[i].append(float(row["batch"] == level))
    return _rank(columns), len(columns)


def load_design(samples_path: str | Path, contrasts_path: str | Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    samples = read_table(samples_path)
    contrasts = read_table(contrasts_path)
    for path, required in (
        (samples_path, SAMPLE_COLUMNS),
        (contrasts_path, CONTRAST_COLUMNS),
    ):
        with Path(path).open(encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        actual = set(header)
        missing = sorted(required - actual)
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    return samples, contrasts


def validate_design(
    samples_path: str | Path,
    contrasts_path: str | Path,
    project_root: str | Path,
    min_replicates: int = 2,
    check_fastq: bool = True,
) -> dict:
    samples, contrasts = load_design(samples_path, contrasts_path)
    root = Path(project_root).resolve()
    inputs_root = (root / "inputs").resolve()
    if min_replicates < 2:
        raise ValueError("min_replicates must be at least 2 for inferential analysis")
    if not samples:
        raise ValueError("Sample table is empty")

    seen_samples: set[str] = set()
    for row in samples:
        sid = row.get("sample_id", "")
        fastq = row.get("fastq", "")
        if not ID_RE.fullmatch(sid) or sid in seen_samples:
            raise ValueError(f"Invalid or duplicate sample_id: {sid!r}")
        seen_samples.add(sid)
        if check_fastq:
            path = (root / fastq).resolve()
            try:
                path.relative_to(inputs_root)
            except ValueError as exc:
                raise ValueError(f"FASTQ must resolve under inputs/: {fastq}") from exc
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing or empty FASTQ: {fastq}")

    batches = [row.get("batch", "") for row in samples]
    if any(batches) and not all(batches):
        raise ValueError("batch must be blank for every sample or populated for every sample")
    assigned_groups = [row.get("group", "") for row in samples]
    if contrasts and not all(assigned_groups):
        raise ValueError("Every sample needs a group when contrasts are requested")
    if not contrasts:
        return {
            "stats_enabled": False,
            "sample_count": len(samples),
            "groups": sorted({g for g in assigned_groups if g}),
            "group_counts": {g: assigned_groups.count(g) for g in sorted(set(assigned_groups) - {""})},
            "batch_enabled": len(set(b for b in batches if b)) > 1,
            "contrasts": [],
        }
    group_counts = {g: assigned_groups.count(g) for g in sorted(set(assigned_groups))}
    group_names = set(group_counts)
    seen_contrasts: set[str] = set()
    for contrast in contrasts:
        cid = contrast.get("contrast_id", "")
        test = contrast.get("test_group", "")
        reference = contrast.get("reference_group", "")
        if not ID_RE.fullmatch(cid) or cid in seen_contrasts:
            raise ValueError(f"Invalid or duplicate contrast_id: {cid!r}")
        seen_contrasts.add(cid)
        if not test or not reference or test == reference:
            raise ValueError(f"Contrast {cid}: test_group and reference_group must be distinct and non-empty")
        unknown = sorted({test, reference} - group_names)
        if unknown:
            raise ValueError(f"Contrast {cid} refers to unknown groups: {', '.join(unknown)}")
        for group in (test, reference):
            if group_counts[group] < min_replicates:
                raise ValueError(
                    f"Contrast {cid}: group {group} has {group_counts[group]} biological replicates; "
                    f"at least {min_replicates} are required"
                )
        subset = [row for row in samples if row["group"] in {test, reference}]
        rank, n = _design_rank(subset)
        expected_columns = 2 + max(0, len({row["batch"] for row in subset if row["batch"]}) - 1)
        if rank < expected_columns:
            raise ValueError(f"Contrast {cid}: group and batch effects are confounded (design matrix is not full rank)")
        if n - rank <= 0:
            raise ValueError(f"Contrast {cid}: design has no residual degrees of freedom")

    return {
        "stats_enabled": True,
        "sample_count": len(samples),
        "groups": sorted(group_counts),
        "group_counts": group_counts,
        "batch_enabled": len(set(b for b in batches if b)) > 1,
        "batch_uninformative": bool(any(batches) and len(set(batches)) == 1),
        "contrasts": contrasts,
        "min_replicates_per_group": min_replicates,
        "low_power_groups": sorted(g for g, n in group_counts.items() if n < 3),
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--contrasts", required=True)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--skip-fastq-check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = validate_design(
            args.samples,
            args.contrasts,
            args.project_root,
            min_replicates=args.min_replicates,
            check_fastq=not args.skip_fastq_check,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    serialized = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
