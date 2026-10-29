#!/usr/bin/env python3
"""Resolve the portable analysis-layout template for a concrete checkout."""
import argparse
import json
import re
from pathlib import Path


ROOT_FIELDS = {
    "project_root": ".",
    "inputs_root": "inputs",
    "reference_root": "reference",
    "log_root": "log",
    "outputs_root": "outputs",
    "pipeline_root": "pipeline",
    "tmp_root": "tmp",
    "snakefile": "pipeline/Snakefile",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--template", type=Path, default=Path("config/analysis_layout.yaml"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    template = args.template.resolve()
    output = args.output.resolve()
    try:
        output.relative_to((root / "tmp").resolve())
    except ValueError:
        raise SystemExit("Runtime layout config must be written under project tmp/")
    if not template.is_file():
        raise SystemExit(f"Missing layout template: {template}")

    lines = template.read_text(encoding="utf-8").splitlines()
    seen = set()
    rendered = []
    for line in lines:
        match = re.match(r"^(project_root|inputs_root|reference_root|log_root|outputs_root|pipeline_root|tmp_root|snakefile):\s*.*$", line)
        if match:
            field = match.group(1)
            if field in seen:
                raise SystemExit(f"Duplicate layout field: {field}")
            seen.add(field)
            value = (root / ROOT_FIELDS[field]).resolve()
            rendered.append(f"{field}: {json.dumps(str(value))}")
        else:
            rendered.append(line)
    missing = sorted(set(ROOT_FIELDS) - seen)
    if missing:
        raise SystemExit(f"Layout template is missing fields: {', '.join(missing)}")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(rendered) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
