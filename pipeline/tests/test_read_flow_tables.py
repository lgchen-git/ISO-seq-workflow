#!/usr/bin/env python3
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_read_flow_tables import build


class ReadFlowTableTest(unittest.TestCase):
    def test_emits_separate_sample_views_with_exact_conserved_fates(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2] / "tmp") as temp:
            root = Path(temp)
            sample_table = root / "samples.tsv"
            sample_table.write_text("sample_id\tfastq\nS1\tinputs/S1.fastq\n", encoding="utf-8")
            sample = root / "outputs/S1"
            summaries = {
                "01_raw_qc/qc.json": {"reads": 100},
                "02_trim/trim_summary.json": {
                    "input_reads": 100, "dual_adapter_reads": 80,
                    "adapter_unassigned_reads": 20, "retained_reads": 70,
                    "short_insert_reads": 10,
                },
                "03_alignment/alignment.json": {
                    "input_reads": 70, "primary_q20_reads": 50,
                    "low_mapq_reads": 12, "unmapped_reads": 8,
                },
            }
            for name, value in summaries.items():
                path = sample / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(value), encoding="utf-8")
            destination = root / "flow"
            build(sample_table, root / "outputs", destination, 20)
            with (destination / "read_flow_nodes.input.tsv").open(encoding="utf-8") as handle:
                nodes = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["view_id"] for row in nodes}, {"S1_reads"})
            self.assertEqual({row["unit"] for row in nodes}, {"reads"})
            self.assertEqual({row["scope"] for row in nodes}, {"sample_library"})
            values = {row["node_id"]: int(row["value"]) for row in nodes}
            self.assertEqual(values["raw"], 100)
            self.assertEqual(values["adapter_unassigned"], 20)
            self.assertEqual(values["short_insert"], 10)
            self.assertEqual(values["primary_q20"], 50)

    def test_rejects_nonconserved_source_summaries(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2] / "tmp") as temp:
            root = Path(temp)
            sample_table = root / "samples.tsv"
            sample_table.write_text("sample_id\tfastq\nS1\tinputs/S1.fastq\n", encoding="utf-8")
            sample = root / "outputs/S1"
            summaries = {
                "01_raw_qc/qc.json": {"reads": 100},
                "02_trim/trim_summary.json": {
                    "input_reads": 100, "dual_adapter_reads": 79,
                    "adapter_unassigned_reads": 20, "retained_reads": 70,
                    "short_insert_reads": 9,
                },
                "03_alignment/alignment.json": {
                    "input_reads": 70, "primary_q20_reads": 50,
                    "low_mapq_reads": 12, "unmapped_reads": 8,
                },
            }
            for name, value in summaries.items():
                path = sample / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fates do not conserve"):
                build(sample_table, root / "outputs", root / "flow", 20)


if __name__ == "__main__":
    unittest.main()
