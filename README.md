# Bulk ISO-seq reproducible workflow

Snakemake runs the analysis stages in order and in parallel where appropriate. The workflow accepts any number of independent biological samples and groups. Runtime analysis is offline: Docker runs with `--network=none`; there are no LLM, API, or external-service calls.

This repository contains workflow code, tests, documentation, and configuration templates. Raw sequencing data, genome references, run outputs, logs, and local run tables are intentionally kept outside version control.

## One-time project setup

Keep the project tree in this layout:

```text
config/  inputs/  log/  outputs/  pipeline/  reference/  tmp/
```

After cloning, create local run configuration from the templates:

```bash
cp config/run.example.yaml config/run.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/contrasts.example.tsv config/contrasts.tsv
```

Place FASTQ files below `inputs/`. Place the matching genome and annotation under `reference/fasta/genome.fa` and `reference/genes/genes.gtf.gz`. The reference is mounted read-only and is not included in the image. Use the same reference release across samples in one run.

Build and export the runtime image on a connected build host:

```bash
bash pipeline/docker/build_export.sh --export
```

This creates image `iso-seq-bulk:1.1.0` and archive `outputs/_docker/iso-seq-bulk_1.1.0.tar.gz`. On an offline analysis server, transfer the project, reference, inputs, and archive, then load the image:

```bash
docker load -i outputs/_docker/iso-seq-bulk_1.1.0.tar.gz
```

## Per-run tables and command

Edit the local copies of these tables for each analysis:

- `config/samples.tsv`: one row per independent biological replicate; columns are `sample_id`, `fastq`, `group`, and optional `batch`. FASTQ paths are project-relative and must resolve under `inputs/`. Leave all groups blank when requesting descriptive analysis only. When using batch, fill it for every sample.
- `config/contrasts.tsv`: one row per requested comparison with `contrast_id`, `test_group`, and `reference_group`. Comparisons are explicit; the workflow does not generate unrequested pairwise tests. Keep only the header for descriptive analysis.

Example:

```tsv
sample_id	fastq	group	batch
S01	inputs/S01.fastq.gz	Control	Batch1
S02	inputs/S02.fastq.gz	Control	Batch2
S03	inputs/S03.fastq.gz	Treatment	Batch1
S04	inputs/S04.fastq.gz	Treatment	Batch2
```

```tsv
contrast_id	test_group	reference_group
Treatment_vs_Control	Treatment	Control
```

Then run:

```bash
bash pipeline/docker/run_snakemake.sh
```

The runner validates the tables and FASTQ paths before analysis. It requires at least two independent replicates in each contrasted group, warns in the report when a group has only two, and stops if the group/batch design is confounded or has no residual degrees of freedom. Paired or repeated-measures designs are not supported. Each command writes to a new timestamped directory under `outputs/`; it does not overwrite earlier runs.

Gene and transcript differential expression use edgeR robust quasi-likelihood negative-binomial GLMs. DTU uses DRIMSeq with stageR stage-wise adjustment. When batch is supplied and estimable, it is included as a fixed effect. The workflow uses a 0.05 FDR threshold. Differential splice-event testing is not included; SUPPA2 PSI remains descriptive.

## Outputs and verification

The run writes per-sample QC and alignment results, shared transcript catalog and quantification, event summaries, and a Chinese offline report under the timestamped output root. The final report and `RUN_COMPLETE.tsv` are under `_project/09_report/` and `_project/`, respectively. Differential tables are under `_project/08_differential/` (`gene_de.tsv`, `transcript_de.tsv`, `dtu_genes.tsv`, and `dtu_transcripts.tsv`). Logs are in `log/`; scratch data and Snakemake state are under `tmp/`.

The report includes an offline interactive read-flow Sankey with one view per sample. Each view conserves raw reads through adapter detection, minimum insert length, and primary alignment fates; the sample views have separate denominators and must not be summed. Canonical flow tables, static SVGs, a standalone HTML, and validation output are stored at `_project/09_report/read_flow/`. The project-local renderer and Plotly asset are included under `pipeline/`, so the runtime does not use a personal skill installation or network resource. Workflow execution uses a network-disabled Docker container and makes no LLM, AI, API, or external-service calls.

The self-tests can be run in the project environment with:

```bash
python3 -m unittest discover -s pipeline/tests -p 'test_*.py'
python3 -m py_compile pipeline/scripts/*.py
bash -n pipeline/docker/run_snakemake.sh pipeline/docker/build_export.sh
```

To run the isolated end-to-end fixture, deterministically copy a small subset from selected local FASTQs, then point the runner at a read-only reference bundle:

```bash
python3 pipeline/scripts/make_fixture.py --reads 20000 --samples S02 S06 --destination tmp/fixture_manual
ISOSEQ_PROJECT_ROOT="$PWD/tmp/fixture_manual" \
ISOSEQ_REFERENCE_ROOT=/path/to/reference \
ISOSEQ_OUTPUTS_ROOT=outputs/fixture_e2e \
ISOSEQ_SCRATCH_ROOT=tmp/fixture_e2e \
bash pipeline/docker/run_snakemake.sh
```

The fixture clears group labels and contrasts, so it exercises descriptive processing without implying a statistical comparison. Completion is accepted only when the workflow writes `outputs/fixture_e2e/_project/RUN_COMPLETE.tsv` and the final report validator passes.

For a dry run, use:

```bash
bash pipeline/docker/run_snakemake.sh config/run.yaml --dry-run
```

Default resource limits are 32 CPU and 128 GiB. Override with `ISOSEQ_CORES` and `ISOSEQ_MEMORY_MB`. To choose a specific output root, set `ISOSEQ_OUTPUTS_ROOT=outputs/<run-id>`; the matching scratch root defaults to `tmp/<run-id>`.

## Interpretation

- Counts are unique-consistent reads, not UMI-deduplicated molecules.
- Gene-level differential expression uses the reference gene count matrix; transcript-level differential expression and DTU use the shared, structure-filtered transcript catalog.
- Statistical comparisons treat samples as independent biological replicates. Two replicates per group meet the workflow minimum but provide limited power.
- Sample-level PSI and transcript-end candidates are descriptive. The workflow does not infer differential splice events, validated TSS/APA sites, or RNA synthesis/degradation rates.

Project-specific run records, sample sheets, and generated results are intentionally kept outside this repository.
