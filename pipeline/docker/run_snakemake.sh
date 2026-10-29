#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)
project=${ISOSEQ_PROJECT_ROOT:-$root}
project=$(cd "$project" && pwd -P)
config=${1:-config/run.yaml}
if [[ $# -gt 0 ]]; then shift; fi
cores=${ISOSEQ_CORES:-32}
mem=${ISOSEQ_MEMORY_MB:-131072}
image=${ISOSEQ_IMAGE:-iso-seq-bulk:1.1.0}
reference=${ISOSEQ_REFERENCE_ROOT:-$project/reference}
samples_tsv=${ISOSEQ_SAMPLES_TSV:-config/samples.tsv}
contrasts_tsv=${ISOSEQ_CONTRASTS_TSV:-config/contrasts.tsv}
run_id=$(date +%Y%m%d_%H%M%S)
outputs=${ISOSEQ_OUTPUTS_ROOT:-outputs/iso-seq_bulk_${run_id}}
scratch=${ISOSEQ_SCRATCH_ROOT:-tmp/$(basename "$outputs")}

case "$config" in
  /*|../*|*/../*) echo "Config must be inside the project: $config" >&2; exit 2 ;;
esac
case "$samples_tsv" in
  /*|../*|*/../*) echo "Sample sheet must be inside the project: $samples_tsv" >&2; exit 2 ;;
esac
case "$contrasts_tsv" in
  /*|../*|*/../*) echo "Contrast table must be inside the project: $contrasts_tsv" >&2; exit 2 ;;
esac
case "$outputs" in
  outputs/*) ;;
  *) echo "Output root must be under outputs/: $outputs" >&2; exit 2 ;;
esac
case "$scratch" in
  tmp/*) ;;
  *) echo "Scratch root must be under tmp/: $scratch" >&2; exit 2 ;;
esac
outputs_abs=$(realpath -m "$project/$outputs")
case "$outputs_abs" in
  "$project"/outputs/*) ;;
  *) echo "Output root escapes outputs/: $outputs" >&2; exit 2 ;;
esac
marker="$outputs_abs/_project/RUN_COMPLETE.tsv"
if [[ -s "$marker" ]] && awk -F '\t' 'NR==2 && $1=="completed" {ok=1} END {exit !ok}' "$marker"; then
  echo "Output root already contains a completed run; choose a fresh ISOSEQ_OUTPUTS_ROOT: $outputs" >&2
  exit 2
fi
if [[ -z "${ISOSEQ_OUTPUTS_ROOT+x}" && -e "$outputs_abs" && -n "$(find "$outputs_abs" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Output directory is not empty; choose a fresh ISOSEQ_OUTPUTS_ROOT: $outputs" >&2
  exit 2
fi
scratch_abs=$(realpath -m "$project/$scratch")
case "$scratch_abs" in
  "$project"/tmp/*) ;;
  *) echo "Scratch root escapes tmp/: $scratch" >&2; exit 2 ;;
esac
[[ "$cores" =~ ^[1-9][0-9]*$ ]] || { echo "ISOSEQ_CORES must be a positive integer" >&2; exit 2; }
[[ "$mem" =~ ^[1-9][0-9]*$ ]] || { echo "ISOSEQ_MEMORY_MB must be a positive integer" >&2; exit 2; }

for d in config inputs log outputs pipeline reference tmp; do
  [[ -d "$project/$d" ]] || { echo "Missing directory: $project/$d" >&2; exit 2; }
done
[[ -f "$project/$config" && -f "$project/$samples_tsv" && -f "$project/$contrasts_tsv" && -f "$project/pipeline/Snakefile" && -f "$project/pipeline/scripts/render_layout_config.py" ]] || {
  echo "Missing config, sample/contrast table, Snakefile or layout renderer" >&2; exit 2;
}
reference=$(cd "$reference" 2>/dev/null && pwd -P) || {
  echo "Reference directory does not exist: ${ISOSEQ_REFERENCE_ROOT:-$project/reference}" >&2; exit 2;
}
[[ -s "$reference/fasta/genome.fa" && -s "$reference/genes/genes.gtf.gz" ]] || {
  echo "Reference must contain fasta/genome.fa and genes/genes.gtf.gz: $reference" >&2; exit 2;
}

python3 "$project/pipeline/scripts/design.py" \
  --samples "$project/$samples_tsv" --contrasts "$project/$contrasts_tsv" \
  --project-root "$project" --min-replicates 2

mkdir -p "$project/$outputs" "$project/$scratch" "$project/log" "$project/tmp/snakemake" \
  "$project/tmp/cache" "$project/tmp/docker-client"
python3 "$project/pipeline/scripts/render_layout_config.py" \
  --project-root "$project" \
  --template "$project/config/analysis_layout.yaml" \
  --output "$project/tmp/analysis_layout.runtime.yaml"

# The analysis scripts re-read the run config from disk, so the effective run paths must
# live in a config file rather than only in --config overrides.
runtime_config="tmp/run.${run_id}.yaml"
awk -v o="$outputs" -v s="$scratch" -v sm="$samples_tsv" -v ct="$contrasts_tsv" '
  $1=="outputs_root:"  {print "outputs_root: " o; next}
  $1=="scratch_root:"  {print "scratch_root: " s; next}
  $1=="samples_tsv:"   {print "samples_tsv: " sm; next}
  $1=="contrasts_tsv:" {print "contrasts_tsv: " ct; next}
  {print}
' "$project/$config" > "$project/$runtime_config"
for key in outputs_root scratch_root samples_tsv contrasts_tsv; do
  grep -q "^${key}: " "$project/$runtime_config" || {
    echo "Runtime config is missing ${key}: $runtime_config" >&2; exit 2;
  }
done
mkdir -p "$project/$outputs/_project/inputs_used"
cp "$project/$runtime_config" "$project/$outputs/_project/inputs_used/"
cp "$project/$samples_tsv" "$project/$contrasts_tsv" "$project/$outputs/_project/inputs_used/"

if [[ -d "$project/.snakemake" && ! -L "$project/.snakemake" ]]; then
  mv "$project/.snakemake" "$project/tmp/snakemake-previous-$(date +%Y%m%d%H%M%S)"
elif [[ -L "$project/.snakemake" && ! -e "$project/.snakemake" ]]; then
  mv "$project/.snakemake" "$project/tmp/snakemake-broken-$(date +%Y%m%d%H%M%S)"
fi
if [[ ! -e "$project/.snakemake" && ! -L "$project/.snakemake" ]]; then
  ln -s tmp/snakemake "$project/.snakemake"
fi

git_commit=$(git -C "$project" rev-parse HEAD 2>/dev/null || true)
git_worktree=$(git -C "$project" status --short 2>/dev/null || true)
export DOCKER_CONFIG="$project/tmp/docker-client"
name="isoseq-$(basename "$project")-$run_id"
logfile="$project/log/runner-$run_id.log"
dry_run=false
for arg in "$@"; do
  if [[ "$arg" == "--dry-run" || "$arg" == "-n" ]]; then dry_run=true; fi
done

docker_args=(
  run --rm --name "$name" --network=none --user "$(id -u):$(id -g)"
  --cpus "$cores" --memory "${mem}m" --memory-swap "${mem}m"
  -e TMPDIR="$project/tmp" -e TMP="$project/tmp" -e TEMP="$project/tmp"
  -e XDG_CACHE_HOME="$project/tmp/cache" -e MPLCONFIGDIR="$project/tmp/cache/matplotlib"
  -e HOME="$project/tmp/cache/home" -e PYTHONDONTWRITEBYTECODE=1
  -e ISOSEQ_GIT_COMMIT="$git_commit" -e ISOSEQ_GIT_WORKTREE="$git_worktree"
  -e ISOSEQ_CORES="$cores" -e ISOSEQ_MEMORY_MB="$mem" -e ISOSEQ_IMAGE="$image"
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -e NUMEXPR_NUM_THREADS=1
  -v "$project:$project:ro"
  -v "$project/log:$project/log:rw"
  -v "$project/outputs:$project/outputs:rw"
  -v "$project/tmp:$project/tmp:rw"
)
if [[ "$reference" != "$project/reference" ]]; then
  docker_args+=(-v "$reference:$project/reference:ro")
fi
docker_args+=(-w "$project" "$image")

if docker "${docker_args[@]}" \
  --snakefile pipeline/Snakefile --configfile "$runtime_config" \
  --config "outputs_root=$outputs" "scratch_root=$scratch" "samples_tsv=$samples_tsv" "contrasts_tsv=$contrasts_tsv" \
  --cores "$cores" --resources "mem_mb=$mem" \
  --rerun-incomplete --printshellcmds --directory "$project" "$@" all \
  2>&1 | tee "$logfile"; then
  :
else
  status=$?
  echo "Snakemake container failed with exit code $status; see $logfile" >&2
  exit "$status"
fi

if [[ "$dry_run" == false ]]; then
  report="$project/$outputs/_project/09_report/report.html"
  marker="$project/$outputs/_project/RUN_COMPLETE.tsv"
  [[ -s "$report" ]] || { echo "Missing or empty final report: $report" >&2; exit 3; }
  awk -F '\t' 'NR==2 && $1=="completed" {ok=1} END {exit !ok}' "$marker" || {
    echo "Missing or incomplete run marker: $marker" >&2; exit 3;
  }
fi
