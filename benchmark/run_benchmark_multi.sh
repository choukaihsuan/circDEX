#!/bin/bash
# run_benchmark_multi.sh -- orchestrate the 2 additional RNase R +/- pairs
# (PRJNA789110: NCI-H23, SW480) one dataset at a time, with disk-space
# logging before/after each dataset and intermediate-file cleanup once a
# dataset's accuracy_summary.tsv is safely produced.
#
# Usage (from circDEX/ root, inside conda env ciriquant):
#   setsid nohup bash benchmark/run_benchmark_multi.sh > logs/bench_multi/run_all.log 2>&1 < /dev/null &
#   disown
#
# Datasets are processed strictly one at a time (Snakemake --cores 8 with
# every rule requesting 8 threads already serializes jobs; this script adds
# an explicit outer loop so cleanup happens between datasets too).

set -uo pipefail
cd "$(dirname "$0")/.."   # repo root

LOGDIR=logs/bench_multi
mkdir -p "$LOGDIR"

DATASETS=(prjna789110_ncih23 prjna789110_sw480)

log_disk() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOGDIR/disk_usage.log"
    df -h /home3 >> "$LOGDIR/disk_usage.log"
    echo >> "$LOGDIR/disk_usage.log"
}

for ds in "${DATASETS[@]}"; do
    echo "=============================================================="
    echo "Dataset: $ds   (started $(date))"
    echo "=============================================================="
    log_disk "BEFORE $ds"

    # 2026-10-06: cores raised 8->16 so the two per-dataset downloads (each
    # now threads:8, aria2c -x12 -s12) can run in parallel instead of being
    # serialized -- halves download wall time. 2x12=24 S3 connections is
    # still well inside CLAUDE.md's documented <=48-connection safe ceiling.
    snakemake \
        --snakefile benchmark/Snakefile_multi \
        --configfile benchmark/config_benchmark_multi.yaml \
        --cores 16 \
        --resources mem_gb=60 \
        --keep-going \
        --rerun-incomplete \
        --rerun-triggers mtime \
        -- "results/benchmark_multi/${ds}/accuracy_summary.tsv" \
        2>&1 | tee -a "$LOGDIR/snakemake_${ds}.log"

    status=$?
    log_disk "AFTER $ds (snakemake exit=$status)"

    if [ $status -ne 0 ]; then
        echo "[run_benchmark_multi] $ds FAILED (exit=$status) -- stopping, intermediate files kept for debugging." \
            | tee -a "$LOGDIR/run_all.log"
        exit 1
    fi

    # ── Cleanup: remove raw/trimmed/BAM/star_tmp for this dataset only ──────
    # (config_benchmark_multi.yaml raw_dir/trimmed_dir), now that
    # accuracy_summary.tsv exists. Per-sample circRNA/{srr}/ working dirs
    # (BAM, mate1/mate2, find_circ unmapped fastq) are also cleared; the
    # small detection outputs (.gtf/.bed/DCC/CIRI2/CIRCexplorer2/find_circ
    # splice_sites.bed) needed for provenance are kept.
    python3 - "$ds" <<'PYEOF'
import sys, yaml, shutil, glob, os
ds = sys.argv[1]
cfg = yaml.safe_load(open("benchmark/config_benchmark_multi.yaml"))
d = cfg["datasets"][ds]
for p in (d["raw_dir"], d["trimmed_dir"]):
    if os.path.isdir(p):
        shutil.rmtree(p, ignore_errors=True)
        print(f"[cleanup] removed {p}")
res_base = f"results/benchmark_multi/{ds}"
for p in (f"{res_base}/.raw", f"{res_base}/.trim", f"{res_base}/.sub"):
    if os.path.isdir(p):
        shutil.rmtree(p, ignore_errors=True)
        print(f"[cleanup] removed {p}")
for srr in [d["total_rna"]] + d["rnaser"]:
    for sub in ("Aligned.sortedByCoord.out.bam", "Aligned.sortedByCoord.out.bam.bai",
                "star_tmp", "mate1", "mate2",
                "find_circ/unmapped_1.fastq", "find_circ/unmapped_2.fastq"):
        p = f"{res_base}/circRNA/{srr}/{sub}"
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
            print(f"[cleanup] removed dir {p}")
        elif os.path.isfile(p):
            os.remove(p)
            print(f"[cleanup] removed file {p}")
PYEOF

    log_disk "AFTER CLEANUP $ds"
    echo "Dataset $ds done (finished $(date))"
done

echo "=============================================================="
echo "All datasets complete. Run the Stage 4 aggregation script next:"
echo "  python scripts/aggregate_benchmark_multi.py"
echo "=============================================================="
