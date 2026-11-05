# CircDEX — circRNA Differential Expression Pipeline

An end-to-end, Snakemake-driven pipeline for circular RNA (circRNA) discovery and
differential expression analysis, from raw GEO/SRA sequencing data to an
interactive, self-contained HTML report — all driven through a browser-based
Web UI, no command-line experience required.

Given a GEO accession (or your own FASTQ files), the pipeline downloads/aligns
the data, detects circRNAs with a dual-tool consensus strategy (CIRIquant +
DCC), runs differential expression analysis (tumor vs. normal, or any two
conditions), ranks candidate biomarkers, and produces a single interactive
HTML report you can open in any browser — complete with volcano plots,
clustering heatmaps, per-circRNA structure diagrams, and miRNA/RBP binding
site predictions.

---

## Table of Contents

- [Key Features](#key-features)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
  - [1. Starting the Web UI](#1-starting-the-web-ui)
  - [2. Choosing an input method](#2-choosing-an-input-method)
  - [3. Configuring the analysis](#3-configuring-the-analysis)
  - [4. Monitoring progress](#4-monitoring-progress)
  - [5. Reading the report](#5-reading-the-report)
  - [6. Cross-dataset comparison](#6-cross-dataset-comparison)
- [Command-Line Usage](#command-line-usage)
- [Configuration Reference](#configuration-reference)
- [Pipeline Architecture](#pipeline-architecture)
- [Choosing a Good Dataset](#choosing-a-good-dataset)
- [Container Deployment (Docker / Singularity)](#container-deployment-docker--singularity)
- [Troubleshooting](#troubleshooting)
- [Project Structure](#project-structure)
- [Methodology Notes](#methodology-notes)

---

## Key Features

- **One-click GEO analysis** — paste a GEO/SRA accession (`GSEnnnnnn`, `PRJNAnnnnnn`,
  or `SRPnnnnnn`) and the pipeline downloads the data and infers tumor/normal
  labels automatically.
- **Dual-tool circRNA consensus** — CIRIquant + DCC, with confidence scoring
  and an adaptive fallback for datasets where one tool underperforms (e.g.
  short-read libraries).
- **Three differential-expression methods** run side-by-side: `edgeR_ciriquant`
  (a BSJ/FSJ-ratio test unique to this pipeline, distinguishing circRNA-specific
  regulation from linear mRNA co-regulation), `DESeq2`, and `limma-voom`.
- **Composite biomarker ranking** combining significance, fold-change, detection
  confidence, circBase novelty, and predicted miRNA/RBP binding partners.
- **Interactive self-contained HTML report** — Plotly volcano/PCA/heatmap
  charts, a hierarchical clustering heatmap with a real dendrogram, per-circRNA
  circular structure diagrams, sortable/exportable tables, and circBase
  hyperlinks — all in one file you can email or archive.
- **Cross-dataset comparison page** — once you've analyzed several datasets,
  see which circRNAs recur across cancer types (or within the same cancer
  type across independent cohorts).
- **Job queue** — submit several datasets back-to-back; they run one at a time
  with email notifications on completion.
- **Three ways to bring in data**: a GEO accession, a manual SRR list/CSV
  upload, or FASTQ files you already have on the server.

---

## Requirements

The pipeline is built around a set of bioinformatics tools that are easiest to
obtain via the provided **Docker/Singularity image** (see
[Container Deployment](#container-deployment-docker--singularity)) or a conda
environment (`envs/circrna.yaml`):

| Tool | Purpose |
|------|---------|
| Snakemake ≥ 7 | Workflow orchestration |
| CIRIquant | Primary circRNA detection (HISAT2 + BWA) |
| DCC | Secondary circRNA detection (STAR chimeric alignment) |
| STAR, HISAT2, BWA | Read alignment |
| samtools | BAM handling |
| fastp, FastQC, MultiQC | QC and trimming |
| sra-tools, aria2c | SRA/GEO data download |
| R ≥ 4.2 (edgeR, DESeq2, limma, qvalue) | Differential expression |
| Python ≥ 3.7 (pandas, plotly, scipy, Flask) | Analysis scripts + Web UI |

You will also need a **reference genome** (hg19 by default) with pre-built
BWA, HISAT2, and STAR indices, and a GTF annotation file.

The pipeline is designed to run on a Linux server or HPC cluster — a modest
laptop is not a realistic target for anything beyond inspecting the code.

### Minimum Hardware Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| OS | Linux (x86_64) | CentOS 7+ / Ubuntu 20.04+ |
| CPU | 8 cores | 16–36+ cores (parallel samples finish faster) |
| RAM | 32 GB | 64–128 GB+ |
| Free disk space | 200 GB | 500 GB – 1 TB+ |
| Network | Stable broadband (SRA downloads can be tens of GB per sample) | — |

These numbers assume human genome (hg19/hg38) alignment with **one sample at
a time**; running several samples concurrently multiplies the RAM and disk
requirements accordingly. A few things drive the requirement up in practice:

- **RAM** — STAR's human genome index alone needs ~30 GB resident in memory;
  CIRIquant (HISAT2 + BWA + CIRI2, run per sample) has been observed to peak
  around 50 GB for a single large total-RNA sample.
  32 GB is survivable for small/short-read datasets but will be tight or fail
  outright on deep total-RNA samples.
- **Disk** — beyond the ~20–30 GB reference genome + indices, intermediate
  files for a *single* sample (STAR BAM, CIRIquant's internal alignment
  files, chimeric junction files) can transiently reach 100–200 GB before
  cleanup, especially on network-attached storage where large SAM/BAM writes
  amplify usage. Raw + trimmed FASTQ for a typical dataset (6–12 samples,
  100–150 bp paired-end) commonly totals 100–300 GB by itself. Plan for
  enough headroom to hold several samples' worth of intermediates
  simultaneously if running with `--cores` high enough to parallelize.
- **CPU** — most steps (STAR, HISAT2, BWA, DCC) scale reasonably well with
  thread count; 8 cores will work but a full dataset (10+ samples) will take
  considerably longer than on a 24–36 core machine.

If your hardware falls short of the recommended column, running the pipeline
one sample/one dataset at a time (rather than relying on Snakemake to
parallelize across samples) is the main lever to keep RAM/disk peaks
manageable — set `--cores` low and let the job queue in the Web UI serialize
multiple datasets instead of running them side by side.

---

## Quick Start

The fastest way to see the pipeline work end-to-end is the Web UI's one-click
GEO launcher:

```bash
git clone <this-repo> circDEX
cd circDEX

# 1. Activate an environment with the dependencies above (conda or container)
conda activate ciriquant   # or: your equivalent environment

# 2. Edit config.yaml and config/ciriquant.yaml — see the warning below
#    (see Reference genome / Configuration Reference for exactly what to change)

# 3. Start the Web UI
python scripts/web_ui.py --host 0.0.0.0 --port 5000
```

> ⚠️ **`config.yaml` and `config/ciriquant.yaml` are checked into this repo
> as the maintainer's own working copies**, not blank templates — every path
> in them (reference genome files, index directories, tool binaries, cache
> directories) points to a specific machine and does not exist on yours. The
> pipeline will fail with confusing `FileNotFoundError`/`ConfigError`
> messages until you replace them. See
> [Reference genome](#reference-genome) below for the exact fields to edit
> before running anything.

Then open `http://<your-server-ip>:5000` in a browser, paste a GEO accession
(e.g. `GSE113230`) into the "GEO Dataset" card, pick a core count, and click
start. The status page will show live progress through download → QC →
detection → differential expression → report generation, and you'll be able
to open the finished HTML report directly from the browser once it's done.

---

## Installation

### Option A — Docker / Singularity (recommended, most portable)

A pre-built image with every dependency already installed is available on
Docker Hub. See [Container Deployment](#container-deployment-docker--singularity)
for the full walkthrough — this is the easiest way to get a working
environment on a new HPC cluster without managing conda yourself.

### Option B — Conda environment

```bash
mamba env create -f envs/circrna.yaml   # or: conda env create -f envs/circrna.yaml
conda activate ciriquant                # note: the env is named "ciriquant", not "circrna"
```

`envs/circrna.yaml` pulls everything from `bioconda` + `conda-forge` (no
`defaults` channel, to avoid Anaconda's commercial-use terms). After creating
the environment, point `config/ciriquant.yaml` at the absolute paths of the
tools inside your environment (CIRIquant validates tool paths literally, so
relative names like `bwa` will not resolve) — see the comments in that file
for the exact keys expected.

### Reference genome

**1. Get a FASTA + GTF pair.** Any matched hg19/GRCh37 pair works — common
sources are the UCSC Genome Browser's "goldenPath" downloads, GENCODE
(release 19 is the last GRCh37-based annotation), or a pre-bundled Illumina
iGenomes `hg19` package (which conveniently also ships pre-built indices).
Search for whichever of these you prefer to use and download both files.

**2. Build the aligner indices** (skip this step if you used an iGenomes
bundle that already includes them):

```bash
# BWA
bwa index -p /path/to/index/hg19 /path/to/hg19.fa

# HISAT2
hisat2-build /path/to/hg19.fa /path/to/index/hg19_hisat2_index

# STAR (--sjdbOverhang should be read_length - 1; 149 assumes 150bp reads —
# adjust to match your dataset)
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles /path/to/hg19.fa \
     --sjdbGTFfile /path/to/genes.gtf \
     --sjdbOverhang 149 \
     --runThreadN 8
```

**3. Point both config files at your paths.** `config.yaml` and
`config/ciriquant.yaml` ship in this repo pre-filled with the maintainer's
own paths (see the warning above) — every field below needs to be replaced
with your own:

| File | Field(s) | What it should point to |
|------|----------|--------------------------|
| `config.yaml` | `genome.fasta` | the FASTA file from step 1 |
| `config.yaml` | `genome.gtf` | the GTF file from step 1 |
| `config.yaml` | `genome.bwa_index` | the `-p` prefix you gave `bwa index` in step 2 |
| `config.yaml` | `genome.hisat2_index` | the prefix you gave `hisat2-build` in step 2 |
| `config.yaml` | `genome.star_index` | the `--genomeDir` you gave `STAR` in step 2 |
| `config.yaml` | `download.sra_cache_dir`, `download.tmp_dir` | any writable scratch directories on your machine |
| `config/ciriquant.yaml` | `tools.bwa`, `tools.hisat2`, `tools.samtools`, `tools.stringtie`, `tools.java` | absolute paths to each binary **inside your environment** (e.g. `which bwa` after activating it) |
| `config/ciriquant.yaml` | `reference.fasta`, `reference.gtf`, `reference.bwa_index`, `reference.hisat_index` | same values as the `config.yaml` genome paths above |

```yaml
# config.yaml — genome: block, for reference
genome:
  fasta:        /path/to/hg19.fa
  gtf:          /path/to/genes.gtf
  bwa_index:    /path/to/index/hg19          # prefix, not the .fa file itself
  hisat2_index: /path/to/index/hg19_hisat2_index
  star_index:   /path/to/star_index
  species:      hg19
```

CIRIquant checks every path in `config/ciriquant.yaml` literally (it does not
search `$PATH`), so a bare command name like `bwa` will fail — always use the
full `/path/to/env/bin/bwa`-style path there.

### Do I need to edit `config.yaml` myself?

It depends on which role you're in:

- **Setting up the environment (once, ever)** — yes. Whoever installs the
  pipeline needs to manually fill in two files before the first analysis can
  run at all:
  - `config.yaml`'s `genome:` block (reference FASTA/GTF/index paths, shown
    above)
  - `config/ciriquant.yaml` (absolute paths to `bwa`/`hisat2`/`samtools`/
    `java`/`perl` inside your environment — CIRIquant checks these paths
    literally, so bare command names like `bwa` will fail)

  These are environment-level settings, independent of which dataset you
  analyze, and normally only need to be done once per server/cluster.

- **Running an analysis (every time)** — no. Once the environment is set up,
  everything else — `project_id`, `raw_dir`/`trimmed_dir`/`results_dir`,
  `metadata`/`groups` paths, circRNA detection tool selection, DE method,
  and threshold parameters (Steps 1–3 in the Web UI) — is generated and
  saved automatically when you submit a dataset through the Web UI. Each
  submission writes its own snapshot to `config/projects/{GSE_ID}.yaml`, so
  you never need to hand-edit `config.yaml` to launch or switch between
  analyses.

In short: someone configures the reference genome and tool paths once when
setting the pipeline up; after that, day-to-day users only interact with the
Web UI's forms.

---

## Usage Guide

### 1. Starting the Web UI

```bash
cd circDEX
conda activate ciriquant   # or your environment name
python scripts/web_ui.py --host 0.0.0.0 --port 5000
```

Open `http://<host>:5000` in a browser. If email login is configured
(`PIPELINE_ALLOWED_EMAILS` environment variable + either `RESEND_API_KEY` or
SMTP credentials), you'll be asked to sign in via a one-time magic link sent
to your email; the link is valid for 30 minutes and the resulting session
lasts 7 days. If no email is configured, the console will print the magic
link directly so you can copy it into the browser.

### 2. Choosing an input method

The main page offers three ways to start an analysis:

**Method 1 — GEO accession (one-click)**
Paste a `GSExxxxxx`, `PRJNAxxxxxx`, or `SRPxxxxxx` accession and a core count,
then submit. The pipeline fetches sample metadata automatically and tries to
detect case/control (tumor/normal) labels from the sample descriptions —
review the detected labels before starting, since automatic detection can be
ambiguous for unusually-named samples.

**Method 2 — Manual SRR list / CSV upload**
If automatic GEO detection doesn't work for your dataset (common for
`PRJNA`/`SRP` accessions, whose metadata endpoints are sometimes blocked on
restricted networks), enter SRR IDs and conditions by hand, or upload a CSV
with `srr_id,condition` columns (a `patient_id` column can be added for
paired tumor/normal designs). A worked example table is shown on the page.

**Method 3 — Local FASTQ files**
If you already have FASTQ files on the server (e.g. from a private
sequencing run), point the Web UI at the directory; it auto-detects paired
files (`_1/_2`, `_R1/_R2`, `_R1_001/_R2_001` conventions) and lets you assign
conditions per sample in a table. The pipeline symlinks these into its
working directory and skips the download step entirely.

Before picking a dataset, it's worth reviewing the built-in **Dataset
Selection Guide** card on the main page — it walks through the read-length,
library-prep, and sample-size criteria that most affect result quality (see
also [Choosing a Good Dataset](#choosing-a-good-dataset) below).

### 3. Configuring the analysis

Three configuration steps appear after you pick an input method:

- **Step 1 — circRNA detection tools**: CIRIquant only, DCC only, or both
  (recommended; the pipeline requires agreement between both tools by
  default, improving specificity at some cost to sensitivity).
- **Step 2 — differential expression**: all three methods (`edgeR_ciriquant`,
  `DESeq2`, `limma-voom`) always run; you choose which one the report
  displays by default (you can switch between them interactively in the
  finished report without re-running anything).
- **Step 3 — advanced parameters**: minimum BSJ read count, coordinate
  tolerance between tools (`slop`), pseudo-circRNA filtering ratio, FDR/log2FC
  cutoffs, and thread count.

Saving the configuration writes both the global `config.yaml` and a permanent
per-project snapshot at `config/projects/{GSE_ID}.yaml`, so re-launching the
same accession later reuses its own settings rather than whatever project ran
most recently.

### 4. Monitoring progress

Submitting an analysis takes you to a **status page** that auto-refreshes
every 5 seconds, showing a progress bar and a grid of pipeline stages (each
turns green as its output appears). If you submit multiple datasets, they're
placed in a **job queue** (`/queue`) and run one at a time — you'll get an
email when a job starts and when it finishes (or fails), if email is
configured.

### 5. Reading the report

The finished report (`report.html`) is fully self-contained — download it or
open it directly in the browser. It includes, top to bottom:

- **Sample overview** — read counts, trimming stats, and QC summary per
  sample (tumor/normal color-coded).
- **DE method switcher** — toggle between `edgeR_ciriquant` / `DESeq2` /
  `limma-voom`; every chart and table below updates in place.
- **PCA and volcano plots** (interactive Plotly charts).
- **Clustering heatmap** — all significant circRNAs, hierarchically clustered
  with a real dendrogram; always expanded, no need to click to reveal it.
- **Top differentially-expressed circRNA tables**, split into up/down-regulated,
  with gene annotation, exon structure, and clickable circBase links.
- **Type I / Type II classification** (only for the `edgeR_ciriquant` method)
  — distinguishes circRNA-specific regulation (Type I) from cases where the
  linear mRNA also changed significantly (Type II).
- **3-method Venn diagram** — click any region to see which circRNAs are
  shared or unique across the three DE methods.
- **Biomarker candidate table**, ranked by a composite score (significance +
  fold-change + detection confidence + circBase novelty + predicted miRNA/RBP
  binding partners); the display count is adjustable.
- **Isoform switching table** — genes with multiple circRNA isoforms whose
  relative usage shifts between conditions.
- **Per-circRNA detail view** — click any circRNA ID to open a modal with its
  circular exon structure, predicted miRNA sponge sites, RBP binding sites,
  a mini volcano plot, and a zoomed-in clustering-heatmap view of its nearest
  neighbors by expression pattern.

Every table can be sorted by clicking its headers and exported to CSV; the
whole report can be printed or saved as PDF from the toolbar at the top.

### 6. Cross-dataset comparison

Once you've analyzed more than one dataset, visit `/cross_dataset` from the
navigation bar to see which circRNAs are differentially expressed in more
than one dataset — including an option to compare only datasets that share
the same cancer type/tissue, rather than mixing everything together.

---

## Command-Line Usage

For scripting, HPC job schedulers, or CI, you can run the pipeline directly
with Snakemake instead of the Web UI:

```bash
conda activate ciriquant
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --resources mem_gb=300 \
    --keep-going \
    --rerun-incomplete
```

`scripts/agent.py` offers a thinner CLI wrapper around the same workflow:

```bash
python scripts/agent.py --gse GSE113230 --cores 8
python scripts/agent.py --gse GSE113230 --cores 8 --dry-run     # preview the DAG
python scripts/agent.py --gse GSE113230 --cluster "sbatch ..." --jobs 50
```

> **Note:** `agent.py`'s GEO metadata lookup depends on `pysradb`, which may
> not be installed in every environment. If it's unavailable, prepare
> `metadata/library_info.csv` and `metadata/sample_groups.csv` yourself (or
> via the Web UI's manual CSV upload, which writes exactly these files) and
> invoke `snakemake` directly as shown above.

---

## Configuration Reference

`config.yaml` (or a per-project snapshot under `config/projects/`) controls
every parameter of a run:

```yaml
project_id: GSE113230                 # used to name output directories
metadata:   metadata/library_info.csv # SRR/sample metadata
groups:     metadata/sample_groups.csv # condition labels

raw_dir:     /path/to/raw_fastq
trimmed_dir: /path/to/trimmed
results_dir: /path/to/results

genome:
  fasta:        /path/to/hg19.fa
  gtf:          /path/to/genes.gtf
  bwa_index:    /path/to/hg19.fa
  hisat2_index: /path/to/hisat2_index
  star_index:   /path/to/star_index
  species:      hg19

ciriquant_config: config/ciriquant.yaml

consensus:
  tools:              [ciriquant, dcc]  # or a single tool
  min_tools:           2                # required agreeing tools
  slop:                10               # coordinate tolerance (bp)
  min_bsj_reads:        2
  max_junction_ratio:   1.0             # pseudo-circRNA QC (CIRIquant-only loci)

de:
  method:              edgeR_ciriquant  # report's default display method
  fdr_cutoff:          0.05
  log2fc_cutoff:       1.0
  de_sig_by:           pvalue           # pvalue (nominal) or padj (BH-corrected)
  tumor_label:         tumor
  normal_label:        normal

threads: 8
```

Key parameters worth understanding before your first run:

- **`consensus.tools`** — using both `ciriquant` and `dcc` gives higher-confidence
  calls (both tools must agree within `slop` base pairs) at the cost of some
  sensitivity; single-tool mode is used automatically as a fallback when one
  tool's detection count is drastically lower than the other's (e.g. short
  reads that STAR's chimeric alignment handles poorly).
- **`de.de_sig_by`** — small studies (3 vs. 3 samples is common for circRNA
  studies) rarely survive Benjamini-Hochberg correction across thousands of
  tests; `pvalue` (nominal, uncorrected) is the common practical choice for
  exploratory work, but should be reported as such in any resulting
  publication.

---

## Pipeline Architecture

```
SRA/GEO accession
      │
      ▼
 [download]  aria2c (S3) → ascp → prefetch (fallback chain)
      │
      ▼
   [QC/trim]  FastQC + fastp + MultiQC
      │
      ├─────────────────────────┐
      ▼                         ▼
 [CIRIquant]              [STAR chimeric] ──▶ [DCC]
 (HISAT2 + BWA)                 │
      │                          │
      └──────────┬───────────────┘
                 ▼
        [consensus filter]
        coordinate-tolerant voting between tools
        + confidence scoring
                 │
                 ▼
        [count matrix + annotation]
        BSJ/FSJ counts, host gene, circBase ID
                 │
                 ▼
   [differential expression]
   edgeR_ciriquant / DESeq2 / limma-voom (all three, always)
                 │
        ┌────────┴────────┐
        ▼                 ▼
 [biomarker ranking]  [isoform switching]
        │                 │
        └────────┬────────┘
                 ▼
         [HTML report]
```

---

## Choosing a Good Dataset

Not every public RNA-seq dataset is suitable for circRNA analysis. Before
committing compute time to a large download, check:

| Factor | Look for | Avoid |
|--------|----------|-------|
| Read length | ≥ 100 bp paired-end | < 50 bp, single-end |
| Library prep | Total RNA (rRNA-depleted) or RNase R-enriched | poly-A selected (circRNAs lack poly-A tails and are barely detectable) |
| Replicates | ≥ 3 per group (≥ 5 preferred) | n = 1 or 2 |
| Design | Paired tumor/adjacent-normal from the same patient | unrelated case/control cohorts (less statistical power) |
| Sequencing depth | ≥ 50M reads/sample | < 30M reads/sample |

Squamous-cell carcinomas and other cancer types with globally low circRNA
expression can pass every criterion above and still yield very few testable
circRNAs after expression filtering — this becomes apparent only after
running the pipeline, at which point switching the report's default DE
method to `limma-voom` (more robust to sparse count matrices) is a reasonable
fallback.

---

## Container Deployment (Docker / Singularity)

A pre-built image with every dependency (CIRIquant, DCC, STAR, HISAT2, BWA,
R/Bioconductor, Snakemake) is published to Docker Hub, and can be converted
to a Singularity/Apptainer image for HPC clusters that don't allow Docker.

### 1 — Build and test locally (optional — only needed if you're modifying the image)

```bash
bash containers/build_and_deploy.sh
```

This builds `circrna-pipeline:1.0.1` from the `Dockerfile`, runs a smoke test
(CIRIquant, Snakemake, R/edgeR), and pushes it to Docker Hub.

### 2 — Pull the image on your HPC cluster

```bash
singularity pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.1
# or, with Apptainer (identical syntax):
apptainer pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.1
```

> Your cluster must have Singularity/Apptainer with user namespaces enabled.
> If pulling fails with a namespace-related permission error, ask your
> cluster administrator to run `echo 10000 > /proc/sys/user/max_user_namespaces`,
> or use a setuid-installed Singularity build.

### 3 — Run the pipeline with Singularity

Set `use_container: true` in `config.yaml` (this switches CIRIquant to
`config/ciriquant_container.yaml`, which uses in-container tool names instead
of absolute host paths), then:

```bash
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --use-singularity \
    --singularity-args "--bind /path/to/your/data:/path/to/your/data" \
    --keep-going --rerun-incomplete
```

---

## Troubleshooting

- **Web UI page auto-translated into the wrong language by your browser** —
  the pages set `translate="no"`, but some browsers still offer a translation
  prompt on first visit; decline it once and it won't reoffer.
- **A pipeline run appears to hang** — check the status page's log viewer
  first; long CIRIquant runs (multiple hours for large total-RNA samples) are
  normal, not necessarily stuck. If a run was interrupted (e.g. the process
  was killed), you may need `snakemake --unlock` before restarting.
- **`ConfigError` from CIRIquant about missing tools** — CIRIquant checks
  tool paths literally; `config/ciriquant.yaml`'s `tools:` section needs
  full absolute paths (e.g. `/opt/conda/envs/ciriquant/bin/bwa`), not bare
  command names.
- **Very few circRNAs pass expression filtering** — this is common for
  squamous-cell carcinomas, short-read libraries, or small/unbalanced sample
  sizes; try switching the report's DE method to `limma-voom`, which is more
  robust to sparse count matrices than `edgeR`/`DESeq2` defaults.
- **PRJNA/SRP accession metadata lookup fails** — some networks block the
  NCBI eUtils endpoints those accessions rely on; use the Web UI's manual
  CSV upload instead, after fetching metadata from a machine that does have
  access (`python scripts/download_geo.py --gse PRJNAxxxxxx`).

---

## Project Structure

```
circDEX/
├── config.yaml                  # active run configuration
├── config/
│   ├── ciriquant.yaml           # CIRIquant tool paths
│   └── projects/{GSE_ID}.yaml   # per-project saved configuration
├── workflow/
│   ├── Snakefile
│   └── rules/                   # download / QC / detection / DE rule definitions
├── scripts/
│   ├── agent.py                 # CLI entry point
│   ├── web_ui.py                # Flask Web UI
│   ├── consensus_filter.py      # dual-tool circRNA consensus voting
│   ├── analysis.R               # differential expression (3 methods)
│   ├── rank_biomarkers.py       # composite biomarker scoring
│   ├── generate_report.py       # interactive HTML report builder
│   └── templates/               # Web UI + report HTML templates
├── metadata/                    # per-project sample metadata (auto-generated)
├── envs/circrna.yaml            # conda environment definition
└── containers/                  # Docker/Singularity build scripts
```

---

## Methodology Notes

- **Consensus detection**: coordinate-tolerant voting between CIRIquant and
  DCC output (default ±10 bp), following the general consensus-detection
  rationale established by CirComPara2 (Gaffo et al. 2022) and Hansen (2018).
  A confidence score (log-scaled read support × coordinate agreement,
  averaged per supporting tool) is reported per circRNA but is a heuristic
  ranking signal, not a statistical probability.
- **`edgeR_ciriquant`**: tests whether the BSJ/FSJ ratio (back-splice vs.
  linear-splice junction usage) changes between conditions, using FSJ counts
  as a normalization offset — this distinguishes circRNA-specific regulation
  (Type I) from cases where the change is really in overall gene expression
  and the circRNA simply followed along (Type II).
- **Isoform switching**: significance is assessed with within-gene
  Benjamini-Hochberg correction (DEXSeq-style), not a global correction
  across all isoforms — appropriate given that switching is a gene-local
  phenomenon.

For dataset-specific analysis notes, known issues, and their resolutions,
see `CLAUDE.md` in the repository root.
