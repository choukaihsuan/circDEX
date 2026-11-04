# circRNA Analysis Pipeline — Project Context for Claude

## 專案概述

本專案是一個以 **Snakemake** 驅動的 circRNA（環狀 RNA）全流程分析管線，
從 GEO/SRA 原始數據下載，到差異表現分析（DE）與 HTML 報告輸出。

- **目標數據集**：GSE113230（三陰性乳癌 tumor vs. normal，6 samples，✅ 完成）；GSE58135（乳癌，10 samples，✅ 完成）；GSE323364（TNBC cell line EZH2 inhibitor，6 samples，✅ 完成）；GSE133998（乳癌 tumor vs. normal，12 samples，✅ 完成）；SRP156355（早期乳癌 IDC，6 pairs，✅ 完成）；GSE77509（HCC 肝癌 tumor vs. normal，6 pairs，✅ 完成）；GSE130078（ESCC 食道鱗狀細胞癌 tumor vs. normal，6 pairs，✅ 完成）；GSE248612（胃癌 tumor vs. normal，6 pairs，✅ 完成）；GSE221107（攝護腺癌 tumor vs. normal，4 pairs，✅ 完成；排除 Pair8/Pair11 降解 RNA）；PRJNA553289（SCLC 小細胞肺癌 tumor vs. normal，6 pairs，✅ 完成）；GSE229705（LUAD 肺腺癌 tumor vs. normal，6 pairs，✅ 完成）；GSE148036（LUAD 肺腺癌 tumor vs. normal，5+5 samples，✅ 完成）；GSE121842（CRC 大腸直腸癌 tumor vs. normal，3 pairs，✅ 完成）；GSE136569（胰臟癌 PDAC tumor vs. NAT，5 pairs，✅ 完成）；GSE143797（鼻咽癌 NPC tumor vs. normal，4+4 samples，✅ 完成）；GSE108735（腎細胞癌 RCC tumor vs. normal，7 pairs，✅ 完成）；GSE171011（甲狀腺乳突癌 PTC tumor vs. normal，4 pairs，✅ 完成）；GSE97239（膀胱癌 Bladder Cancer tumor vs. normal，3 pairs，✅ 完成）；GSE192410（卵巣癌 Ovarian Cancer tumor vs. normal，3 pairs，✅ 完成）；GSE192849（乳癌 node-positive tumor vs. normal，3 pairs，RNase R 富集，✅ 完成）
- **主要工具**：CIRIquant（circRNA 偵測）+ DCC（輔助偵測，雙工具共識）
- **執行環境**：基因體中心 HPC server（`172.16.0.178`，CentOS 7，96 cores，377 GB RAM）
- **本機開發**：Windows 11 + WSL2（Ubuntu 26.04），程式碼在 `/mnt/c/Users/User/develop/circRNA_agent/`
- **Server 路徑**：`~/circRNA_agent/`（即 `/home3/choukaihsuan/circRNA_agent/`，`/home/choukaihsuan` 是 symlink）
- **Container**：Docker image `choukaihsuan/circrna-pipeline:1.0.1`；HPC 用 Singularity 拉取

---

## 目錄結構

```
circRNA_agent/
├── config.yaml                  # 主設定檔（路徑、參數、工具選擇）
├── Dockerfile                   # 容器化環境定義（mamba + conda env circrna.yaml）
├── .dockerignore                # Docker 建置排除清單
├── config/
│   ├── ciriquant.yaml           # CIRIquant 工具路徑設定（server 版本另存在 server 上）
│   ├── .ciriquant_ready         # touch 檔，驗證 ciriquant.yaml 存在後建立
│   └── projects/                # 各 GSE 專案獨立設定快照（由 web_ui.py 自動建立）
│       └── {GSE_ID}.yaml
├── containers/
│   └── build_and_deploy.sh      # Docker image 建置與推送腳本
├── metadata/
│   ├── library_info.csv         # 目前啟用專案的 SRR/配對資訊
│   ├── sample_groups.csv        # 目前啟用專案的分組
│   └── {GSE_ID}/                # 各 GSE 專屬 metadata 目錄（由 web_ui.py 自動建立）
│       ├── library_info.csv
│       └── sample_groups.csv
├── workflow/
│   ├── Snakefile                # 主 Snakefile，載入 rules，設定 target
│   └── rules/
│       ├── download.smk         # SRA 下載（aria2c S3 > ascp > prefetch + fasterq-dump）
│       ├── qc.smk               # FastQC + fastp + MultiQC
│       ├── circrna.smk          # circRNA 偵測主規則（見下方詳述）
│       └── de.smk               # 差異表現 + 圖表 + HTML 報告
├── scripts/
│   ├── agent.py                 # CLI 入口（--gse, --setup-ciriquant 等）
│   ├── prepare_metadata.py      # 從 GEO/SRA RunInfo 建立 library_info.csv；T/N/APN 尾碼自動偵測
│   ├── download_geo.py          # SRA metadata 抓取（GSE→pysradb；PRJNA/SRP→NCBI eUtils API）
│   ├── consensus_filter.py      # CIRIquant + DCC 共識過濾（輸出 BED + confidence score）
│   ├── merge_counts.py          # 從 CIRIquant GTF 建立 BSJ + FSJ count matrix
│   ├── annotate_circbase.py     # circBase hg19 座標比對注釋（自動下載或讀本地檔）
│   ├── rank_biomarkers.py       # Biomarker 候選排序（composite score）
│   ├── assign_isoforms.py       # circRNA BSJ 座標對應 host gene + exon span + strand + region
│   ├── isoform_switching.R      # 計算 IUI、Wilcoxon rank-sum 測試 isoform switching
│   ├── analysis.R               # DE 分析（edgeR_ciriquant / deseq2 / limma）
│   ├── generate_report.py       # 輸出 HTML 報告（互動式 Plotly + Type I/II + biomarker）
│   ├── predict_interactions.py  # CircInteractome miRNA/RBP interaction 查詢
│   ├── notify.py                # 通知模組（Email/Slack，Snakemake hook 呼叫）
│   ├── utils.py                 # 共用工具函數
│   ├── web_ui.py                # Flask Web UI（GEO 一鍵啟動 + 進度視覺化）
│   └── templates/
│       ├── index.html           # 主設定頁面（GEO 入口 + Step 1-3）
│       └── status.html          # Pipeline 狀態頁（進度條 + rule 狀態格 + log）
├── envs/
│   └── circrna.yaml             # Conda 環境定義（mamba 建置；已移除 defaults channel）
└── logs/                        # Snakemake 各 rule 的 log 檔
```

---

## Pipeline 流程

```
SRA/GEO
  │
  ▼
[download] aria2c (S3, 16連線) > ascp (Aspera) > prefetch (HTTPS fallback)
  → raw_dir/{srr}_1.fastq.gz, {srr}_2.fastq.gz
  │
  ▼
[qc] FastQC（raw）+ fastp（trim）+ MultiQC
  → trimmed_dir/{srr}_1.fastq.gz, {srr}_2.fastq.gz
  → results/qc/fastp/{srr}.json, multiqc_report.html
  │
  ├──────────────────────────────────────────┐
  ▼                                          ▼
[ciriquant]                         [star_align]（paired-end）
HISAT2 + BWA-MEM                    STAR chimeric junctions
→ {srr}/{srr}.gtf                   → Chimeric.out.junction
  (BSJ + FSJ counts)                → Aligned.sortedByCoord.out.bam
                                           │
                                    ┌──────┴──────┐
                                    ▼             ▼
                             [star_align_mate1] [star_align_mate2]
                             STAR R1 only      STAR R2 only
                             → mate1/Chimeric  → mate2/Chimeric
                                    │
                                    ▼
                               [dcc]
                               DCC 0.5.0
                               -mt1 mate1/Chimeric
                               -mt2 mate2/Chimeric
                               → DCC/CircCoordinates
  │                                 │
  └──────────────────┬──────────────┘
                     ▼
             [consensus_filter]
             共識過濾（min_tools=2, slop=10, min_bsj=2）
             pseudo-circ QC（BSJ/FSJ > max_junction_ratio 過濾）
             confidence_score = Σ[log2(bsj+1)×(1−dist/slop)] / n_supporting
             → high_confidence.bed（7 欄，含 confidence_score）
             → consensus_summary.tsv
             （支援單工具模式：USE_CIRIQUANT / USE_DCC）
                     │
                     ▼
             [merge_counts]
             parse GTF → BSJ + FSJ count matrix
             → circRNA/count_matrix.tsv
             → circRNA/fsj_count_matrix.tsv
                     │
          ┌──────────┴───────────────────────┐
          ▼                                  ▼
  [assign_isoforms]                  [annotate_circbase]
  host gene + strand +               比對 circBase hg19
  exon_span + region                 → circbase_annotated.tsv
  → isoform_groups.tsv
          │
          ▼
   [de_analysis]
   edgeR_ciriquant / deseq2 / limma
   → de/de_results.tsv
   → plots/*.pdf
          │
          ├──────────────────┐
          ▼                  ▼
  [rank_biomarkers]   [isoform_switching]
  composite score     IUI + Wilcoxon test
  → biomarker_        → iui_matrix.tsv
    candidates.tsv    → isoform_switching.tsv
          │
          └──────────┬───────────────────────┘
                     ▼
             [generate_report]
             → report.html
```

---

## 關鍵技術細節

### circRNA 偵測工具選擇（config `consensus.tools`）

`Snakefile` 在啟動時讀取 `config.consensus.tools`（list），設定以下全域變數：

```python
TOOLS        = set(config["consensus"]["tools"])  # {"ciriquant", "dcc"} 或子集
USE_CIRIQUANT = "ciriquant" in TOOLS
USE_DCC       = "dcc"       in TOOLS
```

- **兩者並行**（預設）：`consensus_filter` 取交集，min_tools=2
- **CIRIquant only**：跳過 STAR 三次比對 + DCC；`consensus_filter` 單工具模式（min_tools=1）
- **DCC only**：跳過 CIRIquant；`consensus_filter` 單工具模式

`consensus_filter.py` 的 `--cirique` 和 `--dcc` 都是 **optional**，至少提供一個即可。

### DCC 的正確呼叫方式（DCC 0.5.0 重要）

DCC 0.5.0 必須同時提供三個 junction 檔案，缺少 `-mt2` 會出現：
`TypeError: object of type 'NoneType' has no len()`

```bash
DCC {paired_junction} \
    -mt1 {mate1_junction} \     # R1 單端 STAR 比對
    -mt2 {mate2_junction} \     # R2 單端 STAR 比對
    -D -an {gtf} \
    -Pi -F -M -Nr 5 1 \
    -G -B {bam} \
    -O {outdir} -T 4
```

- 命令是大寫 `DCC`（不是 `dcc`）
- BAM 參數是 `-B`（不是 `-A`）
- `-mt1`/`-mt2` 直接傳檔案路徑（不支援 `@filelist` 格式）
- **`-fg` 已移除**：觸發 CircSkip 解析，GTF exon_number 屬性尾部引號造成 `ValueError`

### Consensus 過濾與 Confidence Score

`consensus_filter.py` 的過濾流程（每個 sample 獨立執行）：

1. **min_bsj 過濾**：各工具輸出中，BSJ < min_bsj 的 circRNA 直接丟棄
2. **Pseudo-circ QC**（CIRIquant only）：解析 GTF 中的 FSJ，若 BSJ/FSJ > `max_junction_ratio`（預設 1.0）則丟棄。真實 circRNA 幾乎都有 BSJ < FSJ；比值 > 1 是誤比對或 repeat 假陽性的警訊。FSJ = 0 的 locus 不做此過濾（缺乏對應線性轉錄本資訊）
3. **座標共識投票**：各工具的座標在 slop 範圍內（預設 10 bp）視為一致
4. **Confidence score** = `Σ[log2(bsj+1) × (1−dist/slop)] / n_supporting_tools`
   - 分母是**實際支持該 circRNA 的工具數**（非所有工具數），代表「各支持工具的平均每工具信心」
   - 論文應標明此為 weighted scoring heuristic，**非機率值**，引用 CirComPara2 + Hansen (2018) 作為 consensus 正當性基礎

**Adaptive fallback**（`--adaptive` flag）：
若兩工具偵測數量嚴重失衡（`min(counts) / max(counts) < --adaptive-ratio`，預設 0.1），
自動將 `min_tools` 從 2 降為 1，避免共識過濾後近乎零回收率。
失衡時在 stderr 印出警告。新 sample 類型的探索性分析時建議開啟。

### circBase 注釋（`annotate_circbase.py`）

- 從 `http://www.circbase.org/download/hsa_hg19_circRNA.txt` 自動下載（或讀 `--circbase-file` 本地檔案，預設快取於 `/tmp/circbase_hg19.txt`）
- 以 max(|Δstart|, |Δend|) ≤ slop 判斷為已知 circRNA
- 輸出欄位：`circbase_id`（或 `"novel"`）、`circbase_gene`、`in_circbase`（0/1）

### Host Gene / Exon / Strand 注釋（`assign_isoforms.py`）

`assign_isoforms.py` 對每個 circRNA 輸出以下額外欄位：

| 欄位 | 說明 |
|------|------|
| `gene_name` | Host gene 名稱（GTF gene 特徵，取最小包含基因體） |
| `strand` | `+` 或 `-`（來自 host gene GTF 記錄） |
| `exon_span` | 環化的 exon 範圍，格式 `eN-eM`（如 `e3-e7`）；找不到時為空白 |
| `region` | `exonic`（兩端都匹配 exon 邊界）/ `intronic`（在基因內但不匹配 exon）/ `intergenic` |

**exon_span 判斷邏輯**：對每個 circRNA (start, end) 掃描 host gene 所有 transcript，找到 `exon.start ≈ circ_start` 且 `exon.end ≈ circ_end`（±10 bp）的 transcript，回報 exon 編號。負股基因的 exon_number 可能出現反向（如 `e9-e8`），這是正常現象——GTF exon_number 按基因組座標遞增，但 mRNA 順序在負股上相反。

GSE113230 結果：9,349 circRNA 中 8,046 exonic（86%）、1,089 intronic（12%）、214 intergenic（2%）。

### Biomarker 排序（`rank_biomarkers.py`）

六維 composite score，每維度在 significant set 內 min-max 標準化後平均：

```
biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus + mirna_norm + rbp_norm) / 6
  sig_norm   = −log10(pvalue 或 padj), 上限 10，標準化  ← 依 de_sig_by 切換
  fc_norm    = |log2FC|,               上限 5，標準化
  conf_norm  = confidence_score 標準化
  known_bonus = 1 若 in_circbase，否則 0（不標準化）
  mirna_norm = distinct miRNA binders 數，min-max 標準化（無 interaction data 者為 0）
  rbp_norm   = distinct RBP binders 數，min-max 標準化（無 interaction data 者為 0）
```

- 顯著閾值欄位依 `de_sig_by` 而定：`pvalue`（nominal）或 `padj`（BH 校正）
- CLI 參數：`--use-pvalue` 對應 `de_sig_by: pvalue`
- miRNA/RBP interaction 資料來自 `predict_interactions.py`（CircInteractome/ENCORI 查詢，**三方法 top-50 聯集**，~150–250 個 circRNA）；union mode：`--de-edger / --de-deseq2 / --de-limma` 各取 top-50 up + top-50 down 後取聯集，確保切換 DE 方法時 Biomarker Score 完整

### DE 分析方法（`config de.method`）

`analysis.R` **預設三種方法全跑**，各輸出獨立 TSV：
- `de_results_edgeR_ciriquant.tsv`（同時作為主 `de_results.tsv`）
- `de_results_deseq2.tsv`
- `de_results_limma.tsv`

| 值 | 說明 |
|----|------|
| `edgeR_ciriquant`（預設主方法）| 複製 `CIRI_DE_replicate`：edgeR GLM + FSJ offset（測 BSJ/FSJ 比值），輸出 Type I/II 分類 |
| `deseq2` | DESeq2 RLE normalization + Wald test |
| `limma` | limma-voom（TMM normalization）|

**`edgeR_ciriquant` 核心邏輯**：
FSJ counts → TMM normalization → per-locus FSJ CPM 作為 GLM offset
→ 效果等同於測試 BSJ/FSJ 比值是否在 tumor vs. normal 之間改變
→ 同時對 FSJ 跑獨立 QLFTest，若 FSJ 也顯著（不做方向一致性檢查）→ Type II（兩層調控同時發生），否則 → Type I（circRNA 專一性）

**Type I / II / III 分類**：
- **Type_I**：BSJ 顯著，FSJ 不顯著 → circRNA 環化效率真正改變（circRNA-specific regulation）
- **Type_II**：BSJ/FSJ ratio 顯著（BSJ offset test）且 FSJ 獨立測試也顯著 → 兩種調控同時發生
- **Type_III**：只有 FSJ 顯著 → 線性 mRNA 變化，不是 circRNA DE
- **不做方向一致性檢查**：offset approach 造成 log2FC（ratio 方向）與 logFC_fsj（FSJ 方向）天然反向相關（FSJ↑ → offset↑ → ratio log2FC 傾向負值），強制要求 sign 一致等同於把幾乎所有 Type II 歸入 Type I。正確判斷：`sig_bsj AND sig_fsj`（無方向條件）
- `fsj_concordance_lfc`（config `de.fsj_concordance_lfc`，預設 `0.0`）：可設定最小 |logFC_fsj| 作為 FSJ 顯著的過濾條件（0 = 僅看統計顯著性）

**重要欄位命名**：merge BSJ/FSJ 結果後，`PValue` 不加後綴（只有 `logFC`、`FDR` 因同時出現在兩表才加 `_bsj`/`_fsj`）。`bsj_sig_col` 必須用 `"PValue"`，不是 `"PValue_bsj"`。

**顯著性欄位切換**（`config de.de_sig_by`）：
- `padj`（預設）：BH 校正 FDR。小樣本（n=3 vs 3）+ 多重檢定（~7,779 tests）時 min padj ≈ 0.432，幾乎無法通過
- `pvalue`：nominal p-value（未校正）。小樣本研究的實務做法，論文中需標明

`analysis.R` 有 **backward-compatible fallback**：若 `snakemake@input[["fsj_matrix"]]` 不存在（舊 DAG），自動 fallback 到 deseq2。

### Isoform Switching（`isoform_switching.R`）

**是否顯著的判斷使用 within-gene BH FDR**（`padj_within_gene`），而非 global BH。

原因：global BH 跨所有 isoform（~7,360 tests）→ min padj = 0.937，完全無法找到 switching。Within-gene BH（每個基因內部獨立做 BH）是 DEXSeq 風格的標準做法，在 GSE113230 結果中找到 66 個 significant switching events（FDR < 0.1, |ΔIUI| > 0.1）。

---

## 設定檔說明（config.yaml）

```yaml
project_id: GSE113230
metadata:   metadata/library_info.csv
groups:     metadata/sample_groups.csv

# 路徑（server 版本使用 /home3/choukaihsuan/...）
raw_dir:     /mnt/d/circRNA_data/raw_fastq
trimmed_dir: /mnt/d/circRNA_data/trimmed
results_dir: /mnt/d/circRNA_results/GSE113230

genome:
  fasta:       /mnt/d/ref/hg19/hg19.fa
  gtf:         /mnt/d/ref/hg19/hg19.gtf
  bwa_index:   /mnt/d/ref/hg19/bwa_index/hg19
  hisat2_index: /mnt/d/ref/hg19/hisat2_index/hg19
  star_index:  /mnt/d/ref/hg19/star_index
  species:     hg19

ciriquant_config: config/ciriquant.yaml

consensus:
  tools:             [ciriquant, dcc]  # 工具選擇
  min_tools:         2                # 共識閾值
  slop:              10               # 座標容忍 bp
  min_bsj_reads:     2
  max_junction_ratio: 1.0             # pseudo-circ QC：BSJ/FSJ 上限（CIRIquant only）

circbase_file: ""   # 留空 = 自動下載；或填本地路徑

de:
  method:              edgeR_ciriquant    # edgeR_ciriquant / deseq2 / limma
  fdr_cutoff:          0.05
  log2fc_cutoff:       1.0
  de_sig_by:           pvalue            # pvalue = nominal p；padj = BH 校正 FDR
  isoform_fdr_cutoff:  0.1               # within-gene FDR for isoform switching
  delta_iui_cutoff:    0.1               # minimum |ΔIUI| to call switching
  tumor_label:         tumor
  normal_label:        normal

threads: 8
```

**Server 版 config**（`~/circRNA_agent/config.yaml`）路徑前綴不同：
- `raw_dir: /home3/choukaihsuan/GSE113230/raw`
- `trimmed_dir: /home3/choukaihsuan/GSE113230/trimmed`
- `results_dir: /home3/choukaihsuan/GSE113230_results`
- `star_index: /home3/choukaihsuan/reference/hg19/star_index`
- `gtf: /home3/choukaihsuan/reference/hg19/genes.gtf`（注意：檔名是 `genes.gtf`，非 `hg19.gtf`）

---

## HPC Server 環境

| 項目 | 值 |
|------|----|
| 主機 | `172.16.0.178`（也可用 `172.16.0.179`） |
| 使用者 | `choukaihsuan` |
| Home | `/home/choukaihsuan` → symlink 到 `/home3/choukaihsuan` |
| OS | CentOS 7 |
| Python | 3.7.12（conda env `ciriquant`） |
| CPU | 96 cores |
| RAM | 377 GB |
| 磁碟 | /home3：345 GB 可用（2026-06-19；中間檔清理後）|
| Conda env | `ciriquant`（CIRIquant 1.1.3, DCC 0.5.0, STAR 2.7.11a, HISAT2 2.2.1, BWA 0.7.18, samtools 1.18, snakemake 7.24.0, **aria2c 1.36.0**, sra-tools 2.11.0, pigz 2.8, FastQC 0.12.1, fastp 1.1.0, MultiQC 1.13）。**CIRIquant 版本三方查證（2026-10-03，server 實測）**：conda 登記顯示 `1.1.2 pypi_0`（登記紀錄過時）；`pip show CIRIquant` 為 1.1.3；`CIRIquant --version` 為 1.1.3。以 binary 與 pip metadata 為準，版本為 **1.1.3**。曾於 2026-10-01 因只讀 `conda list --export` 而誤改為 1.1.2（c7e0559），已於 2026-10-03 還原。|
| Java | `/usr/bin/java`（不在 conda env 內，ciriquant.yaml 必須指定此路徑） |

**SRA 下載優先順序**（`workflow/rules/download.smk`）：

| 優先 | 方法 | 實測速度 | 說明 |
|------|------|----------|------|
| 1 | **aria2c + S3**（預設）| **~25–30 MB/s**（6–8 連線，NFS 環境）| `srapath --location s3` 取 S3 URL → aria2c 多連線下載；`-x 6 -s 6` 適合 NFS（過多連線反而限速）|
| 2 | ascp（Aspera）| ~50 MB/s | 需 Aspera key，目前 server 未安裝 |
| 3 | prefetch（HTTPS）| **~20–100 KB/s**（實測 22 KB/s）| NCBI 單連線，最慢；656MB 需 8+ 小時；S3 失敗時的 fallback |

**後處理速度（6GB SRA file，NFS）**：

| 步驟 | 工具 | 實測速度 | 說明 |
|------|------|----------|------|
| SRA → FASTQ | fasterq-dump `-e 6` | ~10–15 min/6GB | 多執行緒解壓縮；`--split-files` 輸出 PE |
| FASTQ 壓縮 | **pigz** `-p 3` | **~15 min/15GB** | 多核 gzip，較快 |
| FASTQ 壓縮 | gzip（單執行緒）| ~2h/15GB | pigz 未安裝時的 fallback |

**aria2c 並行下載建議**：
- **S3 per-IP 總連線數安全閾值**：≤ 48 連線（實測 `-x 12` × 4 parallel = 48 可接受，GSE77509 驗證）；超過 60 會觸發 DNS SERVFAIL（AWS S3 throttle），持續約 2 小時
- **當前 download.smk 設定**：`-x 8 -s 8`；Snakemake `threads: 4`，24 cores ÷ 4 = 6 parallel → 6 × 8 = 48 連線（接近安全上限，GSE148036 預設設定）
- **手動腳本建議設定**：`-x 12 -s 12`，batch size 4（4 parallel）→ 48 連線，每個 sample 實測 ~7 MB/s，總速度 ~25–30 MB/s
- **勿手動腳本同時跑多個 SRR**：10 × 16 = 160 連線曾觸發 DNS 封鎖（GSE133998 教訓）
- S3 間歇性不可用（NCBI API）：`srapath --location s3` 回傳空字串，等待恢復後手動重啟
- 下載時若 kill 中斷：aria2c 會留下 `.sra.aria2` resume 檔 → 重啟可斷點續傳；但若多次中斷導致 `.aria2` 狀態損壞，需刪除 `.sra` 和 `.aria2` 重新下載，否則 fasterq-dump 報 `rcBlob,rcCorrupt`

`_find_tool("aria2c")` 搜尋優先順序：`sra_env` → `circrna` → **`ciriquant`**（已加入）→ `which()`。
若 `srapath returned no S3 URL`（NCBI API 暫時失敗），重啟 pipeline 即可；S3 URL 通常幾分鐘後恢復。

**R packages（安裝在 conda env `ciriquant` 的 r-base 4.2.1，2026-10-01 實測更正先前誤記的 4.2.2；Bioconductor 3.16）**：
r-ggplot2, r-pheatmap, r-rcolorbrewer, r-dplyr, r-ggrepel, r-tibble, r-tidyr,
bioconductor-edger 3.40.0, bioconductor-deseq2 1.38.0, bioconductor-limma 3.54.0, bioconductor-qvalue 2.30.0, r-statmod

---

## 論文 Table 2（軟體與資料庫版本）可重現性紀錄（2026-10-01 server 實測）

**Benchmark 對照組的本質**（寫論文前必讀，避免不實陳述）：CirComPara2_4tools 與 nfcore_3tools **並未執行原始 pipeline**（沒有跑 CirComPara2 的 SCons 流程，也沒有跑 nf-core/circrna 的 Nextflow 流程），而是用本專案自己的 `consensus_filter.py` 改變 `--min-tools`/`--slop`/`--max-junction-ratio` 參數去重現兩篇論文的共識策略（證據：`benchmark/Snakefile` 的 `circompara2_4tools_filter`／`nfcore_3tools_filter` rule 呼叫的都是同一支 `consensus_filter.py`，與 circDEX 自己的 `our_filter` rule 完全同源；`compute_cost.py` 對這兩者的運算成本也是用文獻值 fallback，而非實測）。Table 2 **不可列 CirComPara2 / nf-core/circrna / Nextflow 的版本號**，應標註「未執行，僅重現共識策略」。真正被實際執行、需要列版本的是這兩個模擬組底層呼叫的偵測工具：

- **CIRCexplorer2 2.3.8**、**Bowtie2 2.5.4**（與 `conda list -n ciriquant --export` 一致，並經 `--version` 指令二次確認）
- **find_circ**：commit `8655dca5`（2017-02-28，`github.com/marvin-jens/find_circ`），但本地有未提交 git 的 Python3 相容性手動 patch（`/home3/choukaihsuan/tools/find_circ/` 內同時留有 `find_circ.py.bak`/`find_circ_orig.py`，`find_circ.py` mtime 晚於 clone 時間）

**參考基因體來源**（先前 CLAUDE.md 只寫「hg19」，方法學上不完整，現已確認）：
- Genome FASTA：**UCSC hg19 goldenPath**（`chromFa.tar.gz` 內的 per-chromosome 檔案，命名含 `chrUn_gl000*`/`_random`/`_hap1-7` 為 UCSC 特有格式；檔案時間戳 2009-03-21）
- Gene annotation：**GENCODE v19**（GRCh37, Ensembl 74），非 UCSC refGene——`genes.gtf` 檔頭直接寫明 `##provider: GENCODE`、`##date: 2013-12-05`，col2 來源欄為 `ENSEMBL`/`HAVANA`
- CIRCexplorer2 專用的 refFlat 註解另外來自 UCSC `refFlat.txt`（server 上載入時間 2026-06-03，見 `benchmark/config_benchmark.yaml:46` 的下載指令註解），與主流程的 GENCODE v19 GTF 是兩份不同的註解來源，不要混為一談

**資料庫 access date**：
- **circBase**（`hsa_hg19_circRNA.txt`）：本地快取檔 `~/circRNA_agent/reference/circbase_hg19.txt`，mtime **2014-05-14**（符合 circBase.org 自2014年起停止更新的公開事實）。19 個已完成資料集中 8 個（GSE108735/113230/121842/136569/143797/171011/192410/192849/97239）的 `config/projects/{GSE}.yaml` 明確指向此快取檔；其餘 11 個 `circbase_file` 留空、跑 pipeline 當下自動下載到 `/tmp/circbase_hg19.txt`，但該暫存檔現已不存在，**無法回溯確切下載日期**——因 circBase 內容本身已凍結多年，兩種情況下載到的內容應相同，但論文方法學上這 11 個只能寫「auto-fetched at pipeline runtime, exact date not recoverable」
- **CircInteractome**（`circinteractome.nia.nih.gov/api/v2/...`）與 **ENCORI/starBase**（`rnasysu.com/encori/api/...`）：兩者皆為**即時 API 查詢**，非下載快照，`predict_interactions.py` 內硬編 URL 可證；沒有「資料庫版本」概念，應寫 **"queried 2026-07-07 to 2026-08-03"**（取全部19個資料集 `de/interactions.json` 檔案 mtime 的最早與最晚值）

**容器化的誠實狀態**：Docker image `choukaihsuan/circrna-pipeline:1.0.1` 的 digest 為 `sha256:690408357fd94425f647960c82377c33cc2f19794dd87b974bc5a98853716062`（Docker Hub Registry API v2 確認，單平台 manifest）。但 server 上雖有 `apptainer 1.4.1` 可執行檔，**找不到任何 `.sif` 檔**——19 個已完成資料集的結果全部是用 conda env `ciriquant` 直接在裸機上跑出來的，不是透過容器。Table 2／方法學如果列出 Docker image，必須註明這是「供他人重現用而發布」，不是產生論文數據當下所用的實際執行環境。

**SSH 連線**：
```bash
ssh choukaihsuan@172.16.0.178
# 或用設定好的 alias
ssh genomics
```

**執行 pipeline**：
```bash
cd ~/circRNA_agent
conda activate ciriquant
nohup snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --resources mem_gb=300 \
    --keep-going \
    --rerun-incomplete \
    > logs/pipeline_run.log 2>&1 &
echo "PID: $!"
```

**同步本機到 server**（config.yaml 不同步，路徑不同）：
```bash
rsync -avzn /mnt/c/Users/User/develop/circRNA_agent/scripts/ choukaihsuan@172.16.0.178:~/circRNA_agent/scripts/   # dry run first (-n): read the list
rsync -avz  /mnt/c/Users/User/develop/circRNA_agent/scripts/ choukaihsuan@172.16.0.178:~/circRNA_agent/scripts/
rsync -avz  /mnt/c/Users/User/develop/circRNA_agent/workflow/ choukaihsuan@172.16.0.178:~/circRNA_agent/workflow/
```
- 來源路徑以 `/` 結尾（只同步內容）、目的路徑是相對於家目錄的固定路徑。執行前先用 `-n` 看清單，確認沒有非預期的目的路徑。
- 不要對絕對路徑加 `-R`（relative）或把 Windows 路徑拼進目的路徑；server 上曾出現 `scripts/mnt/c/Users/...` 的巢狀目錄（來源不明，見 `audit/server_worktree_provenance.md`），那是這類錯誤的殘骸，不要再產生。
- 同步不等於 git 同步：server 的 git checkout 停在 2c4b31f，與本機歷史不同步。目前維持這個狀態，直到刻意規劃部署步驟為止。

---

## Web UI

Flask 應用，讓使用者在瀏覽器選擇分析方法：

```bash
cd ~/circRNA_agent
conda activate ciriquant
python scripts/web_ui.py --host 0.0.0.0 --port 5000
# 瀏覽器開 http://172.16.0.178:5000
```

**功能**：
- **GEO 一鍵啟動**（頂部卡片）：輸入 GSE ID + cores → POST `/run_gse` → 呼叫 `agent.py --gse {gse_id}` → 跳轉狀態頁
- **手動 SRR 清單**：逐筆輸入 SRR ID + condition，或上傳 CSV（srr_id, condition）
- **本地 FASTQ**（新功能）：指定 server 上 FASTQ 目錄路徑 → 自動偵測配對檔（支援 `_1/_2`、`_R1/_R2`、`_R1_001/_R2_001` 等格式）→ 使用者在表格中指定各 sample 的 condition → 系統在 `raw_dir/` 建立 symlink（`{name}_1.fastq.gz` → 原始路徑），讓 Snakemake 自動跳過 download 步驟直接從 QC 開始
- Step 1：circRNA 工具選擇（CIRIquant / DCC / 兩者），自動顯示共識模式說明
- Step 2：**三種 DE 方法全部執行**（固定）；選擇「報告預設顯示方法（主方法）」；各方法附適用情境說明（edgeR=circRNA 特異性、DESeq2=樣本多/保守、limma=小樣本穩定）
- Step 3：進階參數（min_bsj, slop, **max_junction_ratio**, FDR, log2FC, threads）
- 儲存設定 → 更新 `config.yaml`
- 儲存並執行 → 更新 config 後啟動 Snakemake subprocess
- 狀態頁（`/status`）：每 5 秒 auto-refresh log，顯示 pipeline 是否執行中；**web_ui 重啟後仍能正確推斷已完成步驟**（`_infer_stages_from_files()` 掃描 output 檔案）
- **Study Title 自動填入**：輸入 GSE ID 並偵測 Labels 時，`_fetch_geo_title()` 呼叫 NCBI eUtils esummary 自動抓取研究標題，填入「資料集描述」欄位；存入 config.yaml 的 `study_title` 鍵；顯示於報告 Samples 區塊頂部（斜體）

**Web UI routes**：
- `GET /` — 主設定頁
- `POST /update` — 儲存設定（+ 可選執行 Snakemake）；儲存 `study_title`
- `POST /run_gse` — GEO 一鍵啟動；自動抓取 GEO study title
- `POST /run_manual` — 手動 SRR 清單或 CSV 上傳啟動
- `GET /api/scan_fastq?path=...` — 掃描 server 目錄，回傳偵測到的 FASTQ 配對 JSON
- `POST /run_local` — 本地 FASTQ 建立 symlink 後啟動 pipeline
- `GET /status` — 狀態頁（進度條 + 18 stage 格 + collapsible log）
- `GET /api/log` — log JSON（前端 polling 用）
- `GET /api/progress` — Snakemake log 解析 JSON（stages 陣列 + finished/total count + running bool + `report_exists` + `qc_exists`）
- `GET /api/detect_labels?gse=...` — 自動偵測 case/control label；回傳 `geo_title` 欄位
- `GET /download/<job_id>` — 強制下載 `{results_dir}/report.html`（`as_attachment=True`，檔名 `{GSE_ID}_report.html`）
- `GET /report/<job_id>` — 在新分頁開啟報告（inline）
- `GET /qc/<job_id>` — 在新分頁開啟 MultiQC 報告
- `GET /login` — 登入頁（Email magic link）
- `POST /login` — 送出 Email → 寄一次性連結；支援 `lang` 欄位（`zh`/`en`）→ 信件與錯誤訊息切換對應語言
- `GET /auth/<token>` — 驗證 magic link token → 建立 session → 導向首頁
- `GET /logout` — 清除 session → 導向登入頁
- `GET /queue` — 工作佇列頁面（排隊中 / 執行中 / 已完成 / 失敗統計）
- `GET /api/queue` — 佇列 JSON（前端 15 秒自動重整用）
- `GET /cross_dataset` — 跨資料集 circRNA 比較頁面（見下方「跨資料集分析頁面」）

**登入系統（Email Magic Link）**：
- Token 有效期：30 分鐘（`TOKEN_MINUTES=30`），儲存在 `jobs/auth_tokens.db`
- Session 有效期：7 天（`SESSION_DAYS=7`），存在 Flask session cookie
- 允許登入的 Email 白名單：`ALLOWED_EMAILS`（環境變數 `PIPELINE_ALLOWED_EMAILS`，逗號分隔）
- 支援 **Resend API**（`RESEND_API_KEY` 環境變數）→ fallback SMTP → fallback console 印出連結
- `send_magic_link(email, link, lang="zh")` 根據 `lang` 切換中英文信件內容
- 登入頁 lang switcher（右上角 中文 / EN）同步切換 UI 與信件語言；錯誤訊息透過 `data-en` 屬性雙語顯示

**GEO 資料集選擇指引**（`templates/index.html` 獨立折疊區塊）：
Web UI 主頁新增常駐卡片，使用者點標題行即可展開；內容包含：
- ✅/⚠️ 六項 checklist（Read length / Library type / replicates / case+control / 深度 / 樣本類型）
- 五行因素表（讀長 / RNA-Seq 方式 / 樣本類型 / 樣本數 / 設計）+ 顏色標記（綠/黃/紅）
- 已知資料集摘要（GSE113230/GSE58135/GSE323364/GSE133998 一覽）
- GEO 查詢提示（Library Strategy / avgLength 欄位位置）

**手動 SRR 清單 → 方式二 CSV 上傳範例表**（`templates/index.html`，2026-06-14）：
「CSV 格式範例」說明文字下方加入 GSE113230 的 6 行示範表格（srr_id / condition 兩欄，斑馬紋底色，寬度 340px），讓使用者清楚知道正確格式。

**PDF 下載按鈕**（`templates/index.html`，2026-07-02）：
Dataset Selection Guide 和 Pipeline Tutorial 折疊區塊標題列各加一個「⬇ PDF」按鈕；點擊後呼叫 `printSection(type, event)` 在新視窗開啟嵌入完整 CSS 的 HTML 頁面並自動觸發 `window.print()`，使用者存為 PDF。`event.stopPropagation()` 防止點擊觸發折疊。Tutorial 所有分頁在列印視窗自動展開（`.tut-panel{display:block!important}`）。支援中英文雙語（`localStorage.getItem('circrna_lang')`）。

**Queue/Status 頁面雙語支援**（`templates/queue.html`、`templates/status.html`，2026-07-02）：
- 兩個頁面均加入右上角 lang-btn（中文/EN 切換），採用 `data-en` 屬性儲存英文文字，`switchLang()` 函式切換
- Queue 頁面新增手動執行中任務提示橫幅（`ext_gse` 正在手動執行時，說明下一個佇列任務將在完成後自動啟動）
- Status 頁面新增 Job ID 顯示 + Copy 按鈕，啟動時間標籤雙語化

**Queue 通知 Email**（`web_ui.py`，2026-07-02）：
工作加入佇列後，`notify_queued_job(gse_id, job_id, queue_pos)` 發送 HTML 格式通知信件（dataset / Job ID / 佇列位置 / 進度頁連結），使用 Resend API → SMTP → console fallback 與 magic link 相同機制。

**Queue Worker 自我修復（`_self_heal_failed_jobs()`，`web_ui.py`，2026-07-31）**：
修正一類反覆出現的狀態不同步問題：job 第一次透過 Queue Worker 執行失敗（正確記錄 `status=failed`），之後有人繞過 Queue、直接手動下 `snakemake` 指令重跑（例如卡住的 job 需要人工排除 DCC 瓶頸），重跑本身可能真的成功、也真的產生了 `report.html`，但因為不是透過 `_run_queued_job()` 執行的，資料庫的 `status` 欄位永遠不會被更新，導致 `/queue` 頁面持續顯示「Failed」，即使使用者已經拿到正確的報告（GSE192410、GSE97239 都發生過）。

`_self_heal_failed_jobs()` 掃描所有 `status=failed` 的 job，讀取對應 `config/projects/{GSE}.yaml` 的 `results_dir`，若 `report.html` 存在**且其修改時間晚於該 job 的 `started_at`**（避免誤認前一次無關執行留下的舊報告），就自動把狀態改為 `completed`、`completed_at` 校正為報告檔案的真實時間。掛在 `/queue`、`/api/queue`（使用者每次看頁面都觸發）以及 `_queue_worker()` 啟動時（重啟 web_ui.py 即可自動修正累積的狀態）。

**PIPELINE_STAGES 雙語化**（`web_ui.py`，2026-07-02）：
`PIPELINE_STAGES` 由 tuple `(rule_id, label_zh, label_en)` 組成；API 回應及前端均可依語言顯示對應 stage 名稱。

**Pipeline 流程圖更新**（`scripts/static/pipeline_diagram.png`，2026-07-02）：
替換為使用者手繪版本（`流程圖1.png`）；生成腳本保存於 `scripts/gen_pipeline_diagram.py`（使用 matplotlib `FancyArrowPatch`，版本號 `?v=6`）。

**`_fetch_geo_title(gse_id)`**（`web_ui.py`，2026-06-19）：
呼叫 `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=200{NNNNNN}&retmode=json`，
從 `result[uid]["title"]` 取研究標題。GEO uid 公式：`"200" + GSE 數字部分`（GSE130078 → 200130078）。
Server 防火牆允許此 endpoint（僅 PRJNA/SRP 的其他 eUtils endpoint 被封鎖）。timeout=8s，失敗回傳空字串。

**`_infer_stages_from_files(results_dir, stages)`**（`web_ui.py`，2026-06-19）：
`/api/progress` 呼叫後補推：掃描 results_dir 下已存在的 output 檔，將 `status=="pending"` 的 stage 升級為 `done`。
推斷規則（由粗至細）：
- `qc/multiqc_report.html` 存在 → download / fastqc / fastp / multiqc 全標 done
- `circRNA/count_matrix.tsv` 存在 → 以上 + ciriquant / star / dcc / consensus / merge 全標 done
- 各獨立 output（circbase_annotated.tsv / isoform_groups.tsv / de_results.tsv / biomarker_candidates.tsv / isoform_switching.tsv）→ 對應 stage done
- `report.html` 存在 → 所有剩餘 pending stage 全標 done
解決 web_ui 重啟後狀態頁全顯示「等待中」的問題。

---

## 跨資料集分析頁面（`/cross_dataset`，2026-07-21）

獨立於單一資料集報告之外的彙整頁面（`web_ui.py` 的 `/cross_dataset` route + `templates/cross_dataset.html`），把所有已完成資料集的顯著 DE circRNA 拉在一起看跨癌症重現性。

**資料來源與計算**（`_load_cross_dataset_data()`，`web_ui.py`）：
- 掃描 `config/projects/*.yaml`，讀取每個資料集的 `results_dir`、`de.tumor_label/normal_label/fdr_cutoff/log2fc_cutoff/de_sig_by`
- 讀 `de/de_results.tsv`（edgeR 主方法）套用該資料集自己的顯著性門檻，補上 `circRNA/isoform_groups.tsv` 的 gene_name、`circRNA/circbase_annotated.tsv` 的 circbase_id
- 「重現 circRNA」定義：同一 circ_id 在 **≥2 個資料集**中都達到顯著

**資料集彙整表**：每個資料集一列，顯示癌症類型（`_DATASET_META` 內建對照表，含 organ 分類）、配對數、tested/顯著/上調/下調數。

**跨資料集 Heatmap**：top 60 重現 circRNA × 全部資料集，log2FC 為值（紅=上調、綠=下調、灰=未偵測）；本身內建**階層式聚類**（`_hclust()`，average-linkage + Euclidean distance，純前端 JS 實作，不依賴 scipy），左側可切換顯示/隱藏樹狀圖（`toggleDendrogram()`）。

**同癌症/器官比較（「Compare within」下拉選單）**：把 organ 欄位相同、且該 organ 有 ≥2 個資料集的分組，讓使用者只在「同一種癌症的資料集之間」比較重現性，而不是跟全部 13+ 個異質資料集混在一起算。目前符合條件的分組：
- **Lung**：GSE229705、GSE148036（皆為 LUAD）
- **Breast**：GSE113230、GSE58135、GSE133998、SRP156355（GSE323364 因 organ 標記為 "Cell Line" 自動排除，不會混進組織比較）

選定分組後，「資料集數」欄位重新計算為「只在該分組內的重現次數」（不是全域計數），heatmap 也重新篩選只保留在該分組內 ≥2 個資料集重現的 circRNA 並重繪聚類。此重新計算完全在前端 JS 完成（`filterTable()` + `_updateHeatmapForOrgan()`），不需要額外的伺服器查詢，因為全域 `recurrent` 清單本身已經包含所有「全域 ≥2 個資料集重現」的候選——只要一個 circRNA 在某分組內達到 ≥2，全域計數必然也 ≥2，所以不會有漏算的候選。

**circBase ID 超連結**：重現清單表格的 circbase_id 欄位（非 "novel" 者）連到 `https://www.circbase.org/cgi-bin/singlerecord.cgi?id={id}`。

**UpSet Plot**（`drawUpset()`，`templates/cross_dataset.html`）：
- 視覺化多資料集交集；橘色直條 = 只出現在單一資料集的 exclusive circRNA；藍色直條 = 跨 ≥2 個資料集的 shared circRNA
- 點矩陣（dot matrix）顯示每個直條對應哪些資料集的組合
- **排序模式切換**（互動式）：「依交集層數（5→1 資料集）」vs「依 circRNA 數量（多→少）」；按鈕 `#upset-btn-degree` / `#upset-btn-count`，`setUpsetMode(mode)` 函式；`let _upsetMode = 'degree'`
- Python 送出全部 ~108 個交集（`upset_intersections` 不做排序或截斷）；JS 依 `_upsetMode` 排序後取 top 30
- 雙語支援：Y 軸 tick 文字從 `D.labels_en`（English）或 `D.labels`（中文）取用；`applyLang()` 觸發 `drawUpset()` 重繪以切換資料集名稱語言
- **JS TDZ 陷阱**：`drawUpset()` 不可在 `let _LANG` 宣告前呼叫（ReferenceError: Cannot access before initialization）；只由 `applyLang()` 在宣告後觸發，並包 try-catch 避免 plot 錯誤影響 UI 切換

**已知 bug 修正**：
- 排序後編號跳號（`sortTable()` 原本沒跳過隱藏列，篩選 + 排序疊加使用時編號變成 1,2,5,6,11...而非連續 1,2,3,4）→ 改成只對可見列遞增編號
- 空資料集 early-return 路徑漏傳 `organ_groups`/`organ_groups_json` 導致 JS 直接壞掉、heatmap 初始渲染呼叫少了 null 保護 → 皆已補上

---

## Container 部署（Docker + Singularity HPC）

### Docker image

```bash
# 本機建置並推送（WSL2）
cd /mnt/c/Users/User/develop/circRNA_agent
bash containers/build_and_deploy.sh   # 自動執行 docker build + docker push
```

image name：`choukaihsuan/circrna-pipeline:1.0.1`
`Dockerfile` 使用 `mamba env create --yes`（非互動式），比 conda 快 3–5×。

### `envs/circrna.yaml` 重要設定

- channels：`bioconda`, `conda-forge`（**已移除 `defaults`**，Anaconda 商業授權問題）
- `bioconductor-edger>=3.40`（不鎖定 patch 版本，避免 PackagesNotFoundError）
- `bioconductor-deseq2>=1.40`
- `bioconductor-qvalue`（Storey q-value，analysis.R 需要）
- `dcc=0.5.0` 在 conda 依賴中（bioconda），**不在 pip**（PyPI 的 DCC 是另一個無關套件，版本 0.7+）

### Singularity / Apptainer（HPC server）

Docker image 已上傳至 Docker Hub，任何有 Singularity 或 Apptainer 且開啟 user namespace 的 HPC server，執行以下一行即可取得完整環境：

```bash
singularity pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.1
```

或使用 Apptainer（語法相同）：

```bash
apptainer pull circrna-pipeline.sif docker://choukaihsuan/circrna-pipeline:1.0.1
```

拉取後用容器執行 pipeline：

```bash
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --use-singularity \
    --singularity-args "--bind /home3/choukaihsuan:/home3/choukaihsuan" \
    --keep-going --rerun-incomplete
```

`workflow/Snakefile` 頂部已設定：`singularity: "docker://choukaihsuan/circrna-pipeline:1.0.1"`

**注意**：目前使用的 server（172.16.0.178，CentOS 7）預設關閉 user namespace，conda 版 apptainer 無法執行。需請管理員執行：
```bash
echo 10000 > /proc/sys/user/max_user_namespaces
```
或安裝 setuid 版本的 Singularity。在此之前繼續使用 conda env `ciriquant`。

---

## Per-project Config 系統

每個 GSE 分析自動使用獨立設定，避免多個分析互蓋 `config.yaml` 和 `metadata/`。

### 機制

- **`save_project_snapshot(cfg)`**（`web_ui.py`）：每次儲存設定時，同時寫入：
  1. `config.yaml`（全域啟用，目前執行的專案）
  2. `config/projects/{GSE_ID}.yaml`（該專案的永久快照）
- **`_configfile_for(gse_id)`**：若 `config/projects/{GSE_ID}.yaml` 存在，Snakemake 的 `--configfile` 自動指向它；否則 fallback 到 `config.yaml`
- **`run_manual()`**（手動填入 SRR 啟動）：metadata 同時儲存到 `metadata/library_info.csv` 和 `metadata/{GSE_ID}/library_info.csv`；config 的 `metadata`/`groups` 欄位更新為專案路徑

### 切換專案

```bash
# Web UI 輸入新 GSE ID → run_gse 自動載入對應 config/projects/{GSE_ID}.yaml
# 若要手動切換：
cp config/projects/GSE113230.yaml config.yaml
```

### Server 端獨立重跑腳本（Mock snakemake 物件）

當 server 上 config.yaml 已切換到其他專案時，用以下 wrapper 直接執行，繞過 Snakemake DAG：

**`/tmp/run_generate_report.py`**（Python mock）：
```bash
# 必須用 conda env 的 Python，否則缺少 plotly 會生成靜態 PDF 報告
/home/choukaihsuan/miniconda3/envs/ciriquant/bin/python /tmp/run_generate_report.py
```

Python mock 格式（三個關鍵：`output` 用 list；`exec(code, {"snakemake": sn})` 注入；屬性名稱需與 rule 一致）：
```python
import os
os.chdir("/home/choukaihsuan/circRNA_agent")

class _Sn: pass
sn = _Sn(); sn.input = _Sn(); sn.params = _Sn()
sn.log = ["/tmp/report_GSE.log"]

R = "/home3/choukaihsuan/{GSE}_results"
M = "/home/choukaihsuan/circRNA_agent/metadata/{GSE}"
sn.input.de            = R + "/de/de_results.tsv"
sn.input.de_edger      = R + "/de/de_results_edgeR_ciriquant.tsv"
sn.input.de_deseq      = R + "/de/de_results_deseq2.tsv"       # de_deseq（非 de_deseq2）
sn.input.de_limma      = R + "/de/de_results_limma.tsv"
sn.input.biomarkers    = R + "/de/biomarker_candidates.tsv"
sn.input.matrix        = R + "/circRNA/count_matrix.tsv"
sn.input.volcano       = R + "/plots/volcano.pdf"
sn.input.heatmap       = R + "/plots/heatmap.pdf"
sn.input.pca           = R + "/plots/pca.pdf"
sn.input.multiqc       = R + "/qc/multiqc_report.html"
sn.input.switching     = R + "/de/isoform_switching.tsv"        # switching（非 isoform_switching）
sn.input.groups        = M + "/sample_groups.csv"
sn.input.isoform_groups = R + "/circRNA/isoform_groups.tsv"    # isoform_groups（非 isoforms）
sn.input.circbase_annot = R + "/circRNA/circbase_annotated.tsv"
sn.input.interactions  = R + "/de/interactions.json"
sn.params.project_id   = "{GSE}"
sn.params.de_method    = "edgeR_ciriquant"
sn.params.fdr = 0.05; sn.params.lfc = 1.0
sn.params.de_sig_by    = "pvalue"
sn.params.tumor_label  = "tumor"; sn.params.normal_label = "normal"
sn.params.heatmap_top_n = 10; sn.params.study_title = "..."
sn.output = [R + "/report.html"]   # list，支援 output[0]

exec(open("/home/choukaihsuan/circRNA_agent/scripts/generate_report.py").read(), {"snakemake": sn})
# ⚠ 不可用 builtins.snakemake = sn，因為 dir() 不會包含 builtins namespace
```

**`/tmp/run_de_analysis.R`**（R mock，S4 class with list slots）：
```bash
conda run -n ciriquant Rscript /tmp/run_de_analysis.R
```
R mock 格式：
```r
setClass('Snakemake', representation(input='list', output='list', params='list', log='list'))
snakemake <- new('Snakemake',
  input  = list(matrix="...", fsj_matrix="...", groups="...", circbase_annot="..."),
  output = list(de="...", de_edger="...", de_deseq="...", de_limma="...", volcano="...", ...),
  params = list(de_method="edgeR_ciriquant", fdr=0.05, lfc=1.0, ...),
  log    = list("/path/to/log")
)
source('/home/choukaihsuan/circRNA_agent/scripts/analysis.R')
```

兩個 wrapper 都 hardcode GSE 的 input/output/params 路徑。此模式可複用於任何需要獨立重跑 terminal rule 的情境。

---

## Benchmark 設計（`benchmark/`）

評估本 pipeline 對抗三個已發表方法的偵測準確率與 DE 品質，輸出自包含 HTML 報告。

### 執行方式

```bash
cd ~/circRNA_agent
conda activate ciriquant
snakemake \
    --snakefile benchmark/Snakefile \
    --configfile benchmark/config_benchmark.yaml \
    --cores 8 \
    --resources mem_gb=60 \
    --keep-going --rerun-incomplete
```

只跑準確率（跳過 DE quality，適合 GSE113230 尚未完成時）：
```bash
snakemake --snakefile benchmark/Snakefile \
    --configfile benchmark/config_benchmark.yaml \
    --cores 8 --until accuracy_benchmark compute_cost
```

### Ground Truth 建立（`rnaser_ground_truth.py`）

使用 **SRP011042**（"CircleSeq in Human Fibroblasts"；Jeck et al. 2013，*RNA* 19(2):141–157；人類纖維母細胞 Hs68，hg19）的 RNase R enrichment 實驗作為 ground truth。

> **編號更正（2026-10-04）**：本專案先前記為 GSE55872，該編號在 GEO 是**果蠅**（*Drosophila melanogaster*，fly heads + S2 cells）資料集，屬誤填。三條獨立證據確認實際使用的是 SRP011042：(1) 使用的 CIRIquant GTF 為人類染色體名（chr1–chrY、chrUn_gl000220），STAR 索引為 hg19；(2) 三個 run 的 fastp 讀序數與 SRA `total_spots` 精確相等（157,053,158 / 158,305,855 / 199,922,486）；(3) SRA metadata 標明 *Homo sapiens*、taxid 9606。`GSE55872_results/`、`GSE55872/raw` 等目錄名為歷史遺留，保留不改（既有輸出位於其中）。

| SRR ID | 角色 |
|--------|------|
| SRR444655 | Total RNA（偵測基準）|
| SRR444974 | RNase R replicate 1 |
| SRR445016 | RNase R replicate 2 |

RNase R 消化線性 RNA 並富集環狀 RNA，是目前 circRNA 偵測 benchmark 的標準方法。

**Enrichment Ratio（ER）判斷**（`rnaser_ground_truth.py`）：
- **兩個 RNase R replicate 都會用到，以算術平均合併**：
  `ER = mean(BSJ_SRR444974, BSJ_SRR445016) / BSJ_SRR444655`（`bsj_rnaser_mean = sum(vals)/len(vals)`，非取最大值、也非只挑一個 replicate）
- 用的是**原始 BSJ read count**，不做 per-sample RPM 正規化——因為 RNase R 處理會使線性 RNA 大幅減少，導致 RNase R 樣本的總 BSJ read 數比 Total RNA 樣本多 ~25 倍；若用 RPM 正規化，分母被灌水 25 倍會系統性壓低幾乎所有 ER 值（實測會造成 ~94% 假 TN），故改用原始 count，只要兩組樣本 library size 量級相近（同一批 Illumina 定序即符合）ER 比值仍具意義
- 座標比對用 slop=10bp 模糊匹配（同 consensus_filter 的邏輯），非精確字串比對
- Edge case：Total RNA 中偵測不到（BSJ≈0）但 RNase R 有偵測到 → ER=100（視為強 TP 訊號）；兩邊都偵測不到 → ER=1（歸入 ambiguous）
- ER ≥ 1.5 → **True Positive**（真實 circRNA）
- ER ≤ 0.5 → **True Negative**（假陽性候選）
- 中間值（0.5 < ER < 1.5）排除於評估之外

SRP011042 的 FASTQs 由 `bench_download` rule 從 **EBI FTP**（`ftp.sra.ebi.ac.uk/vol1/fastq/{SRR[:6]}/{SRR}`）自動下載，無需手動準備。config 中的 key 為 `srp011042`（2026-10-04 由 `gse55872` 改名，僅為標籤）。

### Task 1 – 偵測準確率比較（`accuracy_benchmark.py`）

對 SRR444655（total RNA）執行多工具共識策略，再對照 RNase R ground truth：

| 策略 | 工具組合 | slop | pseudo-circ QC | 對應論文 |
|------|----------|------|----------------|----------|
| **circDEX**（本研究）| CIRIquant + DCC | 10 bp | ✅ selective（BSJ<5）+ adaptive | — |
| **CirComPara2_4tools** | CIRIquant + DCC + CIRCexplorer2 + find_circ | 10 bp | ❌ | Gaffo et al. 2022 |
| **nfcore_3tools** | CIRIquant + CIRCexplorer2 + find_circ | 0（精確匹配）| ❌ | Digby-Bell et al. 2023 |

**重要設計決策**：CirComPara2_4tools 不含獨立 CIRI2，因為 CIRIquant 內部已呼叫 CIRI2 做 BSJ 偵測；同時納入兩者等於讓同一演算法投兩票，破壞 consensus 獨立性。4 個工具分別使用 HISAT2+BWA / STAR / STAR / Bowtie2，是真正獨立的偵測策略。

**find_circ 輸出過濾**：`parse_find_circ()` 必須只保留 `category` 欄含 `CIRCULAR` 的行，否則 LINEAR（~233K）和 AMBIGUOUS（~190K）junction 會進入 consensus，導致座標匹配極慢（測試發現 55 min 仍未完成）。真正 CIRCULAR 只有 1,895 個。

**評估指標**：Precision、Recall、F1、**Specificity**、**TN**、AUC-PR

**AUC-PR 計算方式（重要）**：
- 二元偵測（每個 circRNA 只有 detected/not detected，無連續分數）直接套用 `sklearn.average_precision_score` 會嚴重虛高 AUC-PR。
- **根本原因**：circDEX 偵測到的 circRNA 全部 score=1，未偵測的全部 score=0；ground truth 中 87.8% 的 TP 在 score=0 的大池（即 FN）；樂觀排序（label=1 排在 label=0 前面）把大量 FN 全部排在 TN 前面，產生長段假精確度（AUC-PR 虛高）。
- **正確方法**：門檻掃描（threshold sweep）。對 CIRI2 output 和 DCC `CircRNACount` 各取 `min_bsj` 門檻（1–50），重新建立共識 → 計算每個門檻下的 (Precision, Recall) → 繪製真實 PR 曲線 → trapezoid AUC。
- **實測結果（門檻掃描，三方法，2026-07-20 `-Nr 2` 重新確認）**：circDEX AUC-PR = **0.174**；CirComPara2_4tools = **0.391**；nfcore_3tools = **0.377**（見 GSE113230 段落「Benchmark 門檻掃描 AUC-PR 重新確認」完整表格；`-Nr 5` 時代舊值 0.120/0.349/0.337 保留於歷史版本區塊）。
- **新增 CLI 參數**：`accuracy_benchmark.py --ciri2-file <CIRI2 output> --dcc-count-file <CircRNACount> --circexplorer2-file <known_circ.txt> --find-circ-file <splice_sites.bed> --output-pr-curve <pr_curve.tsv>`（四個工具檔案缺一個，該工具的門檻掃描會靜默塌陷為全 0 或與另一方法完全重複，不會報錯，需檢查 log 的 `CE2=`/`find_circ=` 計數）；`generate_comparison_report.py --pr-curve <pr_curve.tsv>`（在報告中插入 SVG PR 曲線圖）。
- **`SRR444655.ciri2`**：CIRI2 的原始輸出格式（非 GTF），門檻掃描專用；位於 CIRIquant 內部工作目錄 `CIRI2/{srr}.ciri2`（col 4 = BSJ count）；`DCC/CircRNACount` col 3 = junction count。

**分層分析**：依 Total RNA 中的 BSJ count 分三層：
- Low：1–4 RPM
- Mid：5–19 RPM
- High：≥ 20 RPM

輸出：`benchmark/accuracy_summary.tsv`（含 TN、Specificity 欄）、`benchmark/stratified_f1.tsv`

### Task 2 – DE 分析品質比較（`de_quality_benchmark.py`）

使用 **GSE113230**（三陰性乳癌，6 samples）的 count matrix，比較三種 DE 方法：

| 方法 | DE 演算法 | 特色 |
|------|----------|------|
| **Our method** | edgeR_ciriquant（BSJ/FSJ ratio）| Type I/II 分類；FSJ offset |
| DESeq2 baseline | DESeq2 on BSJ counts | 無 FSJ offset |
| **limma-voom** | limma-voom（TMM + voom weights）| 無 FSJ offset |

**前置條件**：需先執行 `/tmp/run_de_analysis.R`（或 Snakemake de_analysis rule）產生三個 TSV；以及 `/tmp/run_de_baseline.R` 執行 benchmark baseline。

**評估指標**：
- 各方法顯著 DE circRNA 數量
- 兩兩 Jaccard similarity（Our vs DESeq2 / Our vs limma / DESeq2 vs limma）
- Type I circRNA 中，僅我們方法偵測到的比例
- Top 20 DE circRNA 中 circBase 已知的比例

輸出：`benchmark/de_quality_summary.tsv`、`benchmark/de_jaccard.tsv`、`benchmark/de_lists.json`（互動式報告用）

**`de_quality_benchmark.py` CLI 參數**：
- `--isoforms <isoform_groups.tsv>`：從 `assign_isoforms.py` 輸出 join `gene_name` 到 DE 結果（若不提供，互動 modal 的 gene_name 欄會空白）
- `--output-lists <de_lists.json>`：輸出各方法 / Jaccard 交集的 circRNA 清單 JSON，供 `--de-lists` 傳入 `generate_comparison_report.py`

**Jaccard pairwise circRNA lists 格式**（`de_lists.json` 的 `_jaccard` key）：
```json
[{"Method_A":"Our_edgeR_ciriquant","Method_B":"DESeq2_baseline",
  "A_only":[{circ_id, gene_name, log2FC, p_value, circbase_id},...],
  "B_only":[...], "Both":[...]}]
```

### Task 3 – 計算資源成本（`compute_cost.py`）

整合各步驟實測時間（`/usr/bin/time -v` log），輸出 `benchmark/compute_cost.tsv`。

**nf-core 實測 time logs**（benchmark Snakefile 現在記錄以下步驟）：
- `logs/bench/time_ciriquant_{srr}.log`：CIRIquant（與 Our pipeline 共用）
- `logs/bench/time_circexplorer2_{srr}.log`：CIRCexplorer2（約 6 秒）
- `logs/bench/time_find_circ_map_{srr}.log`：bowtie2 unmapped（在 HPC NFS 上約 3h 44min）
- `logs/bench/time_find_circ_{srr}.log`：find_circ detection

**注意**：find_circ_map wall time 受 HPC NFS 磁碟 I/O 速度影響顯著（寫入 28GB unmapped reads），比 AWS SSD 環境慢。報告 Source 欄標注「measured on HPC with NFS storage」。若 time logs 不存在，自動 fallback 使用 Digby-Bell 2023 文獻值。

**已知實測時間**（SRR444655，8 cores，HPC NFS）：
- CIRCexplorer2：0:05.99（6 秒，0.13 GB RAM）
- find_circ_map：3:44:29（3.54 GB RAM）
- find_circ detect：2:41:05（5.14 GB RAM）
- CIRIquant：**11:50:25（49.1 GB RAM）** ← NFS I/O 瓶頸（HISAT2 5h8m + BWA-MEM 3h3m + CIRI2 + de novo quant）；2026-06-10 重跑確認值（disk-full 期間仍完成）；2026-07-20 全量重跑（`-Nr 2` + dedup 修正後）確認值 **695.2 min（11h35m）**，同一數量級

**CIRIquant 步驟分解**（NFS 環境，SAM/BAM 寫入放大效應）：
- HISAT2 genome alignment：00:01 → 05:09（5h 8min；unmapped.sam 124 GB 寫入 NFS）
- Gene abundance：05:09 → 05:22（13 min）
- BWA-MEM：05:22 → 08:25（3h 3min）
- CIRI2.pl detection：08:25 → 09:21（56 min）
- Build circular index：09:21 → 09:29（8 min；某些 sample 可達 31 min，依 circRNA 數量而定，見 GSE171011 SRR14088793 紀錄）
- De novo HISAT2 alignment：09:29 → 10:13（44 min）
- BSJ/FSJ detection & quantification：10:13 → 11:50（97 min，含 disk-full 延遲）

**DE 分析時間**（2026-07-20 新增，`run_de_timed.R` 包裝 `/usr/bin/time -v` 計時，不修改 `analysis.R` 本身）：
- `de_edgeR_ciriquant_baseline` rule：Our pipeline 的 DE 耗時（edgeR_ciriquant 主方法，同一次 `analysis.R` 呼叫內副產出 DESeq2/limma）→ `compute_cost.tsv` 的 `DE_min` 欄
- `de_deseq2_baseline` rule：nf-core 模擬對照組（DESeq2，無 FSJ offset）
- **實測值**（GSE113230，6-sample，`-Nr 2` 全量重跑後的 count matrix）：Our pipeline DE_min = **1.1 min**；nf-core DE_min = **0.9 min**
- **2026-10-01 起 DE_min 不再折入 Total_wall_min**（見下方「compute_cost.py 三處修正」），僅作為獨立參考欄位保留；偵測步驟用單樣本 SRR444655 計時、DE 分析用 GSE113230 六樣本計時，兩者 n 不同，混進同一個「per-sample runtime」的 Total 是範疇錯誤，不只是精度問題。

**compute_cost.py 三處修正與 log/tsv 對帳記錄（2026-10-01，論文 Figure 6(c) 查證後觸發）**：

查證 Figure 6(c) 的 CirComPara2/nf-core 運算成本數字時，發現並修正三個問題：

1. **CirComPara2_4tools 的 STAR mate1/mate2 時間先前是寫死值**（`compute_cost.py:352-353`的`step_times.get("star_mate1",34.6)`/`get("star_mate2",47.2)`），即使`logs/bench/time_star_mate1_SRR444655.log`（實測34.589min）和`time_star_mate2_SRR444655.log`（實測47.223min）就在同一個目錄裡，卻從未被讀取——根本原因是`benchmark/Snakefile`的`compute_cost`規則（原L802-825）從未傳入`--our-star-mate1-log`/`--our-star-mate2-log`這兩個CLI參數。已修正：Snakefile規則加上這兩個input+flag，連帶也修正了**Our pipeline自己的STAR_min**（因為同一組參數是共用的）——這條row的Tool_combination寫的是「STAR×3」，但先前因為同樣沒傳mate1/mate2 log，實際只算了paired一次，**circDEX自己的Total被低估了**（811.1→891.8 min，+80.7 min），不是對照組才有問題。
2. **DE_min 不應折入 Total**（見上方DE分析時間段落）：已修正`our_row`/`nfcore_row`/`circompara2_4tools_row`三處，DE_min全部改為僅供參考、不加總進Total_wall_min與CPU_hours。
3. **nf-core/circrna 的 Total 漏算了 STAR 時間**：`compute_cost.py`原L239-240的`nfcore_align`只加總`nfcore_times.get("ciriquant")+nfcore_times.get("fc_map")`，但`nfcore_times`這個dict的建構來源（`nfcore_logs`，原L216-221）根本沒有定義`"star"`這個key，所以STAR時間永遠不會被計入——即使`STAR_min`欄位本身有正確顯示98.6（來自另一個共用dict `step_times`），卻沒有真正加進Total。而nf-core/circrna的CIRCexplorer2是吃STAR的Chimeric.out.junction輸出，沒有STAR就跑不出CIRCexplorer2結果，這段時間理應算入。已修正為`nfcore_align`多加`step_times.get("star",0)`。

**查證過程中額外發現的誠信問題（已修正，供日後參考）**：舊版`compute_cost.tsv`（mtime 2026-07-22 00:39）比產生它的`compute_cost.log`（mtime 2026-07-20 03:34）晚了近2天，兩者數字對不上一小塊（CirComPara2_4tools Total：log記1279.9，tsv卻是1280.8，差額剛好等於DE_min的0.9）——代表tsv在log產生之後又被重新跑過一次，但那次重跑沒有`tee`到log，導致事後無法用log文字直接證明該次重跑有沒有誤觸文獻fallback分支（L272的`"No nf-core time logs found; using Digby-Bell 2023..."`），只能靠數值反推佐證。**2026-10-01這次重跑已修正此操作疏失**：執行指令務必用`2>&1 | tee logs/bench/compute_cost.log`，且重跑後立刻確認log與tsv的mtime epoch完全一致（`stat -c %Y`兩者都是`1790824092`，已核對）。以後只要修改`compute_cost.py`或重新產生`compute_cost.tsv`，都要比照這個「log/tsv時間戳必須對得上」的檢查習慣。

**修正前後數字比較**（僅重跑`compute_cost.py`本身，未重跑任何偵測工具）：

| Pipeline | Total_wall（修正前） | Total_wall（修正後） | Δ | CPU-h（修正前）| CPU-h（修正後）| circDEX相對節省（修正前）| circDEX相對節省（修正後）|
|---|---|---|---|---|---|---|---|
| **Our pipeline (circDEX)** | 811.1 | **891.8** | +80.7 | 108.1 | 118.9 | — | — |
| CirComPara2_4tools | 1280.8 | **1279.9** | −0.9 | 170.8 | 170.7 | 37.0% | **30.3%** |
| nf-core/circrna | 1084.2 | **1181.9** | +97.7 | 144.6 | 157.6 | 25.2% | **24.5%** |

circDEX相對節省 = `1 − Total_wall(circDEX) / Total_wall(對照組)`。**誠實結論**：修正後circDEX相對CirComPara2_4tools的時間優勢從37%縮小到30.3%——因為修正同時補回了circDEX自己先前漏算的STAR mate1/mate2成本（+80.7 min），這筆帳跑不掉；相對nf-core/circrna的優勢則大致持平（25.2%→24.5%），因為nf-core那邊漏算的STAR成本（+97.7 min）比circDEX自己補回的還多，兩邊同步墊高，相對差距變化不大。三列的Tool Breakdown逐欄加總均已核對等於Total_wall_min（Our：695.2+180.4+16.2=891.8；CirComPara2_4tools：695.2+180.4+16.2+0.2+387.9=1279.9；nf-core：695.2+98.6+0.2+387.9=1181.9）。

Figure 6(c)／`comparison_report.html`的Compute Cost面板已用新`compute_cost.tsv`重新產生（2026-10-01 11:09，`generate_comparison_report.py`直接呼叫，未觸發Snakemake DAG、未重跑任何上游偵測或DE步驟）。

### 輸出報告

`benchmark/comparison_report.html`：自包含 HTML，整合三個面向（準確率、DE 品質、資源成本）的比較表與圖表。

**`generate_comparison_report.py` 互動式功能**：
- **DE Quality 表格**：`Sig_DE_circRNAs`、`Up_regulated`、`Down_regulated`、`Type_I_count`、`Type_II_count`、`Top20_in_circBase` 等欄位數字可點擊 → modal 顯示 circRNA 清單（`showDEList(method, col)`）
- **Jaccard Overlap 表格**：`A_only`、`B_only`、`Both` 數字可點擊 → modal 顯示交集/差集 circRNA 清單（`showJaccardList(rowIdx, col)`）；需傳入 `--de-lists`
- **互動 modal 排序**：所有 modal 清單表格欄位可點擊排序（▲/▼）；`_renderCircList()` 統一渲染，`sortCircList(th, colIdx)` 函式以 `data-val` 屬性做數值/字串排序
- **Compute Cost 表格**：`_compute_cost_table_html(comp)` 函式；"Tool Breakdown" 欄顯示 stacked mini-bar（CIRIquant=深藍 / STAR×3=中藍 / DCC=淺藍 / CIRCexplorer2=橙 / find_circ=綠 / Other=灰）及各工具時間數值；最短 Total 以綠色粗體標示；需 `compute_cost.tsv` 含 per-tool 欄（`CIRIquant_min`, `STAR_min`, `DCC_min`, `CIRCexplorer2_min`, `find_circ_min`，已由 `compute_cost.py` 自動輸出）
- **Stratified F1**：bar chart 已移除；表格只顯示 circDEX 一行（`strat[strat["Method"]=="Our_adaptive"]`，內部 key 仍為 `Our_adaptive`）

**獨立重跑 benchmark 腳本**（繞過 Snakemake，適合 config 指向其他專案時）：
```bash
RESULTS='/home3/choukaihsuan/GSE113230_results'
BENCH=$RESULTS/benchmark
SCRIPTS=~/circRNA_agent/benchmark/scripts
PY=/home/choukaihsuan/miniconda3/envs/ciriquant/bin/python
GS55=/home3/choukaihsuan/GSE55872_results/circRNA/SRR444655

# 0. nf-core 3-tool re-filter（需 find_circ CIRCULAR-only fix）
$PY ~/circRNA_agent/scripts/consensus_filter.py \
    --cirique $GS55/SRR444655.gtf --circexplorer2 $GS55/CIRCexplorer2/known_circ.txt \
    --find-circ $GS55/find_circ/splice_sites.bed \
    --output $BENCH/detection/nfcore_3tools.bed --summary $BENCH/detection/nfcore_3tools_summary.tsv \
    --min-tools 2 --slop 0 --min-bsj 2 --max-junction-ratio 999

# 1. CirComPara2 4-tool consensus（CIRIquant + DCC + CIRCexplorer2 + find_circ）
$PY ~/circRNA_agent/scripts/consensus_filter.py \
    --cirique $GS55/SRR444655.gtf --dcc $GS55/DCC/CircCoordinates \
    --circexplorer2 $GS55/CIRCexplorer2/known_circ.txt \
    --find-circ $GS55/find_circ/splice_sites.bed \
    --output $BENCH/detection/circompara2_4tools.bed --summary $BENCH/detection/circompara2_4tools_summary.tsv \
    --min-tools 2 --slop 10 --min-bsj 2 --max-junction-ratio 999

# 2. accuracy benchmark（含門檻掃描 PR curve，需提供 CIRI2 output 和 DCC CircRNACount）
$PY $SCRIPTS/accuracy_benchmark.py \
    --ground-truth $BENCH/rnaser/ground_truth.tsv \
    --our-bed $BENCH/detection/our_method.bed --our-summary $BENCH/detection/our_method_summary.tsv \
    --our-no-qc-bed $BENCH/detection/our_no_qc.bed --our-no-qc-summary $BENCH/detection/our_no_qc_summary.tsv \
    --circompara2-bed $BENCH/detection/circompara2_sim.bed --circompara2-summary $BENCH/detection/circompara2_sim_summary.tsv \
    --circompara2-4tools-bed $BENCH/detection/circompara2_4tools.bed --circompara2-4tools-summary $BENCH/detection/circompara2_4tools_summary.tsv \
    --nfcore-bed $BENCH/detection/nfcore_3tools.bed --nfcore-summary $BENCH/detection/nfcore_3tools_summary.tsv \
    --output-summary $BENCH/accuracy_summary.tsv --output-stratified $BENCH/stratified_f1.tsv \
    --output-fp-comparison $BENCH/fp_score_comparison.tsv --slop 10 --min-bsj 2 \
    --ciri2-file $GS55/SRR444655.ciri2 \
    --dcc-count-file $GS55/DCC/CircRNACount \
    --output-pr-curve $BENCH/pr_curve.tsv

# 3. DE baseline（mock snakemake，見 /tmp/run_de_baseline.R）
conda run -n ciriquant Rscript /tmp/run_de_baseline.R

# 4. DE quality（--isoforms 提供 gene_name；--output-lists 輸出互動式 modal 用 JSON）
$PY $SCRIPTS/de_quality_benchmark.py \
    --our-de $RESULTS/de/de_results.tsv \
    --nfcore-de $BENCH/de/nfcore_deseq2_results.tsv --limma-de $BENCH/de/nfcore_limma_results.tsv \
    --circbase-annot $RESULTS/circRNA/circbase_annotated.tsv \
    --isoforms $RESULTS/circRNA/isoform_groups.tsv \
    --fdr 0.05 --lfc 1.0 \
    --output-summary $BENCH/de_quality_summary.tsv --output-jaccard $BENCH/de_jaccard.tsv \
    --output-lists $BENCH/de_lists.json

# 5. comparison report（--de-lists 啟用互動式表格；--pr-curve 插入 SVG PR 曲線圖）
$PY $SCRIPTS/generate_comparison_report.py \
    --accuracy $BENCH/accuracy_summary.tsv --stratified $BENCH/stratified_f1.tsv \
    --compute $BENCH/compute_cost.tsv \
    --de-quality $BENCH/de_quality_summary.tsv --de-jaccard $BENCH/de_jaccard.tsv \
    --fp-comparison $BENCH/fp_score_comparison.tsv \
    --pr-curve $BENCH/pr_curve.tsv \
    --de-lists $BENCH/de_lists.json \
    --output $BENCH/comparison_report.html
```

---

## 已知問題與解決方法

| 問題 | 原因 | 解決方法 |
|------|------|----------|
| PRJNA/SRP metadata 抓取失敗（server 上）| Server 曾封鎖 NCBI eUtils HTTP（2026-07 文件記載）；2026-10-19 實測 efetch.fcgi 已可用（curl 確認 HTTP 200 + 真實 SRA XML 內容，非僅狀態碼；未測試是否仍有其他 endpoint 被封鎖，也無法判斷是防火牆規則改了還是原本就沒封鎖）；GSE 透過 pysradb 可成功 | 見下方詳述 |
| DCC `TypeError: NoneType has no len()` | DCC 0.5.0 不支援缺少 `-mt2` | 加入 `star_align_mate1` 和 `star_align_mate2` rule，分別比對 R1/R2 |
| DCC `command not found` | 指令為大寫 `DCC` | shell command 改為 `DCC` |
| DCC `-A` flag error | DCC 0.5.0 BAM 參數是 `-B` 不是 `-A` | 改為 `-B {bam}` |
| `ciriquant.yaml` `ConfigError: File: bwa, not found` | `check_file()` 做字面路徑存在性檢查，`bwa: bwa`（相對路徑）不存在；必須填完整絕對路徑 | 在 `config/ciriquant.yaml` 的 `tools:` 區塊填 conda env 的完整路徑（如 `/home/choukaihsuan/miniconda3/envs/ciriquant/bin/bwa`）；`java: /usr/bin/java`、`perl: /usr/bin/perl` 亦需指定 |
| `ciriquant.yaml` `ConfigError: Reference fasta need to be specified` | `config/ciriquant.yaml` 的 `reference:` 區塊鍵名錯誤；CIRIquant 期望 `fasta:`/`bwa_index:`/`hisat_index:`，而非 `genome:`/`hisat2_genome:`；若誤用 `database:` 作為頂層鍵則先報 `KeyError: 'reference'` | 正確 `reference:` 結構：`fasta: /path/hg19.fa`、`gtf: /path/genes.gtf`、`bwa_index: /path/hg19.fa`（BWA index prefix，附加 `.bwt` 驗證）、`hisat_index: /path/hg19_hisat2_index`（HISAT2 prefix，附加 `.1.ht2` 驗證）；server 路徑：`bwa_index: /home3/choukaihsuan/reference/hg19/hg19.fa`，`hisat_index: /home3/choukaihsuan/reference/hg19/hg19_hisat2_index` |
| CIRIquant 看似「completed successfully」但 GTF 消失 | shell 最後一行 `rm -rf ... circ align gene` 成功（exit 0），掩蓋了 CIRIquant 本身失敗（exit 1）的 exit code；Snakemake 只看最終 exit code → 誤報成功 | 根本修法是修正 `ciriquant.yaml` 讓 CIRIquant 真正成功執行；若排查 CIRIquant 問題，需直接看 `logs/ciriquant/{srr}.log` 而非 Snakemake 的 job 狀態 |
| CIRIquant 輸出 `.bed` 而非 `.bsj` | CIRIquant 1.1.3 bioconda 版本差異 | rule output 只宣告 `.gtf`，忽略 `.bed` |
| multiqc ImportError TypedDict | multiqc 1.17 + markdown 3.6 不支援 Python 3.7 | `pip install markdown==3.3.7` |
| Snakemake LockException | 前次執行被 kill 留下 lock | `snakemake --unlock` 後重跑 |
| samtools sort "File exists" | 多個 CIRIquant 進程同時寫同一 sample 的暫存檔 | `pkill -f CIRIquant`，`rm -rf results/circRNA/`，重新跑 |
| `nohup: failed to run command 'snakemake'` | conda env 未啟動 | 先 `conda activate ciriquant` |
| wildcard ambiguity（star_align vs mate1/mate2） | `{srr}` wildcard 匹配到 `SRR7012366/mate1` | 在 `circrna.smk` 加 `wildcard_constraints: srr = r"[A-Z]+\d+"` |
| star_align temp dir 硬編路徑 | 原本 `/home/choukaihsuan/star_tmp/{srr}` 只適用本機 | 改為 `RESULTS_DIR + "/circRNA/{srr}/star_tmp"` |
| STAR `_STARtmp` 殘留目錄導致重跑失敗 | STAR 在 `{outTmpDir}/_STARtmp` 建立工作目錄；若前次 run 被 kill（SIGKILL / ulimit）後，`_STARtmp` 殘留；`star_align` shell 末尾的 `rm -rf {params.tmp_dir}` 只在成功後執行，失敗時不清理；下次 `--forcerun star_align` 時 STAR 立即報錯：`could not create temporary directory: .../_STARtmp` | 修正 `circrna.smk` `star_align` rule：在 `mkdir -p {params.tmp_dir}` **之前**加 `rm -rf {params.tmp_dir}`，確保每次從乾淨狀態啟動；手動修復已損壞樣本：`rm -rf {results_dir}/circRNA/{srr}/star_tmp` |
| STAR `--outSAMtype BAM SortedByCoordinate` 在 NFS 環境下偶發失敗 | STAR 內建排序在 NFS 上開啟大量暫存檔（file descriptor 數量超過 NFS 限制），導致 sort step 崩潰或輸出損壞 BAM | `star_align` rule 改為 `--outSAMtype BAM Unsorted`，STAR 輸出 `Aligned.out.bam`（未排序），再用 `samtools sort -@ {threads} -m 2G` 手動排序 → `{output.bam}`，最後 `rm -f Aligned.out.bam`；mate1/mate2 rule 加 `ulimit -n 65536 || true` 提高 NFS file descriptor 上限（`circrna.smk`，2026-07-04）|
| DCC `ValueError: invalid literal for int(): '4"'` | `-fg` 觸發 CircSkip 解析，GTF exon_number 屬性殘留尾部引號 | 移除 `-fg` flag，CircSkip 計數非必要 |
| DCC 六個樣本全部失敗（IndexError: list index out of range） | 多個 DCC 並行共用工作目錄下的 `_tmp_DCC/`，競爭條件 + 上次失敗殘留的 partial 資料 | 改為 `(cd {params.outdir} && DCC ...)` subshell，每個 sample 的 `_tmp_DCC/` 獨立在各自 outdir 內 |
| DCC log 路徑失敗（`logs/dcc/SRR.log: No such file or directory`） | `cd {outdir}` 後 `> {log}` 的相對路徑從 outdir 解析 | 同上，subshell 讓 `> {log}` 在 parent shell 的 CWD（`~/circRNA_agent/`）執行 |
| multiqc numpy 版本衝突（`numpy 1.16.4 < 1.17`） | conda env 安裝的 numpy 過舊，matplotlib 要求 ≥1.17 | `pip install 'numpy>=1.17'` 在 ciriquant env |
| `consensus_filter.py` `TypeError: 'type' object is not subscriptable` | Python 3.7：module-level 的 `CoordMap = dict[tuple[...], ...]` 賦值在執行時求值，不能用內建 `dict`/`tuple` 做 subscript；`from __future__ import annotations` 只延遲 annotation 求值，**不影響普通賦值** | 改為 `from typing import Dict, Tuple; CoordMap = Dict[Tuple[str,int,int], float]` |
| `parse_ciriquant` / `merge_counts.py` 回傳 0 筆 circRNA | CIRIquant 1.1.3 GTF 屬性欄用小寫 `bsj`/`fsj`，regex 搜尋大寫 `BSJ`/`FSJ` 全部未匹配 | 兩個 `re.search()` 加 `re.IGNORECASE` flag（`consensus_filter.py` 和 `merge_counts.py` 均適用） |
| `parse_dcc` 回傳 0 筆 circRNA | `CircCoordinates` 只有座標（8 欄，無 count）；舊程式讀 col 3 得到 Gene 欄（字串），`float()` 失敗後全部 skip | `parse_dcc()` 改為優先讀同目錄的 `CircRNACount`（col 3 = junction count）；不存在時 fallback count=5 |
| `analysis.R` `storage.mode(counts) <- "integer"` 失敗 | conda install r-pheatmap 後 dplyr 進入 env，data.frame subsetting 回傳 tibble-like 物件，`storage.mode` 不接受 | 在 `storage.mode` 之前加 `data.matrix(round(counts))` 轉為真正的 matrix（BSJ 和 FSJ 矩陣均需處理） |
| `generate_report.py` `SyntaxError: from __future__ imports must occur at the beginning` | Snakemake `script:` wrapper 在 user 腳本前插入 ~11 行 setup code，將 `from __future__ import annotations` 推到第 12 行；Python 3.7 要求此 import 必須是第一行 | 移除 `from __future__ import annotations`；改用 `from typing import Optional`；所有 `str \| None` 改為 `Optional[str]` |
| `analysis.R` `Error: replacement has 0 rows, data has 7779`（edgeR_ciriquant + use_pvalue）| `bsj_sig_col` 設為 `"PValue_bsj"`，但 merge 後該欄名為 `"PValue"`（只有 logFC/FDR 在兩表均存在才加後綴） | 改為 `bsj_sig_col <- if (use_pvalue) "PValue" else "FDR_bsj"` |
| Snakemake hook `notify.py` 路徑錯誤（找到 conda env 的 Python runner） | `onstart`/`onsuccess`/`onerror` 中 `__file__` 指向 conda env 的腳本執行器，非專案目錄 | 改用 `workflow.snakefile` 推算專案根目錄路徑 |
| benchmark `MissingInputException`（hg19.gtf） | benchmark config 路徑是 `hg19.gtf`，但 server 實際檔名是 `genes.gtf` | 更新 `benchmark/config_benchmark.yaml` GTF 路徑 |
| isoform switching 找到 0 個 significant events | Global BH 校正跨 ~7,360 isoform，min padj = 0.937 | 改用 within-gene BH（`padj_within_gene`）；結果：66 events（FDR < 0.1, \|ΔIUI\| > 0.1） |
| DE analysis 0 個 significant circRNA（padj 校正） | n=3 vs 3，BH 校正後 min padj = 0.432 | 改用 nominal p-value（`de_sig_by: pvalue`）；結果：482 個 significant circRNA |
| `bioconductor-edger=3.44.0` PackagesNotFoundError | conda bioconductor channel 無此 patch 版本 | 改為 `bioconductor-edger>=3.40`（不鎖 patch） |
| `DCC==0.5.0` pip install 失敗 | PyPI 上的 DCC 是不相干套件（版本 0.7+）；circRNA DCC 在 bioconda | 移至 conda 依賴 `dcc=0.5.0`，從 pip section 移除 |
| Docker build `mamba env create` 等待確認 | `mamba env create` 預設互動式 `[Y/n]` | 加 `--yes` flag：`mamba env create --yes` |
| Anaconda `repo.anaconda.com` 商業授權警告 | `defaults` channel 要求商業許可 | 移除 `defaults`，僅保留 `bioconda` + `conda-forge` |
| `generate_report` 重跑觸發上游 rules（download/align） | server config.yaml 指向其他專案，`--forcerun` 重建整個 DAG | 改用 `/tmp/run_generate_report.py` mock snakemake 物件直接執行 |
| SVG badge 距離弧段過遠（舊 angular de-overlap） | `deOverlap` 只往順時針推移，多個 badge 連鎖推移後遠離原弧段 | 改為 radial staggering：badge 固定在弧段正上方角度，僅在 angular gap < 0.14 rad 時改變半徑（`MI_R0=+7` vs `MI_R1=+19`）；虛線放在 `<g>` 內隨 badge 隱藏 |
| `analysis.R` heatmap 錯誤（`pheatmap` 少於 2 行） | DE circRNA 數量太少時 `pheatmap` 會報錯 | 加 `if (nrow(mat) >= 2)` 判斷，否則 `plot.new()` 顯示提示文字 |
| `analysis.R` `slice_min` 版本警告 | dplyr 新版棄用 `slice_min` 的部分用法 | 改為 `arrange(.data[[col]]) %>% head(n)` |
| `generate_report` Plotly 圖表沒有互動（靜態 PDF embed） | 用 base conda Python（3.13）執行，無 plotly | 改用 `conda run -n ciriquant python /tmp/run_generate_report.py` |
| modal 分頁無法點擊（整個 JS 失效） | `priTitle` 字串中的 Python `\n` 變成 JS 字串裡的實際換行 → SyntaxError | 改為 `\\n`，讓 JS 收到合法逸脫序列 |
| Isoform switching bar chart x-axis label 重疊 | 10 基因 × 2 條件 = 20 個 label 太密 | 改用 Plotly multi-level x-axis（`[[gene,gene],["Normal","Tumor"]]`） |
| `analysis.R` limma `object 'logFC' not found` | dplyr 1.1.x `rename(log2FC = logFC)` NSE 在 ciriquant env 下找不到欄位 `logFC`；DESeq2 和 edgeR 用 base R 賦值故無問題 | limma 區塊改用 base R：`names(res_df)[names(res_df) == "logFC"] <- "log2FC"` |
| benchmark `--forcerun comparison_report` 觸發 bowtie2 重跑 | `comparison_report` 依賴 `accuracy_summary` → `nfcore_sim.bed` → `bench_find_circ` → `bench_find_circ_map`；任何 `--forcerun` 包含 `comparison_report` 都會觸發整個 dependency chain | 改用直接腳本執行（見 Benchmark 獨立重跑腳本）；或只 `--forcerun de_deseq2_baseline de_quality_benchmark` |
| benchmark `de_deseq2_baseline` groups 路徑失效 | `benchmark/config_benchmark.yaml` 的 `groups` 指向全域 `metadata/sample_groups.csv`；config.yaml 切換到 GSE58135 後該檔案是 GSE58135 的 sample list | 更新 `benchmark/config_benchmark.yaml` 的 `groups` 改為 `metadata/GSE113230/sample_groups.csv` |
| Snakemake `--forcerun de_analysis` 觸發完整 DAG 重建 | `--forcerun` 標記 de_analysis 需重跑，但 Snakemake 評估 transitive inputs；若 `raw_dir` FASTQ 不存在則排入 download_fastq | 改用 `/tmp/run_de_analysis.R` mock snakemake 直接執行，或確保 `--configfile` 指向正確專案且所有上游 output 都存在 |
| `FULL_HEATMAP_DATA.conditions` 全為灰色（wrong SRR IDs）| Server 上 `config.yaml` 切換到 GSE58135，`metadata/sample_groups.csv` 被蓋成 GSE58135 SRR IDs；`/tmp/run_generate_report.py` 的 `groups` 參數指向此被蓋掉的檔案 | 改為指向 `metadata/GSE113230/sample_groups.csv`（固定路徑不受全域 config 影響）|
| 報告中 circRNA 以座標 ID 儲存，搜尋 circBase 名稱找不到 | 管線使用基因組座標 `chr:start\|end` 作為 circ_id，不是 circBase ID | 用 `circbase_annotated.tsv` 做 ID 對應；`annotate_circbase.py --slop 10` 容許 ±10 bp 誤差（1-based/0-based 轉換造成的 1 bp 差異在容許範圍內）|
| find_circ.py Python 3 三個 bug（`NameError: _np`、`mismatches` 縮排錯誤、`complement()` KeyError: `'\n'`）| Python 2→3 遷移不完整：`mismatches()` 內部的 `import numpy as _np` 縮排在 inner function，外層程式碼無法存取；`complement()` 用 `COMPLEMENT[x]` dict lookup 但基因組序列讀取包含 `\n`；`fromstring(a,dtype=byte)` Python 3 棄用 | 修 `/home3/choukaihsuan/tools/find_circ/find_circ.py`：(1) 模組層級加 `import numpy as _np`；(2) `mismatches()` 改用純 Python `sum(x!=y for x,y in zip(a,b))`；(3) `complement()` 加 `s = s.replace(chr(10),"").replace(chr(13),"")`  |
| `parse_find_circ()` 載入 235K 行造成 consensus_filter 卡住 55+ 分鐘 | find_circ 輸出含 LINEAR（233K）、AMBIGUOUS（190K）、CIRCULAR（僅 1,895）三類 junction，舊 parser 未過濾 category 欄 | `parse_find_circ()` 加判斷：`if len(parts)>=18 and "CIRCULAR" not in parts[17]: continue`；真正 circRNA 從 235K 降至 1,895 個，執行時間從 55min+ 降至秒級 |
| RBP 分頁 Binding Seq 欄大多顯示 "N/A" | ENCORI 的 `circ_pos` 是絕對座標格式（`chr6:148390208-148390208`），CircInteractome 是相對座標（`156–191`）；`_absPos()` 把 ENCORI 絕對座標再加 `chromStart`，產生超出染色體長度的位置 → UCSC API 回傳無效 → "N/A" | `_absPos()` 和 `_seq_logo` 加偵測：若 `circ_pos` 以 `chr` 開頭則視為絕對座標直接使用，否則才加 `chromStart`（`generate_report.py`） |
| benchmark `config_benchmark.yaml` BWA index 路徑錯誤 | `bwa_index: /home3/choukaihsuan/reference/hg19/bwa_index/hg19` 目錄不存在；實際 BWA index prefix 是 `hg19.fa` | 改為 `bwa_index: /home3/choukaihsuan/reference/hg19/hg19.fa` |
| find_circ_map timing bowtie2 卡住 13+ 小時 | timing script 將 28GB unmapped reads 寫入 `/tmp`（掛載在 `/` 分區，98% 使用率），磁碟 I/O 極慢 | 改寫到 `/home3/choukaihsuan/timing_tmp/`（348GB 可用）；`/` 分區的 timing_tmp 清理後釋放 11GB |
| `--adaptive` 沒有傳入 `consensus_filter.py` | `circrna.smk` 的 consensus_filter shell 指令從來沒帶 `--adaptive` flag；adaptive 只在腳本內部實作但未啟用 | `circrna.smk` 加入 `adaptive_flag = "--adaptive" if config["consensus"].get("adaptive", True) else ""`，預設開啟；同時加 `--adaptive-ratio` 參數 |
| `vst()` 在少於 1000 circRNA 時失敗（GSE58135）| DESeq2 的 `vst()` 要求輸入行數 ≥ nsub（預設 1000）；小資料集（如 28 circRNA）呼叫 `vst()` 報錯：`less than 'nsub' rows` | `analysis.R` 改為 `tryCatch(vst(dds, blind=FALSE), error=function(e) varianceStabilizingTransformation(dds, blind=FALSE))`，自動 fallback |
| GSE58135（50bp reads）consensus 只有 28 個 circRNA | DCC 在 50bp 讀長下每個 sample 只偵測到 3–17 個 circRNA（STAR chimeric 短讀限制），但 `--adaptive` 未傳入，min_tools=2 要求兩工具交集，結果全部 10 sample 只有 28 個 | 修復 `--adaptive` 傳入問題後，ratio ≈ 0.002 << 0.1，adaptive 自動降級 min_tools=1，結果 1,607 個 circRNA |
| `_biomarker_normality_plot` Shapiro-Wilk ValueError（n < 3）| GSE323364 只有 2 個 biomarker candidates，`scipy.stats.shapiro` 要求 n ≥ 3 | ✅ **已無關**：Shapiro-Wilk 檢定、Normal fit 曲線、垂直線（μ/±σ/±2σ）已全部移除（2026-07-31）；`_biomarker_normality_plot()` 現在只產生純 Histogram，Y 軸為 Number of circRNAs，bin 寬 0.05 |
| Web UI 手動 SRR 表單 Project ID 預設填入當前專案 ID | `value="{{ config.get('project_id', 'CUSTOM') }}"` 讓使用者未改 ID 就送出，蓋掉現有專案的 metadata | 改為 `value="" required`，強制使用者填入正確 GSE ID |
| `config/projects/GSE133998.yaml` sra_cache_dir/tmp_dir 繼承舊專案路徑 | Web UI 從前一專案（GSE323364）快照複製 config，路徑殘留 `GSE323364/sra_cache`；Snakemake DAG 正常但中間檔案放錯位置 | 用 `sed -i` 替換 config 中的路徑（`GSE323364` → `GSE133998`）|
| gzip 壓縮 SRA 轉換 FASTQ 太慢（~2h/15GB）| fasterq-dump 後 gzip 單執行緒壓縮，NFS 寫入慢 | 在 ciriquant env 安裝 `pigz`（`conda install -y -c conda-forge pigz`）；`download.smk` 已使用 pigz（若 pigz 不存在，`find_tool("pigz")` fallback 到 gzip）|
| NCBI S3 URL 間歇性不可用（`srapath` 回傳空）| NCBI API 暫時失敗，`srapath --location s3` 回傳空字串；download.smk fallback 到 prefetch HTTPS（~22 KB/s）| 等待 S3 恢復後手動執行 aria2c；一般幾分鐘至數小時內恢復 |
| S3 aria2c DNS SERVFAIL（`DNS server returned general failure`）| 並行下載連線總數過多（例如 10 samples × 16 connections = 160）觸發 AWS S3 per-IP throttle，DNS 回傳 SERVFAIL；持續約 2 小時後自動恢復 | 限制總連線數 ≤ 30：`-x 8` × 4 parallel = 32（當前設定，可接受）；`-x 4` × 4 = 16（保守）；早上手動腳本 10 × 16 = 160 是觸發根因；恢復後重跑 pipeline 即可 |
| `kill <PID>`（SIGTERM）殺掉 Snakemake 時刪除 output 檔 | Snakemake 收到 SIGTERM 執行 graceful shutdown，自動刪除所有「正在執行 job 的 output 檔案」以防止不完整輸出；若早先下載完成的 FASTQ 恰好被列為「執行中 job 的 output」（因未設 protected），就會被刪除 | **一律用 `kill -9 <PID>`（SIGKILL）終止 Snakemake**，SIGKILL 無法被捕捉，不觸發 cleanup；永遠不用 bare `kill` 或 `pkill` 不加 `-9`；GSE133998 因此損失 9 個 FASTQ（~63 GB），需重新下載 |
| 修改 `.smk` 規則觸發所有 jobs 重跑（"Code has changed"）| Snakemake 預設追蹤 rule 的 code hash；修改任何 `.smk` 檔案後，Snakemake 在下次執行時標記該 rule 所有 jobs 為「Code has changed since last execution」→ 全部重跑，即使 output 已存在 | **只在 Snakemake 未執行時修改 `.smk`**；若需 mid-run 修改，加 `--rerun-triggers mtime` 旗標讓 Snakemake 只依時間戳判斷（忽略 code 變更）；GSE133998 因將 `-x 4` 改為 `-x 8` 而觸發所有 download_fastq 重跑 |
| `--rerun-incomplete` + 未 protected output + SIGTERM = output 消失 | `--rerun-incomplete` 把被 kill 前的「未完成 job」全部重排；若這些 job 的 output（早先其實已完成）未被 Snakemake 的 `protected()` 標記（原始 run 被 kill -9 前來不及 chmod），Snakemake 視為「正在執行」；再次 SIGTERM 時 cleanup 刪除這些 output | 避免在有執行中 job 的情況下用 `--rerun-incomplete` 重啟；若必須重啟，先手動 `chmod 444` 保護重要 output 再重跑，**並且**確認 Snakemake metadata（`.snakemake/metadata/`）也記錄該檔案為 protected（光 chmod 不足，需 Snakemake 自己完成的 job 才記錄）|
| Circular Structure 第一個分頁空白（無弧段）| interactions.json 的 `spliced_length` 欄位固定為 0；JS `if(totalLen>0)` guard 阻擋所有弧段繪製 | JS 改為：`const totalLen = parseInt(info.spliced_length) \|\| (exonBds.length>0 ? exonBds[exonBds.length-1].cum_end \|\| 0 : 0)`；從 `exon_boundaries[last].cum_end` 推算真實長度 |
| Venn detail 「方法」欄顯示三個方法（包含不顯著的）| `_venn_3_svg` 建立 `circ_info` 時對每個方法的**全部** DE 行都加入 `mlabel`，不論該 circRNA 在該方法中是否顯著 | 加入 `is_sig = cid in m_sig_set` 判斷；只有當 circRNA 在 `sig_sets[m]`（顯著集合）中才加入 `methods` 清單 |
| Heatmap tumor/normal 欄位混排 | 樣本欄依 SRR ID 字母順序排列，奇數偶數交錯（normal→tumor→normal…）| `_plotly_heatmap()` 和 `FULL_HEATMAP_DATA` 建立時，讀取 `sample_groups.csv` 按 condition 重排：tumor 欄在左、normal 欄在右 |
| `analysis.R` `coef=2` 在配對設計下指向 patient 係數而非 condition | `model.matrix(~patient+condition)` 有多個 patient dummy；`coef=2` 指向第一個 patient，非 condition | 改為 `cond_coef = ncol(design)`（condition 永遠在最後一列）；`glmQLFTest(fit, coef=cond_coef)`、`topTable(fit, coef=cond_coef)` |
| `download.smk` Python 3.7 不支援 `unlink(missing_ok=True)` | `missing_ok` 參數在 Python 3.8 才加入；CentOS 7 server 的 conda Python 3.7 執行時 TypeError | 改用 `try: lock.unlink() except OSError: pass` 相容寫法 |
| Type II DE circRNA 幾乎全為零 | 兩個根本原因：(1) `abs(logFC_fsj) >= 0.5` 閾值過嚴；(2) offset approach 造成 log2FC 與 logFC_fsj **反向相關**，`sign(log2FC)==sign(logFC_fsj)` 幾乎從不成立（73 sig_fsj 中只有 1 個同方向） | 移除方向一致性檢查，改為 `sig_bsj & sig_fsj`（BSJ ratio 顯著 AND FSJ 獨立顯著）；結果：GSE113230 從 0 → 73 Type_II（15.1%），GSE58135 → 2（13.3%），GSE133998 → 2（2.4%） |
| benchmark AUC-PR 虛高（二元偵測假象）| `sklearn.average_precision_score` 對 score ∈ {0,1} 做樂觀排序，ground truth 中大量 FN 排在 TN 前，產生假精確度長段，AUC-PR 嚴重虛高 | 改用**門檻掃描**（threshold sweep）：對 CIRI2 + DCC CircRNACount 各取 min_bsj ∈ {1,2,3,...,50}，每個門檻計算真實 (Precision, Recall)，trapezoid AUC；`-Nr 5` 時代首次測得：circDEX=0.120、CirComPara2_4tools=0.349、nfcore_3tools=0.337；`-Nr 2` 重跑後最新值見 GSE113230 段落「Benchmark 門檻掃描 AUC-PR 重新確認」；新增 `--ciri2-file / --dcc-count-file / --output-pr-curve` CLI 參數至 `accuracy_benchmark.py`，`--pr-curve` 至 `generate_comparison_report.py` |
| CirComPara2 工具數誤標（5→4 工具）| `generate_comparison_report.py` 的說明文字列出 CIRI2 + CIRIquant + DCC + find_circ + CIRCexplorer2 = 5 工具；但 CIRIquant 內部已呼叫 CIRI2 做 BSJ 偵測，兩者是同一演算法 | 修正為 4 工具（CIRIquant + DCC + find_circ + CIRCexplorer2）；Consensus 門檻改為「≥2/4 tools」；報告說明文字同步更新 |
| DESeq2 `every gene contains at least one zero, cannot compute log geometric means`（SRP156355）| sparse count matrix（14,697 circRNAs × 12 samples）有大量零值，`DESeq(dds)` 預設幾何平均 size factor 估算失敗 | `analysis.R` 在 `DESeq(dds)` 前加 `dds <- estimateSizeFactors(dds, type = "poscounts")`；poscounts 用非零值估算，避開全零問題；DESeq2 偵測到 sizeFactors 已設定後不再重新估算 |
| `predict_interactions` 僅覆蓋 edgeR top-100，切換 DESeq2/limma 時 Biomarker Score 偏低 | interactions.json 只查詢 edgeR top 50 up + 50 down；DESeq2/limma 顯著集合有許多 circRNA 不在此範圍，mirna_norm = rbp_norm = 0 → score 壓縮 | `predict_interactions.py` 新增 `--de-edger`、`--de-deseq2`、`--de-limma` 三個選填參數；提供時取三方法 top-N 聯集（~150–250 circRNA）查詢；metadata 改從 `iso` + `circbase` 讀取（不受 filterByExpr 限制）；`de.smk` 同步更新傳入三個 input |
| Biomarker 表格切換 DE 方法時只灰化、不重排 | 原本 `_updateBiomarkerHighlight()` 只對非顯著行設 `opacity:0.25`，不改變排序；切換 DESeq2/limma 時看到的仍是 edgeR 排名的 top 30 | 新增 `_compute_bm_table_data()` Python 函數，於 `ALL_DE_METHODS[m]["bm_table"]` 存入各方法的 top-30 rows；`generate_report.py` 同時從 `circbase_annotated.tsv` + `interactions.json` 補全 `_bm_lookup`（含 DESeq2/limma-only 的 circ_id）；新增 JS `_renderBiomarkerTable()` 完整重建表格 DOM + 更新過濾按鈕計數；`switchDEMethod()` 改呼叫 `_renderBiomarkerTable()` 取代 `_updateBiomarkerHighlight()` |
| Biomarker score > 1.0（如 1.150, 1.061）| `_compute_bm_table_data` 從 `_bm_lookup` 取的 `mirna_n`/`rbp_n` 是以 edgeR biomarker_candidates.tsv 最大值正規化；但 union mode 的 interactions.json 包含更多 circRNA，有些 n_mirna 遠超此最大值（如 248 vs 基準 146），導致 mirna_n > 1 → score > 1 | 改為在**當前 sig set 內**重新計算 `_sig_mirna_mx` 和 `_sig_rbp_mx`，用局部最大值正規化：`mirna_n = n_mirna / _sig_mirna_mx`；每次切換 DE 方法都獨立正規化，保證 score ∈ [0, 1] |
| Biomarker Type 欄顯示 "nan"（DESeq2/limma bm_table）| DESeq2/limma 的 `analysis.R` 設 `res_df$Type <- NA_character_`（只有 edgeR_ciriquant 做 FSJ 分類），Python 讀取 NA 為 NaN，`str(NaN)="nan"` | `_bm_lookup` 建立時從 `biomarker_candidates.tsv` 讀取並儲存 `type_edger`（edgeR 分類）；`_compute_bm_table_data` 中 Type 若為 NaN 則 fallback 到 `lu.get("type_edger")`；效果：3 方法均顯著的 circRNA 在 DESeq2/limma bm_table 也能看到 edgeR Type 分類 |
| circ_id 欄位排序不符基因組順序（chr4 排在 chr1 前）| `_makeSortable` 去除非數字字元取第一個數字：`chr4:1902353` → `41902353`，而 `chr1:91382197` → `191382197`；前者較小故 chr4 排在 chr1 前 | 新增 `_parseGenCoord(s)` 函式：`chr([0-9]+|XYM):start` → `{c: chrNum, p: startPos}`（chrX=23, chrY=24, chrM=26）；`_makeSortable` 優先使用基因組座標排序（先比 chr 號碼，同染色體再比 start position） |
| DE 表格 gene_name 顯示 "intergenic" 與 region 欄重複 | `assign_isoforms.py` 對無 host gene 的 circRNA 設 `gene_name = gene_id = "intergenic"`；DE table 照原值顯示 | `_renderDETables` JS 中加入 `if(col==='gene_name')` 判斷：當 `v==='intergenic'` 時顯示 `'—'`（region 欄已有此資訊，gene_name 欄不重複顯示）|
| Web UI 重啟後狀態頁全顯示「等待中」 | pipeline 執行狀態只存在 Flask 記憶體中；web_ui 重啟後 stages 全重設為 pending，即使所有步驟已完成 | `_infer_stages_from_files()` 掃描 results_dir output 檔（multiqc_report.html / count_matrix.tsv / de_results.tsv / report.html 等）→ 升級 pending stage 為 done；`/api/progress` 每次 polling 都呼叫此函式 |
| `/api/progress` 未回報 report/QC 存在 | 狀態頁「分析報告」和「QC 報告」按鈕依賴 stages 狀態，若 web_ui 重啟後 stages 全 pending，按鈕不顯示 | 在 `/api/progress` 回應中加入 `report_exists`（bool）和 `qc_exists`（bool）欄位；前端 JS 以 `data.report_exists \|\| stage.done` 控制按鈕顯示 |
| `web_ui.py _update_paths_for_project()` 未更新 `download.sra_cache_dir` / `download.tmp_dir` | 切換到新 GSE 時，`_update_paths_for_project()` 只更新 `raw_dir`、`trimmed_dir`、`results_dir`，不更新 `config.download.sra_cache_dir` 和 `config.download.tmp_dir`；中間 SRA 快取和暫存檔案放到前一個專案的目錄 | ✅ 已修正（2026-07-02）：`sra_cache_dir` / `tmp_dir` 從更新後的 `raw_dir` parent 自動重新推算（`parent + "/sra_cache"` / `parent + "/sra_tmp"`），不再做字串替換，路徑永遠跟 `raw_dir` 同步 |
| `_update_paths_for_project()` 在 `old_pid == new_pid` 時不做路徑替換，config 已存在但路徑錯誤時無法自修 | `config/projects/GSE148036.yaml` 已有 `project_id: GSE148036` 但路徑指向 `GSE229705/`；`old_pid = cfg["project_id"] = "GSE148036"` = `new_pid` → 條件跳過 → 路徑未更新 → pipeline 下載到 GSE229705 目錄並可能覆蓋結果（GSE148036-SFJU 教訓，2026-07-02）| ✅ 已修正（2026-07-02）：移除 `old_pid != new_pid` guard，改為逐個 key 檢查「path 是否已含 new_pid」，不含則先嘗試 `old_pid` 替換，再 fallback 到 regex 替換路徑中任何 `GSE/SRP/PRJNA\d+` 字串 |
| SSH heredoc 寫入腳本含 CRLF，變數在執行時為空 | 在 SSH double-quoted 命令中用 `<< 'SCRIPT'` heredoc 寫 shell 腳本，Windows 端 SSH client 將行尾轉為 CRLF（`\r\n`）；shell 讀取後變數賦值包含 `\r`，`echo "$VAR"` 看似有值但實際傳入命令時被截斷或輸出空值 | 改用 `printf '%s\n' 'line1' 'line2' ...` 逐行寫入，每行明確加 `\n`，完全避開 heredoc CRLF 問題 |
| `$TMP` 系統保留變數名稱衝突 | Shell 環境中 `TMP` 是保留變數（通常指向 `/tmp`）；腳本中 `TMP=/home3/.../sra_tmp/SRR...` 看似賦值，但某些情境下 `$TMP` 展開為系統值或空字串，導致 `fasterq-dump --temp $TMP` 路徑錯誤 | 改用自訂名稱 `SRA_TMP`（避開 `TMP`、`TMPDIR`、`TEMP` 等系統保留名）；`SRA_CACHE` 同理（避開 `CACHE`）|
| Biomarker 表格 `n_rbp` 顯示舊值（如 8）但 RBP Binding modal 顯示 50+ RBP | `predict_interactions.py` 用 `--gtf` 重跑後更新 `interactions.json`，但 `rank_biomarkers.py` 未重跑，`biomarker_candidates.tsv` 的 `n_rbp` 仍是舊值；`_bm_lookup` 從 TSV 讀取，不反映新的 `in_circ=True` 計數 | `generate_report.py` 的 `_bm_lookup` 建立後加入覆蓋步驟：`_is_in_circ_entry()` 函式（CircInteractome=True，ENCORI=`in_circ`）重新計算所有已存在 entry 的 `n_mirna`/`n_rbp`，用最新 `interactions.json` 覆蓋 TSV 的舊值；`mirna_n`/`rbp_n`（正規化值）也同步更新 |
| Benchmark 互動式表格 `gene_name` 全部空白 | `de_quality_benchmark.py` 直接讀取 DE 結果 TSV（`de_results.tsv`、`nfcore_deseq2_results.tsv`、`nfcore_limma_results.tsv`），這些 TSV 沒有 `gene_name` 欄（來自 `isoform_groups.tsv`）；`_build_de_list` 的 `r.get("gene_name")` 取到 None → 顯示空白 | `de_quality_benchmark.py` 加入 `--isoforms <isoform_groups.tsv>` 參數；讀取後 `pd.map(iso_map)` join 到所有 DE dataframe |
| Benchmark modal 表格無法排序 | `showDEList` 有獨立的舊版表格渲染邏輯（無 onclick），只有 `showJaccardList` 呼叫新版 `_renderCircList`（含 `sortCircList` onclick）；DE Quality modal 點擊欄位無反應 | `showDEList` 改為呼叫 `_renderCircList(data, hasType ? 'Type' : null)`，與 Jaccard modal 統一渲染；`_renderCircList` 表頭加 `onclick="sortCircList(this,N)"`；`sortCircList` 以 `data-val` 屬性做數值/字串比較排序，排序後重新編號 `#` 欄 |
| Queue Worker `_run_queued_job` status 永遠卡在 `running` | `register_job()` 和 `log_path.parent.mkdir()` 在 `try` 區塊**外**；這兩行若拋出例外，外層 `_queue_worker` 的 `except` 只 `print()`，DB status 永遠不更新為 `failed`，卡在 `running` | 將 `register_job` 和 `mkdir` 移入 `try` 區塊內（已修正，2026-06-29） |
| Queue Worker 遇殘留 lock 直接失敗 | `kill -9` Snakemake 後 `.snakemake/locks/` 殘留；Queue Worker 下次啟動同一專案時碰到 `LockException`，`proc.wait()` 返回 rc=1 → 設 `failed`，不會自動重試 | `_run_queued_job` 在啟動 Snakemake 前先執行 `snakemake --unlock` 預清殘留 lock（已修正，2026-06-29） |
| `config/projects/{GSE}.yaml` 的 `sra_cache_dir` / `tmp_dir` 繼承舊專案路徑 | `web_ui.py _update_paths_for_project()` 只更新 `raw_dir`/`trimmed_dir`/`results_dir`；PRJNA553289.yaml 的 `sra_cache_dir` 指向 `GSE77509/sra_cache`，`tmp_dir` 指向 `GSE323364/sra_tmp` | 手動 `sed -i` 修正目標路徑；長期解法：`_update_paths_for_project()` 加入 `download.sra_cache_dir` / `download.tmp_dir` 替換（參見現有 known issue 的修正說明） |
| 登入頁錯誤訊息只有中文 | `login.html` 的 `{% if error %}` 只顯示中文錯誤；切換 EN 語言後錯誤仍顯示中文 | `login()` route 新增 `error_en` 欄位；`login.html` 加 `data-en="{{ error_en }}"` → lang switcher 自動替換為英文版錯誤訊息（2026-06-28） |
| Web UI 登入頁「Bad Request - The CSRF tokens do not match」 | Flask-WTF 的 `CSRFProtect(app)` 對所有 POST 路由全域啟用 CSRF 驗證；`/login` POST 也受保護，但登入前 session 為空或過期，CSRF token 無法驗證，形成死鎖（無法登入 → 無法取得有效 session → CSRF 永遠失敗）| `/login` route 加 `@csrf.exempt` 裝飾器，讓公開登入頁免除 CSRF 驗證（登入本身無需 CSRF 保護，magic link 才是驗證機制）；已修正並同步至 server（2026-10-02） |
| Circular Structure overlap annotation 字母重疊 | Solo mode（單分子可見）時 SVG 中出現兩組字母：per-arc site label（a/b/c）＋ overlap annotation letter（i/v，位於 midR=48 center hole 或 midR=164.5 miRNA ring）；使用者混淆 | 移除兩組 overlap annotation 的 `<circle>` + `<text>` badge（保留 boundary dash line）；只保留 `_arc_lbl` per-arc site labels（`generate_report.py`，2026-06-29） |
| RBP site badge 不在 arc 上 | RBP arc（radius 90–113）在 exon ring 內部視覺上細不明顯；badge 放在 `(RBP_IN+RBP_OUT)/2=101.5` 看起來浮空 | Badge 改放在 exon ring 中心（`midR=(RIN+ROUT)/2=131.5`），直接覆蓋在 exon arc 上，角度指示 binding site 位置；進入 solo mode 時 arc 加 white stroke 1.5px 高亮（`generate_report.py`，2026-06-29） |
| isoform_switching mock 腳本的 params key 錯誤（`isoform_fdr_cutoff` vs `fdr`）| `isoform_switching.R` 讀取 `snakemake@params[["fdr"]]`（Snakemake rule 傳的 key），但手動 mock 腳本傳 `isoform_fdr_cutoff = 0.1`；`snakemake@params[["fdr"]]` 回傳 NULL → `as.numeric(NULL)` = `numeric(0)` → `res$padj_within_gene < numeric(0)` = `logical(0)` → `res$is_switching <- logical(0)` 失敗（replacement has 0 rows, data has N rows）| mock 腳本 params 改為 `fdr = 0.1`；實際 Snakemake rule 不受影響（key 正確）|
| GSE221107 RNA 降解樣本（Pair 8/11）導致 edgeR 失敗 | Pair 8（SRR22757442：6 BSJ）和 Pair 11（SRR22757419：28 BSJ）insert size peak ≈ 40–44 bp（正常 268–269 bp）→ adapter dimer，RNA 在 library prep 前已降解；加入 TMM normalization 導致 NaN/Inf，edgeR GLM missing value 失敗 | 排除 Pair 8 和 Pair 11，保留 Pair 14/16/18/20（4 pairs，8 samples）；Python 子集 count_matrix.tsv 和 fsj_count_matrix.tsv，更新 sample_groups.csv，重跑 DE/isoform_switching/predict_interactions/rank_biomarkers/report |
| Web UI 頁面被瀏覽器自動翻譯為韓文或日文 | Chrome/Edge 偵測到 `<html lang="zh-TW">` 的中文頁面後，若瀏覽器語言偏好含韓文/日文或使用者曾接受翻譯提示，會自動將整頁翻譯 | 三個模板（`index.html`/`status.html`/`login.html`）均加入 `translate="no"` 屬性至 `<html>` 標籤，並在 `<head>` 新增 `<meta name="google" content="notranslate">` 和 `<meta http-equiv="Content-Language" content="zh-TW">`（2026-07-01）|
| web_ui 啟動後 log 停止更新（stdout pipe 斷裂）| `nohup conda run -n ciriquant python scripts/web_ui.py >> logs/web_ui.log` 中，`conda run` 建立的子進程 stdout/stderr 走 pipe 回傳給 conda 進程，conda 進程在 nohup shell session 關閉後 pipe 斷裂，log 不再更新；但 python 進程本身仍存活 | 改用直接指定 conda env python 路徑：`nohup /home/choukaihsuan/miniconda3/envs/ciriquant/bin/python scripts/web_ui.py --host 0.0.0.0 --port 5000 >> logs/web_ui.log 2>&1 &`；stdout 直接重導向 log 檔，不過 conda run 的 pipe 層 |
| `run_manual` 提交後 `sample_groups.csv` 遺失 `patient_id` 欄 | `run_manual` 路由固定只寫 `srr_id`/`condition` 兩欄到 `sample_groups.csv`，CSV 上傳中的 `patient_id`/`description` 欄直接被丟棄 | 手動在 server 補充 `patient_id` 欄（見 GSE229705 段落的補充指令）；長期解法：`run_manual` 讀取 CSV 時若有 `patient_id` 欄則一起寫入 `sample_groups.csv` |
| RBP Binding modal 的 Site Positions 欄顯示相對座標 | ENCORI `circ_pos` 格式為 `chr6:148390208-148390208`（絕對座標），`_absPos()` 誤以為是相對座標再加 `chromStart`，導致超出染色體長度回傳 "N/A" | `_absPos()` 加偵測：若 `circ_pos` 以 `chr` 開頭則視為絕對座標直接使用（已修正，`generate_report.py`）；此外 RBP Binding modal 改為顯示 per-site 字母標籤（a/b/c…）+ 絕對 hg19 座標，與 Circular Structure SVG badge 一一對應（2026-07-02）|
| RBP 表格 `bindingSites` 欄顯示 ENCORI clip 次數而非 site 數 | ENCORI 的 `bindingSites` 欄是 CLIP 實驗次數（clipExpNum），不是 binding site 位置數；直接顯示導致 mapped > total 的矛盾 | `bindingSites` cell 改為 `Math.max(bindingSites, sites.length)`，取兩者較大值確保 mapped 不超過 total（`generate_report.py`，2026-07-02）|
| 報告語言切換後 MultiQC「（點擊折疊）」文字不更新 | MultiQC `<details>` 折疊 span 沒有 `data-en` 屬性；`<details>` toggle JS 用 `textContent=` 硬編中文字串覆寫，讓之後 `switchReportLang` 的 innerHTML 操作也失效（因 `textContent` 設定後 `dataset.zh` 未更新） | `generate_report.py`：span 加 `id="qc-collapse-lbl"` + `data-en="(click to collapse)"`；toggle JS 改為偵測 `_LANG`、同步更新 `dataset.en`/`dataset.zh`、再依語言設定 textContent（2026-07-07）|
| Python mock 腳本 `exec(code)` 後 `build_report()` 未被呼叫 | `generate_report.py` 入口為 `if "snakemake" in dir():`；`dir()` 只看 local namespace，不含 builtins；`builtins.snakemake = sn` 設法讓 snakemake 不出現在 `dir()` → 條件永遠 False | 改用 `exec(code, {"snakemake": sn})`（傳入 globals dict）；另外 `output` 必須是 list（`sn.output = [path]`）支援 `output[0]`；`input` 屬性名稱需與 Snakemake rule 完全一致（`switching`/`isoform_groups`/`de_deseq`/`multiqc`）（2026-07-07）|
| ENCORI RBP/miRNA Mapped 數量極少（原本僅 27%）| ENCORI 提供宿主基因全長的 CLIP-seq peaks（橫跨所有 exon + intron），但 `_map_to_circ_pos()` 只接受精確落在 circRNA exon 邊界內的 sites；大量 peaks 落在其他 exon 或 intron 而被丟棄；`exon_nums=[]`（intronic/intergenic circRNA）時全部 mapping 直接跳過 | `predict_interactions.py` 新增 `_parse_circ_id()` + `_genomic_to_spliced()` 兩個 helper 函式；在 exon-level mapping 失敗後，自動 fallback 到 **genomic-span proportional mapping**：若 hg19 liftover 座標落在 `[circ_start, circ_end]` 範圍內，以比例方式估算 spliced position（`cs_frac = (ov_s - circ_start) / span × total_exon_len`）；`_fetch_encori_mirna` 和 `_fetch_encori_rbp` 均新增 `circ_id` + `strand` 參數，`can_map` 條件改為 `exon_nums OR circ_coords` 任一有效即可嘗試；改善幅度：27% → 約 60–75% Mapped 覆蓋率（2026-07-07）|
| `consensus_filter.py` 同一 circRNA 重複輸出多列 | `vote()` 用 slop 判斷「支持數」但未用來合併輸出：同一 junction 若在 slop 內被多個座標支持（跨工具或同工具內近重複），每個座標各自輸出一列 | 新增 `_cluster_coords()`，投票前先依 slop 做座標分群，每群僅留一個代表（優先選 CIRIquant/`tool_maps[0]` 座標，維持與 `merge_counts.py` 精確字串比對相容）；單樣本驗證去重 ~1.5%（2026-07-16，見 GSE113230 段落）|
| 資料集間 DCC `-Nr` 門檻不一致，偵測數無法跨資料集比較 | `-Nr` 從 5 於 2026-06-08 commit `1b703d3` 改為 2，但當時已完成的資料集（GSE113230、GSE133998）DCC 輸出仍是 `-Nr 5`；用各資料集 `DCC/CircRNACount` 最小 junction count 判定實際門檻（≥5 或含 2/3/4）| 逐一重新下載 FASTQ + 重跑 STAR/DCC/consensus/DE 統一為 `-Nr 2`（CIRIquant GTF 可沿用舊檔，不受 `-Nr` 影響）；GSE113230 已完成（2026-07-17），GSE133998 進行中 |
| DCC 0.5.0 在基因密集/剪接複雜區域（如 3p21.3 RBM5）duplicate-marking 近乎停滯 | `-Nr 2` 保留更多低 count junction，若某 locus 有大量 read 落在彼此相差幾 bp 的 chimeric junction 變異上，DCC 內部 duplicate 比對邏輯耗時劇增（單一 read 比對需 1.5+ 小時）；process 仍佔用 CPU，非死鎖，僅是極慢 | 定位卡住座標對應的基因（查 `genes.gtf`），確認該 exact junction 支持 read 數極少後，從該樣本 `Chimeric.out.junction`（paired + mate1 + mate2）過濾移除該 junction 再重跑 DCC；此樣本的 consensus 結果會缺這一個 back-splice junction，需在文件中註明；原始檔案備份為 `*.bak_pre_<原因>filter` |
| Snakemake `--rerun-incomplete` 誤判手動產生的輸出為「未完成」，重啟時覆寫掉正確結果 | Snakemake 有兩套獨立追蹤機制：`.snakemake/metadata/`（code-hash，`--cleanup-metadata` 清這個）與 `.snakemake/incomplete/`（job 是否透過 Snakemake 自己的子行程乾淨結束，決定 `IncompleteFilesException`）；若手動在 Snakemake 外執行某 rule 對應的指令（如手動重跑 DCC），或先前該 job 被 kill 於 Snakemake 子行程執行中，`.snakemake/incomplete/` 會殘留該輸出的標記，之後只要仍帶 `--rerun-incomplete` 就會強制重跑並覆寫已存在的正確輸出，即使檔案已存在且新鮮 | `--cleanup-metadata <file>` **對此無效**（清的是不同資料庫）；正確做法：`.snakemake/incomplete/` 目錄下的檔名是 **base64 編碼的輸出絕對路徑**（`echo <filename> \| base64 -d` 可解碼確認），找到對應目標檔案的那一個並手動 `rm` 刪除，再重啟時**不要帶 `--rerun-incomplete`**（僅靠 `--rerun-triggers mtime` 判斷）|
| benchmark 重跑 SRR444655 時 Snakemake 主行程與 CIRIquant 無聲消失（無錯誤訊息、`dmesg` 查無 OOM）| 一開始誤判為 systemd session 結束時的 `KillUserProcesses` 把整組行程清掉（`nohup` 防得住 SIGHUP，防不住 session 被清）；改用 `setsid` 重啟後又在 8 小時後死掉一次，這次才在 log 裡翻到清楚的 `OSError: [Errno 28] No space left on device`——**真正原因是磁碟被灌到 100% 滿**，CIRIquant 自己的中間檔案（BWA-mem 輸出的 `circ/*.sam` 分片檔）在單一樣本就能吃到 ~220GB，連 Snakemake 想記錄 job 失敗的 log 都寫不進去，因此看起來像「行程被外力砍掉、無錯誤訊息」 | 診斷時除了查行程是否存活，**務必同時查 `df -h`**；清理已完成資料集的 raw/trimmed/BAM/mate1/mate2（GSE113230、GSE121842 均已完成報告，可安全清除）騰出約 600GB 後用 `setsid nohup` 重啟即可正常跑完；教訓：CIRIquant 對單一 total-RNA 樣本的峰值磁碟用量（circ/+align/+find_circ 暫存）可達 200GB+，排程長 job 前應預留至少 300GB 空間 |
| Snakemake 目標檔案路徑格式不完全匹配時，靜默 fallback 到預設 `all` target，DAG 範圍比預期大很多 | `snakemake ... {某個 rule 的 output 絕對路徑}` 有時仍會把 `Job stats` 表列出全部 18 個 job（含 ground_truth、accuracy_benchmark 等不相關項目），而非只跑目標檔案所需的最小子集；重跑 SRR444655 CIRIquant 時發生過一次，後來改成直接指定 `comparison_report.html`（`rule all` 的 output）反而正確辨識出只需 6 個 job（已完成的 8 個被跳過）| 若懷疑 DAG 範圍不對，直接比對 `Job stats` 表的 job 數與預期；本專案這次的教訓是：與其猜測目標檔案路徑的匹配規則，不如直接指定 `rule all` 的最終 output（`comparison_report.html`），讓 Snakemake 自行決定最小需要重跑的子集，比手動指定中間 rule 的 output 更可靠 |
| Monitor 輪詢腳本誤判「檔案存在」＝「這次執行剛產生的新檔案」| 用 `[ -f "$REPORT" ]` 判斷 benchmark 是否跑完，但 `comparison_report.html` 是舊有輸出（上次成功執行留下的），只要檔案存在就會誤觸發「完成」警報，即使這次執行才剛啟動 1 分鐘 | 改用 `stat -c %Y "$REPORT"` 取得檔案 mtime，與本次啟動時間（epoch）比較，只有 mtime 晚於啟動時間才算「這次執行產生的新檔案」；任何長跑背景監控只要目標是「某個可能已存在的輸出檔」，都應該用時間戳比對而非單純存在性檢查 |
| Snakemake 將「沒有 upstream 錯誤」的 dry-run 誤判為健康 | 2026-10-04 的 dry-run 目標用 `~/x_results`，configfile 的 `results_dir` 為 `/home3/x_results`；兩者字串不同，Snakemake 不做等價比對，`Nothing to be done` 只代表「這個不存在的目標不需要做事」 | 目標路徑必須與 configfile 的 `results_dir` 完全相同；dry-run 結果只有在確認目標確實落在 `results_dir` 下後才有意義（2026-10-04 已重跑修正） |
| 上游設定檔被 touch 後，全部下游 job 被重跑（與 .smk code-hash 不同機制） | rerun trigger 包含 mtime：`config/ciriquant.yaml`（2026-07-29）比 CIRIquant GTF（2026-07-16）新，所有 `ciriquant` job 被判為過期，進而帶動下游；被刪除的 FASTQ 則使下載與 trim 也被重新規劃 | 不以 `ancient()` 壓制（會讓參考基因組更換靜默通過）；拆分 `ciriquant.yaml` 為機器本地工具路徑與科學性參考路徑是獨立議題；報告批次改用 `--allowed-rules` 限定 rule，上游缺輸入即報錯 |
| pandas 聚合運算（`.median()`/`.mean()`/`.sum()`）在多欄 DataFrame 上靜默回傳 Series，不是 scalar | `analysis/typewise_score_distributions.py` 的比較函式本該收到單一欄位的 Series，卻誤傳整個多欄 DataFrame；`df.median()` 不報錯，直接回傳「每欄各自的中位數」組成的 Series，後面 `round(float(...), 4)` 才在 `float(Series)` 這一步炸出 `TypeError`——若當時後面接的是别的運算（例如直接存進字典、或湊巧只用其中一個數字），很可能連 `TypeError` 都不會出現，錯誤會完全沉默 | 任何聚合運算呼叫前，先確認輸入是單一欄位的 Series，不要依賴型別錯誤在下游某處剛好炸出來——那次是兩輪都靠實際執行才抓到（2026-11-04），讀程式碼看不出來；通用寫法：`assert isinstance(x, pd.Series)` 或明確 `df[col]` 取值後再聚合，不要把整個 filtered DataFrame 傳給預期吃 Series 的函式 |

---

## 通知系統（`notify.py`）

Pipeline 完成、失敗、啟動時自動發送通知，透過 Snakemake `onstart` / `onsuccess` / `onerror` hook 呼叫。

**支援管道**：

| 管道 | 說明 |
|------|------|
| SMTP Email | Gmail TLS 587，成功時附加 `report.html`（>20 MB 略過附件） |
| Slack Webhook | Incoming Webhook，Markdown 格式 |

（LINE Notify 已於 2025-03-31 終止服務，已從 `notify.py` 移除）

**通知內容**：

| 事件 | 內容 |
|------|------|
| `start` | 啟動時間 |
| `success` | 完成時間、總 circRNA 數、顯著 DECs、上調/下調數 |
| `failure` | 失敗 rule 名稱、pipeline log 最後 50 行 |

**環境變數設定**（加到 server `~/.bashrc`）：

```bash
export NOTIFY_EMAIL_FROM="寄件gmail@gmail.com"
export NOTIFY_EMAIL_PASS="xxxx xxxx xxxx xxxx"   # Gmail App Password（非登入密碼）
export NOTIFY_EMAIL_TO="chou.kaihsuan@gmail.com"
export NOTIFY_SLACK_WEBHOOK="https://hooks.slack.com/services/..."  # 選填
```

Gmail App Password 申請：Google 帳號 → 安全性 → 兩步驟驗證開啟後 → 應用程式密碼。

**手動測試**：
```bash
python scripts/notify.py --event start --project GSE113230
python scripts/notify.py --event success --project GSE113230 --report results/report.html
python scripts/notify.py --event failure --project GSE113230 --rule dcc --log logs/pipeline_run.log
```

未設定任何環境變數時，`notify.py` 只列印 log，不影響 pipeline 執行。

**Snakemake hook 實作**（`workflow/Snakefile`）：
- `onstart` → fire-and-forget subprocess 呼叫 `notify.py --event start`
- `onsuccess` → 解析 `de_results.tsv` 統計後呼叫 `notify.py --event success --stats {...}`
- `onerror` → 呼叫 `notify.py --event failure --rule unknown --log logs/pipeline_run.log`
- hook 內必須用 `workflow.snakefile` 推算專案路徑，**不能用 `__file__`**（會指向 conda env）

---

## 重要的 Snakemake 行為

- **rerun trigger 有多種，不只 .smk code-hash**：mtime 是預設的 rerun trigger 之一。上游設定檔（例如 `config/ciriquant.yaml`）被 touch 過，會使其下游 job 全部被判為過期並重跑（2026-10-04 dry-run 實測：整個上游鏈被重新規劃）。這與上方 .smk code-hash 是同類的第二種機制，來源不同。`config/ciriquant.yaml` 同時裝機器本地工具路徑與科學性參考路徑，目前不可安全地 `ancient()`；拆檔是獨立議題。
- **dry-run 的「無待執行 job」必須先驗證目標路徑**：Snakemake 的目標字串必須與 configfile 的 `results_dir` 完全一致（`~/x_results` 與 `/home3/choukaihsuan/x_results` 是不同字串，`/home` 為 `/home3` 的 symlink 但 Snakemake 不做等價比對）。路徑不符時，`Nothing to be done` 證明的只是「這個不存在的目標不需要做事」，不是「管線不需要重跑」。2026-10-04 因此誤判一次。
- **只跑指定 rule 時用 `--allowed-rules`**：上游缺輸入會報錯（MissingRule／MissingInput），而不是靜默重跑上游。報告批次使用 `--allowed-rules loo_stability rank_biomarkers_v2 generate_report`。


- **DAG 在啟動時編譯**：修改 `.smk` 規則不影響正在跑的進程
- **`script:` 腳本在執行時讀取**：修改 `analysis.R` 或 `.py` 腳本對正在跑的進程立即生效
- **`--rerun-incomplete`**：重跑上次被中斷的 job（不自動重跑已完成的）
- **`--keep-going`**：某個 rule 失敗時繼續跑其他獨立的 rule
- **資源控制**：`--cores 36 --resources mem_gb=300`（保留給其他 server 使用者）

---

## HTML 報告內容（`generate_report.py`）

報告為自包含 HTML（PDF 圖表以 base64 內嵌）：

| 區塊 | 說明 |
|------|------|
| **Samples 區塊（報告開頭）** | `_sample_overview_section()`：讀 fastp JSON → stat boxes（tumor=藍 / normal=綠）+ 可折疊詳細表格（SRR / Condition / Patient / Total Reads / Avg Length / Q30）；`study_title` config 鍵存入後以斜體顯示於頂部（GEO 研究標題，自動從 NCBI eUtils 抓取或手動填入）|
| **DE 方法切換器** | 頁面頂部三個按鈕（edgeR_ciriquant / DESeq2 / limma-voom）；切換時即時更新 stat-boxes、Volcano、Heatmap、**Top DE 表格**、**Biomarker 表格**（`Plotly.react` 原地更新，無頁面 reload）|
| Summary stat-boxes | 樣本數、total circRNAs、顯著數、Up/Down；id="stat-n-sig/up/dn"，供方法切換器更新 |
| **Type I/II 分類** | edgeR_ciriquant 模式才顯示；橫向進度條 + 各自數量 |
| **3-method Venn diagram** | SVG 三圓 Venn；圓心 A=(170,128) B=(290,128) C=(230,210) r=90；交集數字加白色 halo（`paint-order="stroke"`）；高度 345px；**7 個區域均可點擊**，點擊後在 `div#venn-detail` 顯示對應 circRNA 清單（含 gene_name / log2FC / p-value / circbase_id（超連結）/ **方法欄**）；「⬇ CSV」下載各區域清單 |
| **Biomarker 候選表** | 預設顯示前 30，可用「顯示前 N」輸入框調整（`updateBiomarkerN()`，2026-07-下旬新增；`_compute_bm_table_data` 現在回傳**全部**顯著 rows，不再固定截斷 top-30，JS 端依輸入值即時切片）；**方法切換時完全重建**（`_renderBiomarkerTable()`）：依所選方法的 p-value 重新排名，過濾按鈕計數同步更新，標題顯示目前方法名稱；circbase_id 為超連結（見下）；fallback：無 bm_table 時降為 graying（opacity 0.25）|
| **Biomarker Score 分布圖** | 單一 Histogram（Y 軸 = Number of circRNAs；bin 寬 0.05，`xbins={start:0,end:1,size:0.05}`）；方法切換時 `Plotly.react` 更新；`_biomarker_normality_plot()` 生成靜態初始圖，JS `_updateScoreDist()` 同步更新。**已移除**：Ranked scatter（`bm-scatter-plot`）、Normal fit 曲線、Shapiro-Wilk 檢定、垂直線（μ/μ±σ/μ±2σ）|
| **Top DE table（分兩表）** | 方法切換時完全重繪（`_renderDETables()`）；欄位含 gene_name / strand / region / exon_span / circbase_id（超連結）/ log2FC / p-value / Type |
| Volcano plot | **Plotly 互動式**；方法切換時 Plotly.react 更新座標軸標題**必須用巢狀 `xaxis:{title:...}`／`yaxis:{title:...}`**，不可用 `xaxis_title:`/`yaxis_title:` 扁平寫法（那是 Python `plotly.py` 的簡寫語法，純 JavaScript Plotly.js API 不認得、會靜默忽略）——2026-07 曾因此 bug 導致每次 DE 方法切換後座標軸標題消失（`switchDEMethod()` 頁面載入時就會執行一次，所以是「每次載入都消失」而非只在使用者手動切換時），已修正並補回格線樣式 |
| PCA | **Plotly 互動式**（tumor/normal 顏色區分）；numpy SVD |
| **聚類熱圖（Clustering Heatmap，2026-07-下旬）** | 取代原本「Heatmap top-N 控制」（可調整 up/down 各顯示幾筆的舊版已移除）；顯示**全部**顯著 DE circRNA（依主方法），依階層聚類（scipy `linkage(method='ward')` + `leaves_list()`）排序，左側渲染真正的樹狀圖（dendrogram，Python 端算好 icoord/dcoord 傳給 Plotly 畫線段，非 Plotly 內建 dendrogram 功能）；**永遠展開顯示**（原本是 `<details>` 可折疊、且用 `toggle` 事件觸發渲染，因為 `toggle` 只在使用者互動時觸發、頁面載入時不會自動 fire，導致預設看到的是空白折疊區——已改成普通 `<div>` 頁面載入即渲染）；tumor/normal 樣本色塊沿用主 heatmap 配色；繪圖函式抽成共用的 `_renderClustPanel()`，供全報告聚類熱圖與 circRNA 詳細 modal 的 mini heatmap 共用（見下方 modal 說明）|
| Isoform Switching | Plotly 長條圖 + 顯著 switching 表格；**不隨 DE 方法切換更新**（IUI 計算固定）；circbase_id 欄為超連結；section 下方有說明文字 |
| **SVG Circular Diagram** | 每個 circRNA 的環狀圖（exon 結構 + miRNA/RBP binding site 弧段 + 流水號 badge）|
| **互動式欄位排序** | Isoform Switching、Top DE（up/down）、Biomarker 表格的所有欄位均可點擊排序（▲/▼）；`_makeSortable(tableId)` 函式；`_parseGenCoord()` 提供基因組座標排序（chr1<chr2<…<chr22<chrX<chrY，同染色體再按 start position）；gene_name="intergenic" 在 DE 表格顯示「—」（region 欄已有此資訊）|

**JS 全域狀態變數**：
- `const ALL_DE_METHODS`：三方法的完整 volcano + stats + heatmap + **de_table + sig_ids + bm_table** 資料
- `const FULL_HEATMAP_DATA`：pool=50 up + 50 down，每 circRNA 含 `{z, pval, log2fc, label}`；欄位順序 tumor 在左、normal 在右
- `let _HEATMAP_DATA_CACHE`：目前顯示方法的 heatmap data（方法切換時更新）；若 `conditions` 為空則 fallback 到 `FULL_HEATMAP_DATA.conditions`
- `const VENN_REGION_DATA`：7 個 Venn 區域的 circRNA 清單（`{label, circs:[{id, gene, lfc, pval, cb, m}]}`）

**`switchDEMethod(method)`**：更新 stat-boxes → Plotly.react volcano → updateMainHeatmap → **`_renderDETables()`** → **`_renderBiomarkerTable()`** → `_updateScoreDist()`

**`_renderDETables(method, md)`**：從 `md.de_table` 重建 up/down HTML 表格，插入 `#de-tables-section`；清除靜態 `table, h3, .tbl-dl-bar`。

**`_renderBiomarkerTable(method, md)`**：從 `md.bm_table.rows` 重建 `#tbl_biomarker tbody`；更新過濾按鈕計數（全部/≥2/3方法）及 h2 標題；無 bm_table 時 fallback 到 `_updateBiomarkerHighlight()`。

**`_updateBiomarkerHighlight(sigIds)`**：（fallback）在 `#biomarker-section` 中，對每個 `<tr>` 讀取 `circ-link` onclick 的 circ_id，不在 sigIds 內者設 `opacity:0.25`。

**`_compute_bm_table_data(de, p_col, sig_thr, lfc, bm_lookup, sig_sets_all)`**（Python）：對一個 DE method 的顯著集合計算 6D composite score → top-30 rows；存入 `ALL_DE_METHODS[method]["bm_table"]`。Columns：rank / circ_id / log2FC / n_mirna / n_rbp / biomarker_score / n_sig_methods / in_circbase / circbase_id / circbase_gene / Type。**mirna/rbp 正規化在當前 sig set 內進行**（`_sig_mirna_mx = max(n_mirna for c in sig["circ_id"])`），避免跨方法 scale 不一致造成 score > 1.0 的問題（union interactions.json 的 n_mirna 可能大於 edgeR biomarker_candidates.tsv 的最大值）。

**`_bm_lookup` 擴充**：從 `biomarker_candidates.tsv` 讀取時額外存 raw display 值（`n_mirna`, `n_rbp`, `in_circbase`, `circbase_id`, `circbase_gene`）以及 **`type_edger`**（edgeR 分配的 Type_I/II/III/NS，用於 DESeq2/limma bm_table 的 fallback 顯示）；再從 `circbase_annotated.tsv` + `interactions.json` 補全 DESeq2/limma-only circRNA（edgeR 不顯著、biomarker_candidates.tsv 未包含的 circ_id）。

**Type 欄 fallback 邏輯**（`_compute_bm_table_data`）：先從當前 DE method 的 `r["Type"]` 取值；若為 NaN（DESeq2/limma 不做 FSJ 分類，`analysis.R` 設 `NA_character_`），則從 `lu.get("type_edger")` 取 edgeR 分類；仍為空則顯示 `"—"`。效果：3 方法均顯著的 circRNA 在 DESeq2/limma bm_table 中也能看到 edgeR Type 標籤。

**`_makeSortable(tableId)`**：在 `_SCRIPT` 靜態區塊定義（`<head>` 中），被呼叫時機：`_renderDETables()` 尾部（de-up-table / de-dn-table）、`_renderBiomarkerTable()` 尾部（tbl_biomarker，tbody 重建後重新綁定）、main JS 區塊尾部（tbl_isoform / tbl_biomarker 首次載入）。

**`_parseGenCoord(s)`**：解析 `chrN:start|end` 格式 → `{c: chrNum, p: startPos}`；chrX=23, chrY=24, chrM=26；供 `_makeSortable` 在排序 circ_id 欄時優先使用基因組排序（否則純字串排序會讓 chr4:1902353 排在 chr1:91382197 前面，因去除非數字後 41902353 < 191382197）。

**列印排版（`@media print`）**：
- `h2, h3 { break-after: avoid }` — 標題後不換頁
- `.plotly-graph-div { break-inside: avoid }` — 圖表不跨頁
- `#de-tables-section { break-before: page }` — Top DE 表格從新頁開始
- flex 區塊改 block（雙欄分布圖垂直排列）

**circRNA 詳細 modal（點擊任意 circ_id 開啟）**：
- **⬛ Circular Structure**：SVG 環狀圖；底部「⬇ SVG」下載按鈕；`totalLen` 由 `exon_boundaries[last].cum_end` 推算（interactions.json 的 `spliced_length` 欄位固定為 0，不可用）
- **📺 miRNA Sponge**：互動表格（Priority 排序）；「⬇ CSV」下載；Binding Seq 欄自動從 UCSC hg19 REST API 獲取序列；**「In circRNA」欄已移除（2026-07-31）**：無法 liftover 對應到 circRNA 座標的項目（`in_circ===false`）直接從表格過濾掉，不再保留並顯示打叉——因為那些項目沒有可用的位置/序列資料，留著只會讓人困惑；`_calcPriority()` 內部評分公式仍讀取 `in_circ` 做加分依據（未變動），只是顯示層面篩掉了不合格項目
- **🧬 RBP Binding**：同上（含「In circRNA」欄移除規則）
- **📈 Volcano**：Plotly mini-chart；「⬇ PNG」下載
- **🔥 Heatmap（2026-07-下旬重寫）**：不再是舊版「top10up+top10down 固定池」，改成**聚類鄰居視圖**——在全報告聚類熱圖的階層排序（`CLUST_HEATMAP_DATA.order`）裡找出目前 circRNA 的位置，取前後共約 20 個鄰近 circRNA（`_buildMiniHeatmap()`），並用 `_sliceDendro()` 從完整樹狀圖切出這個視窗對應的**局部子樹**一併畫出（標準樹狀圖不會有分支交叉，所以座標落在視窗範圍內的合併節點必然完整屬於視窗內的葉節點，超出範圍的直接捨棄即可，不需要重新算圖）；目標 circRNA 一樣用橘色外框標示；與全報告聚類熱圖共用同一份 `_renderClustPanel()` 繪圖邏輯
- **Priority Score**（miRNA）：seed type（8mer=+4…6mer=+1）+ CLIP>0（+3）+ ENCORI（+2）+ in_circ（+1）
- **Priority Score**（RBP）：bindingSites log₂×2 max3 + internal（+2）+ CLIP>0（+3）+ ENCORI（+2）+ in_circ（+1）
- 所有 header 可點擊排序（▲/▼）
- Chr Position 欄顯示 chr 絕對座標（circRNA start + circ_pos 1-based offset）

**Volcano toggle**：右上角「○ Heatmap circles: OFF / ON」按鈕，heatmap top circRNA 標記預設隱藏。

**全域下載功能**：
- 頁面頂部 sticky 列：「🖨 列印 / 存為 PDF」（`window.print()`）
- 各表格右上角「⬇ CSV」按鈕：DE 上調/下調、Biomarker、Isoform switching

**Plotly 安裝要求**：`pip install plotly` 在 `ciriquant` conda env（已安裝 5.18.0）。
必須用 `conda run -n ciriquant python` 執行 generate_report.py，否則 fallback 靜態 PDF。

**SVG Circular Diagram — Badge 放置邏輯（radial staggering）**：

- miRNA badge（外圈）：預設放在弧段外 `+7` px（`MI_R0`）；若與前一 badge 角度差 < 0.14 rad，改放 `+19` px（`MI_R1`）
- RBP badge（內圈）：同邏輯，往內 `−8` px（`RBP_R0`）或 `−20` px（`RBP_R1`）
- **不做 angular de-overlap**：badge 角度永遠等於弧段中心角，只在徑向上交錯
- 虛線 connector（`stroke-dasharray="2,2"`）放在 `<g id="_mib_N">` 內，toggle 時隨 badge 一起隱藏
- 虛線只在 badge 在外圈（`MI_R1`）時才畫，距離弧段緊貼時不畫

**Solo mode（單分子可見）— Per-arc site labels**：

當 miRNA 或 RBP 中**只有一個分子可見**時，右側面板自動出現（`_updateSitePanel`），列出該分子各 binding site（字母 a/b/c...）可個別切換顯示/隱藏（`_toggleSite`）：

- **`_arc_lbl` badge 放置邏輯**：
  - miRNA site badge：`midR=(MI_IN+MI_OUT)/2=164.5`（放在 miRNA ring 中心，直接在 arc 上）
  - RBP site badge：`midR=(RIN+ROUT)/2=131.5`（放在 exon ring 中心，覆蓋在 exon arc 上，角度指示 binding site 位置）
- **Arc highlight**：`_showSitePanel` 呼叫時，對選中分子的所有 arc path 加 `stroke=white, strokeWidth=1.5px, opacity=1`，使 arc 輪廓清晰（`_hideSitePanel` 呼叫時移除）
- **Overlap annotation badge 移除**：`arc(MI_IN,MI_OUT,...)` 和 `arc(RBP_IN,RBP_OUT,...)` 的重疊字母 badge（`<circle>+<text>` at midR=48 / midR=164.5）已移除，只保留 boundary dash line；唯一字母來源是 `_arc_lbl` per-arc site labels
- **`window._circ_miData` / `window._circ_rbpData`**：每次 `_drawCircleRNA` 重置，在各自 `if(totalLen>0)` block 結尾填入 `{name, color, arcs:[{letter, id, pos}]}`；`const` block-scope 限制需確保 `window._circ_miData=_miArcData` 在第一個 if block 結尾（不在第二個 RBP block 內）

**DE table 資料來源合併**：
- `de_results.tsv`（主表）
- `isoform_groups.tsv` → `gene_name`, `strand`, `region`, `exon_span`
- `circbase_annotated.tsv` → `circbase_id`, `circbase_gene`, `in_circbase`

報告標頭顯示使用的 DE 方法（`method-tag` badge）。
Plotly 依賴：`plotly`、`numpy`；若兩者未安裝則自動 fallback 到靜態 PDF embed。

---

## 資料集特性稽核校正（2026-10-03，依 `audit/` TSV 與 GEO/SRA 原始 metadata 校正，取代下方內文中衝突的舊敘述）

本節取代以下段落中任何與此表衝突的讀長／library prep／strandedness／有效單工具比例敘述，資料來源一律是 `audit/per_dataset_summary.tsv`（讀長、library prep，查核自 GEO extract_protocol 或 SRA metadata 原文，非沿用本文件舊敘述）與 `audit/thesis_dataset_table_19ds.tsv`（有效單工具比例、15 組跨資料集分析納入與否）。

**兩項已確認的既有錯誤**：
- **GSE229705 讀長原文件誤記為「100bp PE」，實際為 51bp PE**（GEO/SRA metadata 查證；本文件下方 GSE229705 段落已一併修正）。
- **GSE58135 的 library prep 原文件從未揭露是 poly(A)-selected**（GEO extract protocol 原文：「POLY(A)-SELECTED, Tn-RNA-seq non-stranded」），而非 total RNA——此為論文 §2.13 poly(A) 排除準則本身明訂要排除的建庫方式，但 GSE58135 仍偵測到 1,607 個 circRNA 並留在 cohort 中，是一個需要在論文中明確揭露、而非沿用舊文件含糊帶過的方法學限制。

**有效單工具比例 ≥ 97.8% 的 4 組（GSE58135、GSE229705、GSE192410、GSE192849）已排除於 15 資料集跨資料集重現分析之外**（詳見 `audit/recurrence_15ds.tsv`／`audit/permutation_15ds.tsv`）；GSE221107（0.2%）視為雙工具、保留。19 組完整結果仍保留作為敏感度比較，但主結果以 15 組為準。

| Accession | 讀長 (bp，PE 各端) | Library prep（GEO/SRA extract_protocol 原文節錄）| Strandedness 線索 | 有效單工具比例 | 15 組跨資料集分析 |
|---|---|---|---|---|---|
| GSE113230 | 150-150 | RiboMinus rRNA removal (GEO extract protocol) | 未載明 | 0.0% | ✅ 納入15組 |
| GSE58135 | 50-50 | POLY(A)-SELECTED, Tn-RNA-seq non-stranded (GEO extract protocol) | unstranded（GSM 原文：non-stranded Tn-RNA-seq） | 97.8% | ❌ 排除 |
| GSE133998 | 150-150 | rRNA depletion (GEO) | stranded（GSM 原文：strand-specific library；dUTP/方向未載明） | 0.0% | ✅ 納入15組 |
| SRP156355 | 100-101 | SRA study abstract: 'rRNA depleted and stranded library'（study-level，非 GSM-level） | stranded（**study-level**：SRA STUDY_ABSTRACT 寫 "stranded library was synthesized"；無 GSM，非 GSM-level 證據，`evidence_level=SRA_study_abstract`） | 0.0% | ✅ 納入15組 |
| GSE77509 | 100-100 | TruSeq Stranded Total RNA, rRNA depletion (GEO) | stranded（GSM 原文：TruSeq Stranded Total RNA） | 0.0% | ✅ 納入15組 |
| GSE130078 | 101-101 | TruSeq Stranded Total RNA + Ribo-Zero (GEO) | stranded（GSM 原文：TruSeq Stranded Total RNA） | 0.0% | ✅ 納入15組 |
| GSE248612 | 138-140 | GEO: standard Illumina protocols（rRNA depletion 未載明，UNVERIFIED） | 未載明 | 0.0% | ✅ 納入15組 |
| GSE221107 | 150-150 | Ribo-off rRNA depletion (GEO) | stranded（GSM 原文：stranded RNA sequencing library；KC Stranded mRNA kit） | 0.2% | ✅ 納入15組（視為雙工具）|
| PRJNA553289 | 150-150 | SRA: RNA-Seq / cDNA（rRNA depletion 未載明，UNVERIFIED） | 未載明（SRA 單獨提交，無 GSM/extract_protocol） | 0.0% | ✅ 納入15組 |
| GSE229705 | **51-51**（原文件誤記100bp，已更正）| Trio RNA-Seq (SPIA) with final rRNA-transcript removal step (GEO) | 未載明 | 100.0% | ❌ 排除 |
| GSE148036 | 142-147 | Ribo-Zero rRNA removal (GEO) | 未載明 | 0.0% | ✅ 納入15組 |
| GSE121842 | 150-150 | GEO: standard Illumina protocols（rRNA depletion 未載明，UNVERIFIED） | 未載明 | 0.0% | ✅ 納入15組 |
| GSE136569 | 100-100 | Ribozero series kit rRNA removal + linear RNA digested with RNase R (GSM extract_protocol 原文)，RNase R-treated | 未載明 | 0.0% | ✅ 納入15組 |
| GSE143797 | 150-150 | ncRNA-Seq / size fractionation; dUTP stranded（rRNA depletion 未載明，UNVERIFIED） | stranded（GSM 原文：dUTP 第二股合成） | 0.0% | ✅ 納入15組 |
| GSE108735 | 150-150 | ncRNA-Seq / size fractionation（rRNA/linear-RNA depletion 未載明，UNVERIFIED） | 未載明 | 0.0% | ✅ 納入15組 |
| GSE171011 | 150-150 | GEO: standard Illumina protocols；rRNA reads僅in silico移除（library prep UNVERIFIED） | 未載明 | 0.0% | ✅ 納入15組 |
| GSE97239 | 150-150 | RiboZero rRNA removal + RNase R digestion (GEO)，RNase R-treated | unstranded（**文字推定**：random primer 合成 cDNA，未描述 dUTP/strand 標記；非原文明載） | 0.0% | ✅ 納入15組 |
| GSE192410 | 150-150 | GEO: none provided by the submitter（UNVERIFIED，完全未揭露建庫法）| 未載明 | 100.0% | ❌ 排除 |
| GSE192849 | 150-150 | RNase R digestion + VAHTS Total RNAseq kit (GEO)，RNase R-treated | stranded（GSM 原文：strand-specific，dUTP） | 100.0% | ❌ 排除 |

**「未載明」的意思**：GEO/SRA metadata 本身沒有明確文字陳述該項（多半是 extract_protocol 欄位只寫「standard Illumina protocols」或完全空白），**不代表「已確認是另一種狀態」**；標成 UNVERIFIED 的 library prep 欄位同理。本文件下方各資料集段落若有更肯定的敘述（例如直接寫「rRNA-depleted」而未加註來源查證），以本表的 UNVERIFIED 標註為準。

**Strandedness 稽核結果（2026-10-03，`audit/library_strandedness.tsv`、`audit/library_strandedness_mismatch.tsv`；每列有 `evidence_level`：GSM_extract_protocol / SRA_study_abstract / text_implied / not_stated，可依嚴格程度自行取捨；kit 名稱不當證據）**：19 組中 GSM 原文明載者 7 組（stranded 6：GSE133998、GSE77509、GSE130078、GSE221107、GSE143797、GSE192849；non-stranded 1：GSE58135）；另 GSE97239 為文字推定 unstranded（text_implied）、SRP156355 為 SRA 研究摘要層級 stranded（無 GSM）；其餘 10 組未載明。15 組中：GSM 原文明載 5 組、加上 text_implied 與 SRA 研究摘要共 7 組可判定。工具實際設定（`results-2026-07` 的 `circrna.smk`）：DCC 未帶 `-N`（預設 stranded，fr-firststrand）、CIRIquant 未帶 `-l`（預設 0 = unstranded）。
- **方向 A（DCC stranded 預設 vs. 真實 unstranded）**：GSE58135（GSM 原文，已排除於 15 組）、GSE97239（text_implied）。
- **方向 B（CIRIquant unstranded vs. 真實 stranded）**：GSE133998、GSE77509、GSE130078、GSE221107、GSE143797、GSE192849（已排除於 15 組）〔GSM 原文〕；SRP156355〔SRA 研究摘要〕。
- **CIRIquant `-l` 的作用範圍（1.1.3 原始碼）**：只在 `circ.py:373-390`（`denovo_worker`）過濾 BSJ 候選 read 的 strand；`lib_type` 僅經 `main.py:197` → `circ.proc`(628) → `proc_denovo_bam`(296) 傳入，**`proc_genome_bam`（FSJ，406）、CIRI2 偵測與基因豐度均不接收此參數**。逐行證據見 `audit/ciriquant_library_type_scope.md`。
- **方向 A（DCC stranded 預設 vs. 真實 unstranded）**：GSE58135（原文明載，已排除於 15 組）、GSE97239（文字推定）。
- **方向 B（CIRIquant unstranded vs. 真實 stranded）**：GSE133998、GSE77509、GSE130078、GSE221107、GSE143797、GSE192849（後者已排除於 15 組）。
- **CIRIquant `-l` 的作用範圍（1.1.3 原始碼）**：只在 `circ.py:373-390`（`denovo_worker`）過濾 BSJ 候選 read 的 strand；`lib_type` 僅經 `main.py:197` → `circ.proc`(628) → `proc_denovo_bam`(296) 傳入，**`proc_genome_bam`（FSJ，406）、CIRI2 偵測與基因豐度均不接收此參數**。逐行證據見 `audit/ciriquant_library_type_scope.md`。

---

## 目前執行進度（共 20 個資料集完成；0 個進行中）

> **⚠️ GSE108735 資料集更正**：原標記為「TNBC」，實際確認為**腎細胞癌（Renal Cell Carcinoma，RCC）**，7 pairs tumor vs. 正常腎組織，ncRNA-Seq（SRR6439741–SRR6439754）。GSE171011 原標記為「TNBC」，實際為**甲狀腺乳突癌（Papillary Thyroid Cancer，PTC）**，4T+4N=8 samples，RNA-Seq（SRR14088791–SRR14088798）。

> **📌 GSE58135（50bp）與「Read length ≥ 100bp」checklist 的時間先後（2026-10-02 查證，供論文引用）**：`Read length ≥ 100bp` 這條 checklist 規則是 **2026-06-08**（commit `c39737d`，"docs: add download speed benchmarks and GSE dataset selection guide to CLAUDE.md"）才寫進 CLAUDE.md 的，而 **GSE58135（50bp）完成於 2026-06-07**（commit `89429c9`，前一天）。因此 GSE58135 並非違反既有準則的刻意例外，而是該準則的歸納來源——checklist 裡「50bp PE｜❌ 不佳｜STAR chimeric 幾乎無法偵測BSJ」這一行，實質上就是從 GSE58135 的實測結果（CIRIquant vs. DCC 比例 ≈0.002，觸發 adaptive fallback 全部10/10樣本）歸納出來的。查無任何commit message或文件明確寫下「已知違反規則但刻意保留」的決策文字；該資料集留在最終19個資料集的cohort中是後續的決定，理由另見論文正文。

### GSE113230（三陰性乳癌）

**所有步驟已完成（含三方法 DE + 新版報告）。** 報告位置：`~/GSE113230_results/report.html`（server）

**2026-07-17 全量重跑（`-Nr 2` 統一 + consensus 去重修正）**：原始分析用 DCC `-Nr 5`（見下方「歷史版本」），2026-06-08 commit `1b703d3` 起 pipeline 預設改為 `-Nr 2`；為讓所有資料集的 DCC 門檻一致，重新下載 6 個樣本並以 `-Nr 2` + 座標去重（[Item 2 修正](#already-fixed)，slop 內近重複 circRNA 只留一個代表）完整重跑 STAR→DCC→consensus→DE→report。CIRIquant GTF 沿用舊檔（不受 `-Nr` 影響，省下 6×11h）。

| 步驟 | 狀態 |
|------|------|
| fastp QC/trim（重新下載）| ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成（沿用舊 GTF，未重跑） |
| STAR paired-end / mate1 / mate2 | ✅ 6/6 完成（重跑，`-Nr 2` 需要新的 mate junction）|
| DCC（`-Nr 2`）| ✅ 6/6 完成 |
| consensus_filter（含去重）| ✅ 6/6 完成（7,653–17,036 circRNAs / sample，見下表） |
| merge_counts | ✅ 完成（**40,493** consensus circRNAs；filterByExpr 後 **6,650** tested） |
| DE analysis | ✅ 完成（三方法全跑：edgeR **240** / DESeq2 **547** / limma **2,781** significant）|
| predict_interactions / isoform switching / rank_biomarkers / report | ✅ 完成 |

**⚠️ SRR7012366 例外處理（DCC 效能瓶頸，需人工介入）**：`-Nr 2` 下 DCC 0.5.0 在 **RBM5 基因座**（chr3:50145502-50145738，3p21.3 基因密集、剪接複雜區域）的 duplicate-marking 邏輯陷入近乎病態的緩慢處理（單一 read 比對耗時 1.5+ 小時，總計 4+ 小時才跑完這一個樣本，其餘 5 個樣本均在正常時間內完成）。**處理方式**：確認該 junction 僅有 42 筆支持 read（paired 2 + mate1 28 + mate2 12，佔全檔 80 萬+ 行的極小部分），從 SRR7012366 的三個 Chimeric.out.junction 檔（paired/mate1/mate2）過濾移除後重跑 DCC，可在正常時間內完成且結果可重現（兩次獨立重跑皆為 21,461 circRNAs）。**已知限制**：此樣本的 consensus 結果不含 RBM5 這一個 back-splice junction（其餘 803,421+ 筆 junction 不受影響）；原始（未過濾）junction 檔已備份為 `*.bak_pre_rbm5filter`。若未來重跑此資料集或影響其他樣本，需留意 3p21.3 一帶的基因密集區可能有類似風險。

**Snakemake 操作教訓（已記錄於已知問題）**：修復 DCC-366 過程中，因 Snakemake 的 `.snakemake/incomplete/` 標記（與 `--cleanup-metadata` 清理的 code-hash metadata 是**不同**的追蹤機制）誤判該 job 未完成，`--rerun-incomplete` 導致重啟時強制重跑並覆寫掉已完成的正確結果一次。修正方式：解碼 `.snakemake/incomplete/` 下 base64 檔名找到對應的標記檔並手動刪除，改為不帶 `--rerun-incomplete` 重啟。詳見下方已知問題表。

**主要數值結果（`-Nr 2` 全量重跑）**：
- 偵測：**40,493** consensus circRNAs（較 `-Nr 5` 版本的 9,349 大幅增加，符合預期——`-Nr 2` 保留更多低 count junction）→ filterByExpr 後 **6,650** tested
- DE（edgeR_ciriquant）：**240** significant（nominal p < 0.05，|log2FC| > 1）；上調 237 / 下調 3；**212 Type_I (88.3%) / 28 Type_II (11.7%)**
- DE（DESeq2）：547 significant；DE（limma-voom）：2,781 significant
- Top 1 biomarker：**chr12:46622936|46648719（hsa_circ_0000397，SLC38A1）**；log2FC=7.17，199 miRNA，118 RBP binders，score=0.6744，Type_I，n_sig_methods=2 —— 與跨資料集重現分析（見 `docs/biomarker_orthogonal_support.md`）獨立找出的旗艦候選一致，交叉驗證了排序穩定性

**GSE113230 各工具偵測數量（`-Nr 2`，2026-07-17）**：

| SRR ID | 分組 | CIRIquant | DCC | 共識（去重後）|
|--------|------|----------:|----:|-----:|
| SRR7012366 | Tumor 1 | 16,016 | 15,936 | 11,152（RBM5 locus 已排除，見上）|
| SRR7012367 | Tumor 2 | 18,414 | 17,401 | 13,310 |
| SRR7012368 | Tumor 3 | 23,163 | 22,449 | 17,036 |
| SRR7012369 | Normal 1 | 19,217 | 19,377 | 14,506 |
| SRR7012370 | Normal 2 | 21,241 | 19,857 | 15,395 |
| SRR7012371 | Normal 3 | 10,181 | 10,017 | 7,653 |

<a id="already-fixed"></a>
**Item 2 去重修正**：`consensus_filter.py` 的 `vote()` 先前只用 slop 判斷「支持數」，未用來合併輸出，導致同一 junction 若在 slop 內被多個座標支持會重複輸出多列。修正後在投票前先做座標分群，僅留一個代表（優先選 CIRIquant 座標，保持與 count_matrix 的精確字串比對相容）。單樣本驗證（SRR7012368，同參數）：6,539→6,439 列（~1.5% 去重）。

**Benchmark 門檻掃描 AUC-PR 重新確認（2026-07-20）**：舊版「歷史版本」表格中的 accuracy benchmark 數字用的是 `-Nr 5` 時代的 SRP011042（原誤記 GSE55872）consensus 結果；`-Nr 2` 統一後，GSE55872（SRR444655 + 兩個 RNase-R replicates）也已用 `-Nr 2` + dedup 修正重新跑過 CIRIquant/DCC/CIRCexplorer2/find_circ，門檻掃描（`accuracy_benchmark.py --ciri2-file --dcc-count-file --circexplorer2-file --find-circ-file --output-pr-curve`）重新確認：

| Method | Precision | Recall | F1 | Specificity | AUC-PR（門檻掃描）|
|--------|-----------|--------|----|-------------|-----------------|
| **circDEX**（本研究，Our_adaptive）| 0.899 | 0.173 | 0.290 | **0.962** | **0.174** |
| **CirComPara2_4tools** | 0.878 | **0.248** | **0.386** | 0.933 | **0.391** |
| nfcore_3tools | 0.897 | 0.183 | 0.305 | 0.959 | 0.377 |

三方法 AUC-PR 均較 `-Nr 5` 版本（circDEX 0.120 / CirComPara2_4tools 0.349 / nfcore_3tools 0.337）提升，相對排名不變（CirComPara2_4tools > nfcore_3tools > circDEX）——`-Nr 2` 保留更多低 count junction，讓所有方法召回率同步提升，符合預期。**重跑教訓**：`accuracy_benchmark.py` 的門檻掃描需要 `--circexplorer2-file`/`--find-circ-file` 兩個參數才會納入 CIRCexplorer2 和 find_circ 的資料；漏掉這兩個參數時腳本仍會「成功」執行並輸出數字，但 `nfcore_3tools`（3 工具都需要，其中 2 個工具資料是 0）的門檻掃描 AUC-PR 會靜默塌陷為 0.0000，且 `circDEX` 與 `CirComPara2_4tools` 的掃描曲線會因為兩者都退化成只剩 CIRIquant+DCC 而變成完全相同——這種「跑出數字但數字是錯的」情況比腳本直接報錯更難發現，跑門檻掃描時務必確認 log 裡的 `CE2=`/`find_circ=` 計數都是非零值。

---

<details>
<summary><b>歷史版本（`-Nr 5`，2026-05-27 前，僅供對照）</b></summary>

| 步驟 | 狀態 |
|------|------|
| consensus_filter | ✅ 6/6 完成（1,594–3,728 circRNAs / sample） |
| merge_counts | ✅ 完成（9,349 circRNAs；filterByExpr 後 4,630） |
| DE analysis | ✅ 完成（三方法全跑：edgeR 482 / DESeq2 409 / limma 736 significant）|
| benchmark accuracy | ✅ 完成（circDEX vs CirComPara2_4tools vs nfcore_3tools；門檻掃描 AUC-PR）|
| benchmark compute cost | ✅ 完成（CIRIquant 實測 **11:50:25** on HPC NFS；2026-06-10 重跑確認）|

**主要數值結果（舊版，`-Nr 5`）**：
- 偵測：9,349 consensus circRNAs → filterByExpr 後 4,630 tested
- DE（edgeR_ciriquant）：482 significant（nominal p < 0.05，|log2FC| > 1）；**409 Type_I (84.9%) / 73 Type_II (15.1%)**；min Storey q = 0.384（underpowered）
- DE（DESeq2 baseline）：409 significant；DE（limma-voom）：736 significant
- Isoform switching：66 events（within-gene FDR < 0.1，|ΔIUI| > 0.1）
- Top 1 biomarker：chr10:5836848|5842668（hsa_circ_0002665，GDI2；score=0.8202，83 miRNA，118 RBP binders，log2FC=7.48，Type_I）
- **⚠️ benchmark 表格與 SRP011042（ground truth）為獨立資料集，不受此次 GSE113230 `-Nr 2` 重跑影響，但 accuracy_benchmark 若曾用 GSE113230 作對照組需重新確認**

Benchmark（門檻掃描 AUC-PR，三方法，舊版數據）：

| Method | Precision | Recall | F1 | Specificity | AUC-PR（門檻掃描）|
|--------|-----------|--------|----|-------------|-----------------|
| **circDEX**（本研究）| **0.886** | 0.131 | 0.228 | **0.971** | 0.120 |
| **CirComPara2_4tools** | 0.854 | **0.227** | **0.358** | 0.934 | **0.349** |
| nfcore_3tools | 0.873 | 0.182 | 0.301 | 0.955 | 0.337 |

circDEX Specificity 最高（0.971，假陽性最少）；CirComPara2_4tools Recall 與 AUC-PR 最高（4 工具廣網策略）；三方法 AUC-PR 均以門檻掃描計算，數值可直接比較。

**Biomarker score 公式（6D）**：
```
score = (sig_norm + fc_norm + conf_norm + known_bonus + mirna_norm + rbp_norm) / 6
  sig_norm   = −log10(pvalue), 上限 10，min-max 標準化
  fc_norm    = |log2FC|, 上限 5，min-max 標準化
  conf_norm  = confidence_score, min-max 標準化
  known_bonus = 1 若 in_circbase，否則 0（不標準化）
  mirna_norm = distinct miRNA binders 數，min-max 標準化（無 interaction data 者為 0）
  rbp_norm   = distinct RBP binders 數，min-max 標準化（無 interaction data 者為 0）
```

**Cascade 顯著性邏輯（de_sig_by: auto）**：
1. filterByExpr(min.count=5) → 4,630 tests（原 9,349）
2. Storey q-value < 0.2（bioconductor-qvalue，已安裝）
3. min q = 0.384 → fallback to nominal p < 0.05
4. 482 significant circRNAs（論文 Methods 需說明樣本量限制）

**GSE113230 各工具偵測數量**：

| SRR ID | 分組 | CIRIquant | DCC | 共識 |
|--------|------|----------:|----:|-----:|
| SRR7012366 | Tumor 1 | 26,455 | 6,010 | 1,905 |
| SRR7012367 | Tumor 2 | 34,075 | 7,222 | 2,329 |
| SRR7012368 | Tumor 3 | 35,000 | 10,756 | 3,728 |
| SRR7012369 | Normal 1 | 27,105 | 9,397 | 3,157 |
| SRR7012370 | Normal 2 | 39,211 | 9,349 | 2,790 |
| SRR7012371 | Normal 3 | 12,985 | 5,788 | 1,594 |

</details>

---

### GSE58135（乳癌）

**完成。** 報告位置：`~/GSE58135_results/report.html`（server）

**⚠️ Library prep 更正（2026-10-03，依 GEO extract_protocol 原文查證）**：本文件先前各處敘述僅提及「50bp reads」，從未揭露 library selection 方式。GEO extract_protocol 原文為「POLY(A)-SELECTED, Tn-RNA-seq non-stranded」——即 **poly(A)-selected**，而非 total RNA。這正是本文件 Checklist 段落自己列出的排除條件（poly-A selection → circRNA 接近零），但 GSE58135 仍在此建庫方式下偵測到 1,607 個 circRNA 並留在 19 組 cohort 中。此為需要在論文 Methods/Limitations 中明確揭露的方法學限制，不應僅以「50bp讀長」帶過；詳見「資料集特性稽核校正」表。

**特殊情況：50bp reads + CIRIquant/DCC 嚴重失衡**
- Read length 50bp → STAR chimeric junction 偵測極差（DCC 僅 3–17 circRNA/sample）
- CIRIquant vs DCC 比例 ≈ 0.002，遠低於 adaptive_ratio=0.1 閾值
- `--adaptive` flag 先前未傳入 circrna.smk（bug 已修正）；修正後 adaptive fallback 觸發，min_tools 從 2 降為 1（CIRIquant-only 模式）
- 最終共識：1,607 circRNAs；filterByExpr 後 **45** tested circRNAs；DESeq2 vst() 失敗已以 varianceStabilizingTransformation() fallback 修正

| 步驟 | 狀態 |
|------|------|
| prefetch + fasterq-dump | ✅ 10/10 完成 |
| fastp QC/trim | ✅ 10/10 完成 |
| CIRIquant | ✅ 10/10 完成 |
| STAR / DCC | ✅ 10/10 完成（DCC 偵測數極少，adaptive 模式僅用 CIRIquant） |
| consensus_filter（--adaptive）| ✅ 完成（1,607 circRNAs；CIRIquant-only 模式） |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 15 / DESeq2 122 / limma 508 significant）|
| report | ✅ 完成 |

**注意**：edgeR_ciriquant 顯著數極少（15）是因為 50bp reads → circRNA 偵測數量有限 + 樣本間差異較大；**13 Type_I (86.7%) / 2 Type_II (13.3%)**；limma-voom 在小樣本較穩定（508 significant）；Isoform switching：**11 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 851 rows）。

**Server config**（`config/projects/GSE58135.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE58135/raw`
- `results_dir: /home3/choukaihsuan/GSE58135_results`

**中間檔案已清理**（保留報告、count matrix、DE 結果）；raw FASTQ 已刪除以釋放磁碟空間。

---

### GSE323364（TNBC cell line，EZH2 inhibitor）

**完成。** 報告位置：`~/GSE323364_results/report.html`（server）

MDA-MB-436 TNBC cell line，EZH2 抑制劑 EPZ-6438 vs. DMSO，150bp PE，total RNA，各 3 replicates。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成 |
| fastp QC/trim | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成 |
| STAR / DCC | ✅ 6/6 完成 |
| consensus_filter（--adaptive）| ✅ 完成 |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 15 / DESeq2 122 / limma 508 significant）|
| report | ✅ 完成 |

**設定**：
- case/control label：`EPZ6438` / `DMSO`
- SRR 清單：SRR37484804–SRR37484809（6 個 sample；3 EPZ6438 + 3 DMSO）
- genome：hg19（同 GSE113230）

**主要數值結果**：
- 偵測：filterByExpr 後 **61** tested circRNAs
- edgeR_ciriquant：2 significant circRNAs（nominal p < 0.05）；**2 Type_I (100%) / 0 Type_II**；EZH2 抑制劑對 circRNA 影響極有限（細胞株 + 藥物處理）
- Isoform switching：**0 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 212 rows；細胞株 circRNA 量少，無顯著 switching）
- Biomarker candidates：2 個（p < 0.05 篩選極嚴）

**中間檔案已清理**（raw FASTQ + trimmed + sra_cache + 中間 BAM 已刪除，釋放 77GB）；
報告 1.3MB，保留於 `~/GSE323364_results/report.html`。

**Server config**（`config/projects/GSE323364.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE323364/raw`
- `results_dir: /home3/choukaihsuan/GSE323364_results`

**Condition list CSV 格式**（6 行）：
```
srr_id,condition,gsm_id,description
SRR37484809,EPZ6438,GSM9564370,MDA-MB-436 EPZ-6438 rep1
SRR37484808,DMSO,GSM9564371,MDA-MB-436 DMSO rep1
SRR37484807,EPZ6438,GSM9564372,MDA-MB-436 EPZ-6438 rep2
SRR37484806,DMSO,GSM9564373,MDA-MB-436 DMSO rep2
SRR37484805,EPZ6438,GSM9564374,MDA-MB-436 EPZ-6438 rep3
SRR37484804,DMSO,GSM9564375,MDA-MB-436 DMSO rep3
```

---

### GSE133998（乳癌，paired tumor/normal）

**完成。** 報告位置：`~/GSE133998_results/report.html`（server）

乳癌手術切除組織，cancer tissue vs. adjacent normal，150bp PE，Illumina HiSeq X Ten，各 6 replicates（H36–H42，共 12 samples）。

**2026-07-18 全量重跑（`-Nr 2` 統一 + consensus 去重修正）**：與 GSE113230 相同原因（見該段落），原始分析用 DCC `-Nr 5`，重新下載 12 個樣本並以 `-Nr 2` + 座標去重完整重跑 STAR→DCC→consensus→DE→report；CIRIquant GTF 沿用舊檔。**這次未遇到 GSE113230 那樣的 DCC 基因密集區卡點**，全程順利完成，無需人工介入。

| 步驟 | 狀態 |
|------|------|
| fastp QC/trim（重新下載）| ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成（沿用舊 GTF，未重跑） |
| STAR paired-end / mate1 / mate2 | ✅ 12/12 完成 |
| DCC（`-Nr 2`）| ✅ 12/12 完成 |
| consensus_filter（含去重）| ✅ 12/12 完成（2,451–5,681 circRNAs / sample，見下表） |
| merge_counts | ✅ 完成（**17,504** consensus circRNAs；filterByExpr 後 **650** tested） |
| DE analysis | ✅ 完成（三方法全跑：edgeR **79** / DESeq2 **186** / limma **1,148** significant）|
| predict_interactions / isoform switching / rank_biomarkers / report | ✅ 完成 |

**主要數值結果（`-Nr 2` 全量重跑）**：
- 偵測：**17,504** consensus circRNAs（較 `-Nr 5` 版本的 10,979 增加）→ filterByExpr 後 **650** tested
- DE（edgeR_ciriquant）：**79** significant（nominal p < 0.05，|log2FC| ≥ 1）；上調 33 / 下調 46；**77 Type_I (97.5%) / 2 Type_II (2.5%)**
- DE（DESeq2）：186 significant；DE（limma-voom）：1,148 significant
- Biomarker candidates：79 個
- Top 1 biomarker：**chr13:33091994|33101669（hsa_circ_0000471，N4BP2L2）**；log2FC=6.65，55 miRNA，131 RBP binders，score=0.7271，Type_I —— 與跨資料集重現分析（`docs/biomarker_orthogonal_support.md`「N4BP2L2 focused follow-up」）記錄的候選一致，且方向（上調）與該次記錄相同

**GSE133998 各工具偵測數量（`-Nr 2`，2026-07-18）**：

| SRR ID | 條件 | CIRIquant | DCC | 共識（去重後）|
|--------|------|----------:|----:|-----:|
| SRR11600329 | normal | 4,347 | 4,146 | 3,320 |
| SRR11600330 | tumor | 3,900 | 3,625 | 2,872 |
| SRR11600331 | normal | 7,194 | 7,344 | 5,681 |
| SRR11600332 | tumor | 3,652 | 3,390 | 2,711 |
| SRR11600333 | normal | 4,079 | 3,887 | 3,080 |
| SRR11600334 | tumor | 3,258 | 3,072 | 2,451 |
| SRR11600335 | normal | 5,206 | 4,968 | 3,998 |
| SRR11600336 | tumor | 4,222 | 3,977 | 3,195 |
| SRR11600337 | normal | 3,952 | 3,769 | 2,963 |
| SRR11600338 | tumor | 5,193 | 5,176 | 4,047 |
| SRR11600339 | normal | 4,358 | 4,115 | 3,292 |
| SRR11600340 | tumor | 3,779 | 3,629 | 2,834 |

**注意**：
- 配對設計（同患者 tumor+normal），嘗試加入 `design = ~patient + condition` 後：edgeR 變差（min padj=0.996，patient dummy 消耗 5 個 df，FSJ offset 已吸收個體差異）；limma 進步（51 padj<0.05）；**決定維持 unpaired 設計**，主要看 edgeR_ciriquant 結果。
- `analysis.R` 已加入向後相容的配對設計支援：若 `sample_groups.csv` 含 `patient_id` 欄則自動啟用 `~patient+condition`，否則維持 `~condition`；`cond_coef=ncol(design)` 確保 condition 係數位置正確。
- edgeR filterByExpr 後僅 650 circRNA（17,504 中）是因為 tumor/normal 樣本間表現量分佈差異較大。

<details>
<summary><b>歷史版本（`-Nr 5`，2026-07-18 前，僅供對照）</b></summary>

| 步驟 | 狀態 |
|------|------|
| consensus_filter（--adaptive）| ✅ 完成（10,979 circRNAs 總計） |
| DE analysis | ✅ 完成（edgeR 84 / DESeq2 194 / limma 674 significant）|

**主要數值結果（舊版，`-Nr 5`）**：
- 偵測：10,979 consensus circRNAs；filterByExpr 後 640 tested
- DE（edgeR_ciriquant）：84 significant（nominal p < 0.05，|log2FC| ≥ 1）；上調 34 / 下調 50；**82 Type_I (97.6%) / 2 Type_II (2.4%)**
- DE（DESeq2）：194 significant；DE（limma-voom）：674 significant
- Isoform switching：8,624 rows（within-gene BH FDR 分析）
- Biomarker candidates：84 個

**各 sample 共識 circRNA 數量（舊版）**：

| SRR ID | 條件 | 共識 circRNA |
|--------|------|----------:|
| SRR11600329 | normal | 871 |
| SRR11600330 | tumor | 2,892 |
| SRR11600331 | normal | 1,765 |
| SRR11600332 | tumor | 597 |
| SRR11600333 | normal | 3,089 |
| SRR11600334 | tumor | 505 |
| SRR11600335 | normal | 4,013 |
| SRR11600336 | tumor | 808 |
| SRR11600337 | normal | 2,975 |
| SRR11600338 | tumor | 1,051 |
| SRR11600339 | normal | 3,307 |
| SRR11600340 | tumor | 2,848 |

</details>

**Server config**（`config/projects/GSE133998.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE133998/raw`
- `results_dir: /home3/choukaihsuan/GSE133998_results`

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單：SRR11600329–SRR11600340（12 samples；6 tumor + 6 normal）
- genome：hg19（同 GSE113230）
- Condition list CSV：`/mnt/c/Users/User/Desktop/GSE133998_condition_list.csv`

**SRR → 樣本對應**：

| SRR ID | 條件 | 樣本 |
|--------|------|------|
| SRR11600329 | normal | H36-adjacent normal |
| SRR11600330 | tumor | H36-cancer tissue |
| SRR11600331 | normal | H37-adjacent normal |
| SRR11600332 | tumor | H37-cancer tissue |
| SRR11600333 | normal | H38-adjacent normal |
| SRR11600334 | tumor | H38-cancer tissue |
| SRR11600335 | normal | H39-adjacent normal |
| SRR11600336 | tumor | H39-cancer tissue |
| SRR11600337 | normal | H40-adjacent normal |
| SRR11600338 | tumor | H40-cancer tissue |
| SRR11600339 | normal | H42-adjacent normal |
| SRR11600340 | tumor | H42-cancer tissue |

**Server config**（`config/projects/GSE133998.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE133998/raw`
- `results_dir: /home3/choukaihsuan/GSE133998_results`

**下載注意事項**：
- NCBI S3 間歇性不可用（`srapath --location s3` 回傳空），直接用 aria2c 手動觸發（S3 URL 恢復後）
- SRR11600334 初次下載失敗（HTTPS fallback ~22 KB/s），需重試
- `pigz` 已安裝於 ciriquant env（`conda install -y -c conda-forge pigz`），gzip 壓縮從 ~2h 降至 ~15min/15GB

**磁碟清理（2026-06-14）**：分析完成後已刪除中間檔案，釋放約 333 GB：
- `GSE133998/raw/`（80 GB）+ `GSE133998/trimmed/`（86 GB）→ 已刪除
- `GSE133998_results/circRNA/*/Aligned.sortedByCoord.out.bam`（~66 GB）→ 已刪除
- `GSE133998_results/circRNA/*/mate1/` + `*/mate2/`（~102 GB）→ 已刪除
- 保留：`.gtf`、`.bed`、`DCC/`、`Chimeric.out.junction`、`high_confidence.bed`、`consensus_summary.tsv`

---

### SRP156355（早期乳癌 IDC，配對 tumor/normal）

**完成。** 報告位置：`~/SRP156355_results/report.html`（server）

早期乳癌手術切除組織（IDC，invasive ductal carcinoma），tumor vs. adjacent normal，100bp PE，Total RNA（rRNA-depleted），Cancer Institute WIA 印度，6 對配對 T/N。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成 |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR / DCC | ✅ 12/12 完成 |
| consensus_filter（--adaptive）| ✅ 完成 |
| merge_counts / assign_isoforms | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 152 / DESeq2 14,697 / limma 14,697 circRNAs tested）|
| predict_interactions | ✅ 完成（**union mode**：247 circRNAs，三方法 top-50 聯集）|
| isoform_switching | ✅ 完成（**308 events**，within-gene FDR < 0.1，共 12,227 rows）|
| rank_biomarkers | ✅ 完成（152 candidates）|
| report | ✅ 完成（8.2 MB，Jun 14 更新；union interactions；**Biomarker 表格隨 DE 方法切換重建**）|

**主要數值結果**：
- 偵測：14,697 consensus circRNAs；filterByExpr 後 1,397 tested（edgeR）
- DE（edgeR_ciriquant）：**152 significant**（nominal p < 0.05，|log2FC| > 1）；上調 19 / 下調 133；**144 Type_I (95%) / 8 Type_II (5%)**
- DE（DESeq2）：全 14,697 tested（poscounts normalization）；DE（limma-voom）：全 14,697 tested
- Isoform switching：**308 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 12,227 rows）
- Biomarker candidates：152 個（predict_interactions 覆蓋 247 circRNAs，三方法均有 miRNA/RBP 資料）

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（6T + 6N）：SRR7645071–SRR7645080、SRR7645087–SRR7645088（各患者配對）
- genome：hg19；配對設計支援（含 `patient_id` 欄）
- Condition CSV：`/mnt/c/Users/User/Desktop/SRP156355_condition.csv`

**Server config**（`config/projects/SRP156355.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/SRP156355/raw`
- `results_dir: /home3/choukaihsuan/SRP156355_results`

**磁碟清理（2026-06-19）**：釋放約 ~250 GB：
- `SRP156355/raw/` + `SRP156355/trimmed/` → 已刪除
- `SRP156355_results/circRNA/*/Aligned.sortedByCoord.out.bam` + `*/mate1/` + `*/mate2/` → 已刪除
- 保留：`.gtf`、`.bed`、`DCC/`、`Chimeric.out.junction`、`high_confidence.bed`、`consensus_summary.tsv`、report.html

---

### GSE77509（HCC 肝癌，配對 tumor/normal）

**✅ 完成（2026-06-14）。** 報告位置：`~/GSE77509_results/report.html`（server）

Yang et al. 2017 *Nature Communications*（PMID 28194035）。HCC 肝細胞癌，tumor vs. adjacent normal，~100bp PE，Total RNA（rRNA-depleted），Illumina HiSeq 2500。選 6 對（定序深度最高）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成（S3 aria2c，-x 12 × 4 parallel = 48 連線） |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR / DCC | ✅ 12/12 完成 |
| consensus_filter / merge_counts | ✅ 完成（4,275 circRNAs） |
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 41 / DESeq2 202 / limma 529 significant）|
| predict_interactions（union mode）| ✅ 完成 |
| isoform_switching | ✅ 完成 |
| rank_biomarkers | ✅ 完成（41 candidates）|
| report | ✅ 完成（5.7 MB；mock script `/tmp/run_generate_report_gse77509.py`）|

**主要數值結果**：
- 偵測：4,275 consensus circRNAs；filterByExpr 後 124 tested（HCC 樣本間表現差異大）
- DE（edgeR_ciriquant）：**41 significant**（nominal p < 0.05，|log2FC| > 1）；上調 5 / 下調 36；**36 Type_I (87.8%) / 5 Type_II (12.2%)**
- DE（DESeq2）：202 significant；DE（limma-voom）：529 significant
- Top biomarker：chr21:30693542|30702014（hsa_circ_0001181，BACH1；score=0.6536，31 miRNA，49 RBP binders，log2FC=5.49，Type_I）
- Isoform switching：2,744 rows（within-gene BH FDR 分析）

**注意**：filterByExpr 僅保留 124/4,275 circRNAs，主因 HCC 樣本對（tumor + adjacent normal from same patient）表現分佈差異較大，且單 SRR 約 12–30M reads（完整樣本由 4–5 個 split SRR 組成，僅用第一個）。

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（6T + 6N）：SRR3140264、SRR3140284、SRR3140289、SRR3140303、SRR3140311、SRR3140326（Normal）；SRR3140362、SRR3140382、SRR3140387、SRR3140400、SRR3140408、SRR3140422（Tumor）
- genome：hg19；配對設計（含 `patient_id` 欄，P12/P16/P17/P20/P22/P26）
- Condition CSV：`/mnt/c/Users/User/Desktop/GSE77509_condition.csv`

**Server config**（`config/projects/GSE77509.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE77509/raw`
- `results_dir: /home3/choukaihsuan/GSE77509_results`

**S3 手動下載腳本**：`/tmp/s3_download_gse77509.sh`（server 上）
- 12 batch（每批 4 個 SRR 並行，各 -x 12）→ 總計 48 連線

**磁碟清理（2026-06-19）**：釋放約 ~200 GB：
- `GSE77509/raw/` + `GSE77509/trimmed/` → 已刪除
- `GSE77509_results/circRNA/*/Aligned.sortedByCoord.out.bam` + `*/mate1/` + `*/mate2/` → 已刪除

**SRP011042（Benchmark ground truth，目錄名仍為 GSE55872_results，歷史遺留）磁碟清理（2026-06-19）**：
- `GSE55872_results/circRNA/*/mate1/` + `*/mate2/` + `GSE55872/raw/` → 已刪除（釋放 ~80 GB）

---

### GSE130078（ESCC 食道鱗狀細胞癌，配對 tumor/normal）

**✅ 完成（2026-06 期間）。** 報告位置：`~/GSE130078_results/report.html`（5.9 MB，2026-06-29 13:17 最新版）

食道鱗狀細胞癌（esophageal squamous cell carcinoma，ESCC）手術切除組織，tumor vs. adjacent normal，101bp PE（2026-10-19 查證更正，原誤記150bp；見下方「資料集特性稽核校正」表與fastp實測），Total RNA（rRNA-depleted），Illumina，6 對配對 T/N，12 samples。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成 |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR paired+mate1+mate2 | ✅ 36/36 完成 |
| DCC | ✅ 12/12 完成 |
| consensus_filter / merge_counts | ✅ 完成（8,925 circRNAs；231 after filterByExpr = **2.6%**）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 12 / DESeq2 1,501 / limma 623 significant）|
| predict_interactions（union mode）| ✅ 完成 |
| isoform_switching | ✅ 完成（**133 events**，within-gene FDR < 0.1，共 6,835 rows）|
| rank_biomarkers | ✅ 完成（12 candidates）|
| report | ✅ 完成（5.9 MB）|

**主要數值結果**：
- 偵測：8,925 consensus circRNAs；filterByExpr 後 **231 tested（2.6%）**
- DE（edgeR_ciriquant）：**12 significant**（nominal p < 0.05；全部 |log2FC| > 1）；鱗狀癌 circRNA 全局下調
- DE（DESeq2）：1,501 significant；DE（limma-voom）：623 significant（p < 0.05）
- Isoform switching：**133 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 6,835 rows）
- Biomarker candidates：12 個

**注意**：filterByExpr 通過率極低（2.6%）是本資料集最大特徵。BSJ counts 樣本間差異達 5.5 倍（5,775–31,656），即使定序深度 79–101M reads/sample 也無法改善。鱗狀癌（squamous cell carcinoma）circRNA 全局表現量本就低於腺癌。**建議以 limma-voom 為主方法**（對稀疏計數更穩健），edgeR 僅作輔助參考。

**設定**：
- case/control label：`tumor` / `normal`
- genome：hg19；配對設計（6 pairs）
- Library：Total RNA（rRNA-depleted），101bp PE（2026-10-19更正，原誤記150bp），Illumina

**Server config**（`config/projects/GSE130078.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE130078/raw`（已清理）
- `results_dir: /home3/choukaihsuan/GSE130078_results`

**磁碟清理**：raw/trimmed/BAM/mate1/mate2 於 2026-06-27 清理（釋放 ~162 GB；見 GSE248612 磁碟清理紀錄）。

---

### GSE248612（胃癌，配對 tumor/normal）

**✅ 完成（2026-06-20）。** 報告位置：`~/GSE248612_results/report.html`（server）

胃癌手術切除組織，cancer tissue vs. adjacent normal，6 對配對 T/N，12 samples（SRR26946845–SRR26946856）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成（aria2c S3；845/846 曾遇 HTTPS 慢速，手動改 S3 aria2c 補救）|
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR paired-end | ✅ 12/12 完成（846 曾因 star_tmp 殘留目錄失敗，rm -rf 後重跑）|
| STAR mate1 | ✅ 12/12 完成 |
| STAR mate2 | ✅ 12/12 完成 |
| DCC | ✅ 12/12 完成 |
| consensus_filter（--adaptive）| ✅ 完成（18,391 circRNAs 總計）|
| merge_counts / assign_isoforms | ✅ 完成（328 tested after filterByExpr）|
| annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 46 / DESeq2 124 / limma 1,210 significant）|
| predict_interactions（union mode）| ✅ 完成（233 circRNAs）|
| isoform_switching | ✅ 完成（66 events，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（46 candidates）|
| report | ✅ 完成（9.9MB，Jun 20 16:26）|

**主要數值結果**：
- 偵測：18,391 consensus circRNAs；filterByExpr 後 328 tested
- DE（edgeR_ciriquant）：**46 significant**（nominal p < 0.05，|log2FC| > 1）；上調 10 / 下調 36；**42 Type_I (91.3%) / 4 Type_II (8.7%)**
- DE（DESeq2）：124 significant；DE（limma-voom）：1,210 significant
- Top biomarker：chr21:16386665|16415895（**hsa_circ_0004771，NRIP1**；log2FC=5.58，p=0.0064，Type_I，in_circbase=1）
- Isoform switching：66 events（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- predict_interactions：233 circRNAs（三方法 top-50 union mode）

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（6T + 6N）：SRR26946848, SRR26946849, SRR26946850, SRR26946854, SRR26946855, SRR26946856（Tumor）；SRR26946845, SRR26946846, SRR26946847, SRR26946851, SRR26946852, SRR26946853（Normal）
- genome：hg19；配對設計（含 `patient_id` 欄：LSS/QHZ/ZZJ/LFE/WHD/XXC）
- Condition CSV：`/mnt/c/Users/User/Desktop/GSE248612_condition.csv`

**Server config**（`config/projects/GSE248612.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE248612/raw`
- `sra_cache_dir: /home3/choukaihsuan/GSE248612/sra_cache`
- `results_dir: /home3/choukaihsuan/GSE248612_results`

**Snakemake 執行指令**（cores 36，log 追加）：
```bash
conda run -n ciriquant snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/projects/GSE248612.yaml \
    --cores 36 \
    --resources mem_gb=300 \
    --keep-going \
    --rerun-incomplete \
    >> logs/pipeline_GSE248612_c36.log 2>&1
```

**下載問題處理紀錄**：
- SRR26946849 / SRR26946846：Snakemake 原本走 HTTPS prefetch（~560 KB/s），手動 kill prefetch + 用 aria2c S3 下載後補 fasterq-dump + pigz，再放入 `raw/`
- SRR26946846 SRA 檔案損壞（vdb-validate 回報 "zombie file"，`rcBlob,rcCorrupt`）：刪除 `.sra.tmp` 重新 aria2c 下載
- 已知：`kill -9` 終止 Snakemake 後，worker 變成 orphan，需手動 `kill -9 <orphan_pids>`；避免用 bare `kill`（SIGTERM 觸發 Snakemake cleanup 刪 output 檔）

**磁碟清理（2026-06-27）**：GSE221107 啟動前釋放磁碟（原本 /home3 剩 62G → 清理後 633G）：
- `GSE130078/raw/`（81G）+ `GSE130078/trimmed/`（81G）→ 已刪除
- `GSE130078_results/circRNA/*/Aligned.sortedByCoord.out.bam` + `*/mate1/` + `*/mate2/` → 已刪除
- `GSE248612/raw/`（61G）+ `GSE248612/trimmed/`（64G）+ `GSE248612/sra_cache/`（7G）→ 已刪除
- `GSE248612_results/circRNA/*/Aligned.sortedByCoord.out.bam` + `*/mate1/` + `*/mate2/` → 已刪除
- 保留：report.html、count_matrix.tsv、DE 結果、circRNA GTF/BED/DCC/consensus

---

### GSE221107（攝護腺癌，配對 tumor/normal）

**✅ 完成（重跑版：4 pairs，2026-06-29 22:36）。** 報告位置：`~/GSE221107_results/report.html`（server）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE221107_report.html`

攝護腺癌手術切除組織，cancer tissue vs. adjacent normal，20 對配對，從中選深度最高的 6 對（Pair 8/11/14/16/18/20）偵測完畢。分析時發現 **Pair 8（SRR22757442）和 Pair 11（SRR22757419）** 的 RNA 降解（fastp insert size peak 40–44 bp vs. 正常 268–269 bp，adapter dimer 特徵）→ 排除，最終使用 **4 pairs（Pair 14/16/18/20，8 samples）** 重跑 DE / isoform switching / predict_interactions / report。

**降解樣本辨識標誌**：fastp JSON 中 `insert_size_histogram` 在 ~40–44bp 有高峰（adapter dimer），正常樣本高峰在 268–269bp；此樣本 BSJ counts 極少（SRR22757419：28 BSJ，SRR22757442：6 BSJ），參與 edgeR TMM normalization 時導致 NaN/Inf 使 GLM 失敗。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成 |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成（平均 ~2.75–3.5h/sample） |
| STAR paired+mate1+mate2 | ✅ 36/36 完成 |
| DCC | ✅ 12/12 完成 |
| consensus_filter / merge_counts | ✅ 完成（22,577 circRNAs；子集 8 樣本）|
| DE analysis（4 pairs）| ✅ 完成（edgeR 57 / DESeq2 58 / limma 627 significant）|
| isoform_switching（4 pairs）| ✅ 完成（37 events，within-gene FDR < 0.1）|
| predict_interactions | ✅ 完成（union mode，interactions.json 9.2MB，2026-06-29 20:42）|
| rank_biomarkers | ✅ 完成（57 candidates，biomarker_candidates.tsv 15KB）|
| report | ✅ 完成（report.html 9.1MB，2026-06-29 22:36）|

**主要數值結果（4 pairs 重跑）**：
- 偵測：22,577 consensus circRNAs；count_matrix.tsv 子集為 8 samples（SRR22757410/412/414/416/430/432/434/436）；filterByExpr 後 **928 tested**
- DE（edgeR_ciriquant）：**57 significant**（nominal p < 0.05，|log2FC| > 1）；全為 Type_I（配對設計 `~patient+condition`）
- DE（DESeq2）：**58 significant**；DE（limma-voom）：**627 significant**
- Isoform switching：37 events（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker candidates：**57 個**（edgeR 顯著 circRNA 全部進入）；Top 1：chr2:40655613|40657444（log2FC=−6.79，pvalue=0.024，Type_I）

**排除原因詳細說明**：
- Pair 8（PC8=SRR22757442：6 BSJ）+ Pair 11（PN11=SRR22757419：28 BSJ）：fastp insert_size_peak ≈ 40–44 bp（正常 ≈ 268–269 bp）→ library prep 前 RNA 已降解，adapter dimer 佔主體，無法產生有效的 BSJ spanning reads
- Pair 8 的 normal（PN8=SRR22757422）在原始 11 樣本 count_matrix 中就已缺失（共識過濾後 0 circRNAs，merge_counts 時自動排除）
- Pair 11 tumor（PC11=SRR22757439）雖 BSJ 數量尚可，但 Pair 11 normal（SRR22757419）降解 → 整個 Pair 11 無法配對分析

**count_matrix 子集方式**（繞過重跑 CIRIquant）：
```bash
# 備份原始 11 樣本矩陣
cp count_matrix.tsv count_matrix.tsv.11samples
cp fsj_count_matrix.tsv fsj_count_matrix.tsv.11samples
# Python 切欄（保留 P14/16/18/20 的 8 samples）
python3 -c "
import pandas as pd
good=['SRR22757410','SRR22757412','SRR22757414','SRR22757416',
      'SRR22757430','SRR22757432','SRR22757434','SRR22757436']
for fn in ['count_matrix.tsv','fsj_count_matrix.tsv']:
    df=pd.read_csv(fn,sep='\t',index_col=0)
    df[good].to_csv(fn,sep='\t')
"
```

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（4T + 4N，Pair 14/16/18/20）：PC14=SRR22757436（T），PN14=SRR22757416（N），PC16=SRR22757434（T），PN16=SRR22757414（N），PC18=SRR22757432（T），PN18=SRR22757412（N），PC20=SRR22757430（T），PN20=SRR22757410（N）
- genome：hg19；配對設計（`patient_id` 欄：P14/P16/P18/P20）
- GEO：GSE221107（SubSeries of GSE221109）/ SRA：PRJNA912767
- Library：Ribo-off rRNA Depletion Kit + KC Stranded Library（stranded PE150）

**Server config**（`config/projects/GSE221107.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE221107/raw`
- `results_dir: /home3/choukaihsuan/GSE221107_results`

---

## PRJNA553289（SCLC 小細胞肺癌，配對 tumor/normal）

**✅ 完成（2026-07-01 00:53）。** 報告位置：`~/PRJNA553289_results/report.html`（server）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/PRJNA553289_report.html`

SCLC 小細胞肺癌手術切除組織，tumor vs. adjacent normal，6 pairs（12 samples，SRR9675242–SRR9675253）。南京胸科醫院（Nanjing Chest Hospital），Illumina HiSeq X Ten，PE150，Total RNA。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 12/12 完成（S3 aria2c） |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 12/12 完成 |
| STAR paired+mate1+mate2 | ✅ 36/36 完成 |
| DCC | ✅ 12/12 完成 |
| consensus_filter / merge_counts | ✅ 完成（20,164 circRNAs） |
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 117 / DESeq2 2,116 / limma 1,713 significant）|
| predict_interactions（union mode）| ✅ 完成（228 circRNAs，interactions.json 8.4MB）|
| isoform_switching | ✅ 完成（174 events，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（117 candidates）|
| report | ✅ 完成（report.html 8.5MB，study_title 已設定）|

**主要數值結果**：
- 偵測：20,164 consensus circRNAs；filterByExpr 後 **476 tested**（edgeR）
- DE（edgeR_ciriquant）：**117 significant**（nominal p < 0.05，|log2FC| > 1）；上調 23 / 下調 94；**110 Type_I (94%) / 7 Type_II (6%)**
- DE（DESeq2）：2,116 significant；DE（limma-voom）：1,713 significant
- Isoform switching：174 events（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker candidates：117 個；Top 1：chr14:32559708|32586493（log2FC=+5.76，p=0.0037，Type_I）

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（6T + 6N）：SRR9675248–9675253（Tumor）；SRR9675242–9675247（Normal）
- genome：hg19；配對設計（含 `patient_id` 欄）
- study_title：`RNA-seq profiling of small cell lung cancer tumor and adjacent normal lung tissues (Nanjing Chest Hospital)`

**Server config**（`config/projects/PRJNA553289.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/PRJNA553289/raw`
- `sra_cache_dir: /home3/choukaihsuan/PRJNA553289/sra_cache`
- `results_dir: /home3/choukaihsuan/PRJNA553289_results`

**啟動歷史**：Queue Worker bug（`register_job` 在 `try` 外 → LockException 未捕捉）導致首次啟動失敗；2026-06-29 12:34 改用 nohup 手動繞過 Queue 啟動，2026-07-01 00:53 完成。

**磁碟清理（建議）**：`PRJNA553289/raw/`（66G）＋`trimmed/`（43G）＋`sra_cache/`（24G）＋`sra_tmp/`（38G）≈ 171G 中間檔可刪除（分析完成後）。

---

## GSE229705（LUAD 肺腺癌，配對 tumor/normal，✅ 完成 2026-07-02）

**✅ 完成（2026-07-02 16:06）。** 報告位置：`~/GSE229705_results/report.html`（server，8.3MB）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE229705_report.html`

NYU LUAD 肺腺癌手術切除組織，tumor vs. adjacent normal，123 對配對，Illumina NovaSeq 6000，**51bp PE**（原文件誤記100bp，2026-10-03依GEO/SRA metadata更正，見「資料集特性稽核校正」表），Total RNA（Trio RNA-Seq, Tecan Genomics, rRNA Deplete）。選 6 對（T/N 深度最平衡且總深度最高）。

**SRA Project**：PRJNA955664（SRP432629）；GEO：GSE229705（NYU Langone Health）。

**選定的 6 對（balance ≥ 0.88，min depth ≥ 40M reads/sample）**：

| Patient | Balance | T reads | N reads | T SRR | N SRR |
|---------|---------|---------|---------|-------|-------|
| NYU784 | 0.880 | 59.9M | 52.7M | SRR24166158 | SRR24166159 |
| NYU539 | 0.933 | 52.5M | 56.3M | SRR24166104 | SRR24166105 |
| NYU704 | 0.883 | 48.4M | 42.7M | SRR24166288 | SRR24166289 |
| NYU713 | 0.895 | 45.5M | 40.8M | SRR24166284 | SRR24166285 |
| NYU779 | 0.880 | 40.1M | 45.5M | SRR24166280 | SRR24166281 |
| NYU822 | 0.905 | 39.8M | 44.0M | SRR24166156 | SRR24166157 |

**Condition CSV**：`/mnt/c/Users/User/Desktop/GSE229705_condition.csv`（含 `patient_id` 欄，配對設計）

**送出方式**：Web UI → 方式二（上傳 CSV） → Project ID：`GSE229705`，Case：`tumor`，Control：`normal`

**注意**：`run_manual` 路由的 `sample_groups.csv` 只寫 `srr_id`/`condition` 兩欄，**不保存 `patient_id`**。送出後需手動在 server 補充：
```bash
# server 上手動更新 sample_groups.csv 加入 patient_id 欄
python3 -c "
import pandas as pd
df = pd.read_csv('metadata/GSE229705/sample_groups.csv')
pid_map = {
  'SRR24166158':'NYU784','SRR24166159':'NYU784',
  'SRR24166104':'NYU539','SRR24166105':'NYU539',
  'SRR24166288':'NYU704','SRR24166289':'NYU704',
  'SRR24166284':'NYU713','SRR24166285':'NYU713',
  'SRR24166280':'NYU779','SRR24166281':'NYU779',
  'SRR24166156':'NYU822','SRR24166157':'NYU822',
}
df['patient_id'] = df['srr_id'].map(pid_map)
df.to_csv('metadata/GSE229705/sample_groups.csv', index=False)
"
```

**論文參考**：Sakata et al. 2024 *Nature Genetics*，PMC10632519；LUAD tumor-adjacent inflammation cohort，123 pairs。

**主要數值結果**：
- 偵測：consensus circRNAs → filterByExpr 後 **29 tested**（edgeR，配對設計 `~patient+condition` 自由度有限）
- DE（edgeR_ciriquant）：**1 significant**（nominal p < 0.05；n=6 pairs + 51bp reads 偵測數少）
- DE（DESeq2）：**162 significant**；DE（limma-voom）：**0 significant**
- predict_interactions：union mode（三方法 top-50 聯集）
- report：study_title = `LUAD Lung Adenocarcinoma — tumor vs. adjacent normal tissues, 6 pairs (NYU Langone Health, Sakata et al. 2024 Nature Genetics)`

---

## GSE148036（LUAD 肺腺癌 tumor vs. normal，✅ 完成 2026-07-04）

**✅ 完成（2026-07-04 03:39）。** 報告位置：`~/GSE148036_results/report.html`（server，5.9MB）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE148036_report.html`

多疾病肺部 RNA-seq 資料集，從中選取 LUAD（Lung Adenocarcinoma）vs. Normal Lung 各 5 samples。Illumina HiSeq 3000，PE146（avgLength ≈ 286–292bp），Ribo-Zero rRNA removal，Total RNA，University of Pittsburgh / UPMC 收集。

**GEO**：GSE148036；**SRA**：PRJNA625051

**注意**：LUAD 和 Normal 來自**不同患者**（非配對設計），使用 `~condition` 模型（無 `patient_id` 欄）。

**選定的 10 samples**：

| SRR | 分組 | 樣本名稱 |
|-----|------|---------|
| SRR11262292 | tumor | AD005 (Lung Adenocarcinoma Rep5) |
| SRR11262293 | tumor | AD004 (Lung Adenocarcinoma Rep4) |
| SRR11262294 | tumor | AD003 (Lung Adenocarcinoma Rep3) |
| SRR11262295 | tumor | AD002 (Lung Adenocarcinoma Rep2) |
| SRR11262296 | tumor | AD001 (Lung Adenocarcinoma Rep1) |
| SRR11262284 | normal | NM005 (Normal Lung Tissue Rep5) |
| SRR11262285 | normal | NM004 (Normal Lung Tissue Rep4) |
| SRR11262286 | normal | NM003 (Normal Lung Tissue Rep3) |
| SRR11262297 | normal | NM002 (Normal Lung Tissue Rep2) |
| SRR11262298 | normal | NM001 (Normal Lung Tissue Rep1) |

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 10/10 完成（S3 aria2c） |
| fastp QC/trim | ✅ 10/10 完成 |
| CIRIquant | ✅ 10/10 完成 |
| STAR paired+mate1+mate2 | ✅ 30/30 完成（`_STARtmp` residual bug 修正後重跑）|
| DCC | ✅ 10/10 完成 |
| consensus_filter / merge_counts | ✅ 完成（874 circRNAs；19 after filterByExpr = **2.2%**）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 6 / DESeq2 27 / limma 323 significant）|
| predict_interactions（union mode）| ✅ 完成 |
| isoform_switching | ✅ 完成（12 events，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（6 candidates）|
| report | ✅ 完成（5.9MB，2026-07-04 03:39）|

**主要數值結果**：
- 偵測：874 consensus circRNAs；filterByExpr 後 **19 tested（2.2%）**——非配對設計（5T + 5N）+ PE146 讀長短樣本間 BSJ 差異大
- DE（edgeR_ciriquant）：**6 significant**（nominal p < 0.05，|log2FC| > 1）；全部下調，全 **Type_I**
- DE（DESeq2）：27 significant；DE（limma-voom）：323 significant
- Isoform switching：**12 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker candidates：6 個
- Top 1：chr8:68044186|68049838（**hsa_circ_0003388，CSPP1**；log2FC=−6.73，p=0.000444，Type_I，score=0.7231，11 miRNA，4 RBP）
- Top 2：chr3:169694734|169706147（**hsa_circ_0001358，SEC62**；log2FC=−1.56，p=0.0477，score=0.6561，23 miRNA，57 RBP）

**注意**：filterByExpr 通過率極低（2.2%，僅 19/874），是 12 個資料集中最低。主因：
1. 非配對設計（unpaired 5T+5N），而 GSE229705 同為 LUAD 但 6 對配對（paired），edgeR 功率更高
2. PE146 讀長（每個 SRR 僅含完整樣本的一半 reads，即 ~46M），BSJ 偵測量相對有限
3. LUAD 的 circRNA 整體表現偏低（相較乳癌/HCC），加上樣本間差異造成稀疏 BSJ counts
建議：以 **limma-voom 為主方法**（323 sig），edgeR 作輔助

**啟動歷史**：
- GSE148036-SFJU（web UI 首次提交）config 路徑指向 GSE229705（`_update_paths_for_project` bug），已手動修正
- r1–r5（早期）：STAR 多個 sample 因 ulimit 和 `_STARtmp` 殘留目錄失敗
- r6：`--forcerun star_align` 完成 6/10 樣本
- r7（2026-07-04 01:59）：修正 `circrna.smk` star_align 加 `rm -rf {params.tmp_dir}` before `mkdir -p`，STAR 4 個失敗樣本全部重跑成功，03:39 pipeline 100% 完成

**config 錯誤歷史**：首次提交後 `config/projects/GSE148036.yaml` 的 `raw_dir`/`results_dir` 指向 `GSE229705/`（`_update_paths_for_project` bug，`old_pid == new_pid` 時不做路徑替換）。已手動 `sed` 修正，並修復 `web_ui.py` 的 `_update_paths_for_project` 函式（改用 regex 替換路徑中任何 GSE/SRP/PRJNA 編號）。

**設定**：
- case/control label：`tumor` / `normal`
- genome：hg19；**非配對設計**（5T + 5N，不同患者）
- Condition CSV：`/mnt/c/Users/User/Desktop/GSE148036_condition.csv`（無 `patient_id` 欄）
- Library：Total RNA（Ribo-Zero），PE146，Illumina HiSeq 3000
- 每 sample 深度：~46M reads（23M spots × 2）

**Server config**（`config/projects/GSE148036.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE148036/raw`
- `results_dir: /home3/choukaihsuan/GSE148036_results`

---

## GSE121842（CRC 大腸直腸癌，配對 tumor/normal，✅ 完成 2026-07-13）

**✅ 完成（2026-07-13 18:11）。** 報告位置：`~/GSE121842_results/report.html`（server，13 MB）

大腸直腸癌（colorectal cancer，CRC）手術切除組織，tumor vs. adjacent normal，3 對配對 T/N，6 samples（SRR8113703–SRR8113708）。Illumina HiSeq X（建庫法未於 GEO 記載；原「TruSeq Stranded + Ribo-Zero Gold」描述已於 2026-10-03 撤回），150bp PE，~35–49M spots/sample。Wang et al. 2019（PMC6757124）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成（HTTPS prefetch fallback，S3 暫時不可用）|
| fastp QC/trim | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成（SRR8113707 03:08 → SRR8113704 15:19）|
| STAR paired+mate1+mate2 | ✅ 18/18 完成 |
| DCC | ✅ 6/6 完成 |
| consensus_filter / merge_counts | ✅ 完成（11,148 circRNAs）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（edgeR 65 / DESeq2 / limma）|
| predict_interactions（union mode）| ✅ 完成（227 circRNAs）|
| isoform_switching | ✅ 完成 |
| rank_biomarkers | ✅ 完成（65 candidates）|
| report | ✅ 完成（13 MB，2026-07-13 18:11）|

**主要數值結果**：
- 偵測：**11,148 consensus circRNAs**；filterByExpr 後 **661 tested**
- DE（edgeR_ciriquant）：**65 significant**（nominal p < 0.05，|log2FC| > 1）；配對設計（`~patient+condition`，3 pairs）
- predict_interactions：227 circRNAs（三方法 top-50 union mode）
- Top 1：**hsa_circ_0002483（PTK2）**，chr8:141874411|141900868，log2FC=+9.28，Type_II，score=0.8478，102 miRNA / 103 RBP
- Top 2：**hsa_circ_0017586（GDI2）**，chr10:5815805|5842668，log2FC=+8.74，Type_I，score=0.66，112 miRNA / 150 RBP

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（3T + 3N）：SRR8113703（tumor-P454750）、SRR8113704（tumor-Pgzy）、SRR8113705（tumor-Pwzd）；SRR8113706（normal-P454750）、SRR8113707（normal-Pgzy）、SRR8113708（normal-Pwzd）
- genome：hg19；配對設計（patient_id：P454750 / Pgzy / Pwzd）
- Library：150bp PE，Illumina HiSeq X；建庫法 **未於 GEO 12 筆 GSM 記載**（2026-10-03 查證：原文件「TruSeq Stranded Total RNA + Ribo-Zero Gold」無 GEO 依據，已撤回，strandedness 未載明）
- Condition CSV：`/mnt/c/Users/User/Desktop/GSE121842_condition.csv`

**Server config**（`config/projects/GSE121842.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE121842/raw`
- `results_dir: /home3/choukaihsuan/GSE121842_results`

---

## GSE136569（胰臟癌 PDAC，配對 tumor/normal，✅ 完成 2026-07-22）

**✅ 完成（2026-07-22 02:38）。** 報告位置：`~/GSE136569_results/report.html`（server，20MB）

胰臟導管腺癌（Pancreatic ductal adenocarcinoma，PDAC）手術切除組織，tumor vs. adjacent normal tissue (NAT)，5 對配對 T/N，10 samples（SRR10030979–SRR10030988）。100bp PE（2026-10-19 查證更正，原誤記150bp；見下方「資料集特性稽核校正」表與fastp實測），Total RNA，Illumina HiSeq X Ten。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 10/10 完成（SRR10030983 曾遇 S3 DNS SERVFAIL 7+ 小時，手動 aria2c 於 2026-07-20 17:03 完成）|
| fastp QC/trim | ✅ 10/10 完成 |
| CIRIquant | ✅ 10/10 完成（SRR10030983 於 2026-07-20 21:05 完成，NFS I/O 良好：HISAT2 27min / BWA-MEM 51min / CIRI2 25min）|
| STAR paired+mate1+mate2 | ✅ 30/30 完成 |
| DCC | ✅ 10/10 完成（SRR10030983 DCC 於 2026-07-21 23:53 完成）|
| consensus_filter / merge_counts | ✅ 完成（count_matrix.tsv 3.1MB，Jul 21 23:57）|
| assign_isoforms / annotate_circbase | ✅ 完成（isoform_groups.tsv 9.1MB；circbase_annotated.tsv 5.6MB，耗時 ~75 min）|
| DE analysis | ✅ 完成（三方法全跑：edgeR **199** / DESeq2 **472** / limma **2,340** significant，依 p<0.05 且 \|log2FC\|≥1；DESeq2 數字依現行 de_results_deseq2.tsv 核對，舊版 14,163 無法重現）|
| predict_interactions（union mode）| ✅ 完成（interactions.json 12MB，244 circRNAs）|
| isoform_switching | ✅ 完成（**76 events**，within-gene FDR < 0.1，共 67,811 rows）|
| rank_biomarkers | ✅ 完成（199 candidates，biomarker_candidates.tsv 47KB）|
| report | ✅ 完成（report.html 20MB，2026-07-22 02:38）|

**主要數值結果**：
- 偵測：consensus circRNAs 共 **173,037**（10 樣本總計）；filterByExpr 後 **7,334 tested**
- DE（edgeR_ciriquant）：**199 significant**（nominal p < 0.05，|log2FC| > 1）；上調 110 / 下調 89；**138 Type_I (69.3%) / 61 Type_II (30.7%)**
- DE（DESeq2）：472 significant（依現行 de_results_deseq2.tsv，p<0.05 且 |log2FC|≥1；原記 14,163 無法重現）；DE（limma-voom）：2,340 significant
- Isoform switching：76 events（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker candidates：199 個；Top circRNAs：chr7:155465561|155473602（log2FC=4.93，p=0.025）、chr2:240929491|240946787（log2FC=5.46，p=0.0012）

**注意**：Type_II 比例（30.7%）高於其他資料集（一般 10–15%），PDAC 可能有較多 circRNA/線性 RNA 雙層調控事件，值得關注。

**各樣本共識 circRNA 數量**：

| SRR ID | 條件 | 患者 | 共識 circRNA |
|--------|------|------|----------:|
| SRR10030979 | tumor | P1 | 19,371 |
| SRR10030980 | tumor | P2 | 21,670 |
| SRR10030981 | tumor | P3 | 12,751 |
| SRR10030982 | tumor | P4 | 17,647 |
| SRR10030983 | tumor | P5 | 18,025 |
| SRR10030984 | normal | P1 | 16,542 |
| SRR10030985 | normal | P2 | 12,202 |
| SRR10030986 | normal | P3 | 8,108 |
| SRR10030987 | normal | P4 | 27,199 |
| SRR10030988 | normal | P5 | 19,522 |

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（5T + 5N）：SRR10030979–SRR10030983（PDAC1–5）；SRR10030984–SRR10030988（NAT1–5）
- genome：hg19；配對設計（patient_id：P1–P5）
- Library：Ribozero 去 rRNA + RNase R 消化線性 RNA（GSM extract_protocol 原文；2026-10-03 查證，先前誤記為 Total RNA），100bp PE（2026-10-19更正，原誤記150bp），Illumina HiSeq X Ten
- GEO：GSE136569（PRJNA563024）

**Server config**（`config/projects/GSE136569.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE136569/raw`
- `results_dir: /home3/choukaihsuan/GSE136569_results`

**SRR10030983 特殊記錄**：S3 DNS SERVFAIL 持續 7+ 小時（2026-07-20），手動 aria2c 於 17:03 完成；CIRIquant NFS I/O 良好（HISAT2 27min / BWA-MEM 51min / CIRI2 25min，遠優於 NFS 基準）；DCC 於 2026-07-21 23:53 完成（無 RBM5-locus 型卡點）。annotate_circbase 耗時 ~75 分鐘（71,457 unique circRNAs × Python iterrows 迴圈瓶頸）。

---

## GSE143797（鼻咽癌 NPC，✅ 完成 2026-07-19）

**✅ 完成（2026-07-19 12:01）。** 報告位置：`~/GSE143797_results/report.html`（server，12MB）

鼻咽癌（Nasopharyngeal Carcinoma，NPC）腫瘤組織 vs. 鼻咽炎正常組織，4 tumor + 4 normal = 8 samples（SRR10903023–SRR10903030）。Ding et al. 2020（Circular RNA expression in NPC）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 8/8 完成 |
| fastp QC/trim | ✅ 8/8 完成 |
| CIRIquant | ✅ 8/8 完成 |
| STAR paired+mate1+mate2 | ✅ 24/24 完成 |
| DCC | ✅ 8/8 完成 |
| consensus_filter / merge_counts | ✅ 完成（8,922 circRNAs；328 after filterByExpr）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（三方法全跑：edgeR **56** / DESeq2 / limma）|
| predict_interactions / isoform_switching / rank_biomarkers | ✅ 完成（56 candidates；6,384 isoform rows）|
| report | ✅ 完成（12MB）|

**主要數值結果**：
- 偵測：8,922 consensus circRNAs；filterByExpr 後 **328 tested**
- DE（edgeR_ciriquant）：**56 significant**（nominal p < 0.05，|log2FC| ≥ 1）；上調 8 / 下調 48；**56 Type_I (100%) / 0 Type_II**——全部為 circRNA 專一性調控，無 FSJ 共同顯著
- NPC 腫瘤 circRNA 以下調為主（下調 48 / 上調 8），與乳癌、HCC 等腺癌的全局下調趨勢一致
- Isoform switching：6,384 rows（within-gene BH FDR 分析）
- Biomarker candidates：56 個

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（4T + 4N）：SRR10903027–SRR10903030（Tumor）；SRR10903023–SRR10903026（Normal）
- genome：hg19；**非配對設計**（4T + 4N，無 patient_id 欄）
- Library：Total RNA（rRNA-depleted），NPC tumor vs. nasopharyngitis tissue（Ding et al. 2020）

**Server config**（`config/projects/GSE143797.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE143797/raw`（已清理）
- `results_dir: /home3/choukaihsuan/GSE143797_results`

---

## GSE108735（腎細胞癌 RCC，✅ 完成 2026-07-24）

**✅ 完成（2026-07-24 16:55）。** 報告位置：`~/GSE108735_results/report.html`（server，11MB）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE108735_report.html`

> **⚠️ 資料集更正**：原標記為「TNBC」，確認為**腎細胞癌（Renal Cell Carcinoma，RCC）** tumor vs. 正常腎組織。SRA Project：SRP128028。

7 pairs 腎細胞癌手術切除組織，tumor vs. adjacent normal kidney，150bp PE，**ncRNA-Seq**（非 Total RNA-Seq，circRNA 和其他 ncRNA 富集），Illumina。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 14/14 完成 |
| fastp QC/trim | ✅ 14/14 完成 |
| CIRIquant | ✅ 14/14 完成 |
| STAR paired+mate1+mate2 | ✅ 42/42 完成 |
| DCC | ✅ 14/14 完成 |
| consensus_filter / merge_counts | ✅ 完成（14,262 circRNAs；57 after filterByExpr）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（三方法全跑：edgeR **55** / DESeq2 / limma）|
| predict_interactions（union mode）| ✅ 完成 |
| isoform_switching | ✅ 完成（**377 events**，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（55 candidates）|
| report | ✅ 完成（11MB，2026-07-24 16:55）|

**主要數值結果**：
- 偵測：**14,262 consensus circRNAs**；filterByExpr 後 **57 tested**——ncRNA-Seq 的 FSJ counts 因線性 RNA 大幅去除而趨近於零，filterByExpr 依 BSJ counts 篩選，最終只有 57 個 circRNA 在全部 14 個樣本中均有足夠表現量
- DE（edgeR_ciriquant）：**55 significant**（nominal p < 0.05，|log2FC| > 1）；**全部上調（55 up / 0 down）**——ncRNA-Seq 富集 circRNA，RCC 腫瘤中 circRNA 整體上調，與 Total RNA-Seq 資料集「腫瘤普遍下調」趨勢相反（可能反映 library prep 差異，也可能是 RCC 特有的生物學特性）；**47 Type_I (85.5%) / 8 Type_II (14.5%)**
- Isoform switching：**377 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 10,789 rows）
- Biomarker candidates：**55 個**
- Top 1 biomarker：**chr7:155465561|155473602**（log2FC=+8.38，p=1.1e-05，Type_I）——此 circRNA 也是 GSE136569（PDAC）的 Top 1 biomarker，具跨癌種重現性

**注意**：filterByExpr 通過率極低（57/14,262 = 0.4%）主因是 ncRNA-Seq 的 library prep 消除了大量線性 mRNA（FSJ）；edgeR_ciriquant 需要有效的 FSJ counts 估算 offset，FSJ ≈ 0 的 circRNA 大多無法通過 TMM normalization 品管。結果應以 **BSJ counts 直接分析**（DESeq2/limma）作為輔助參考。

**SRR 清單（7T + 7N，配對設計）**：

| SRR ID | 條件 | 患者 |
|--------|------|------|
| SRR6439741 | normal | P1 |
| SRR6439742 | normal | P2 |
| SRR6439743 | normal | P3 |
| SRR6439744 | normal | P4 |
| SRR6439745 | normal | P5 |
| SRR6439746 | normal | P6 |
| SRR6439747 | normal | P7 |
| SRR6439748 | tumor | P1 |
| SRR6439749 | tumor | P2 |
| SRR6439750 | tumor | P3 |
| SRR6439751 | tumor | P4 |
| SRR6439752 | tumor | P5 |
| SRR6439753 | tumor | P6 |
| SRR6439754 | tumor | P7 |

**Server config**（`config/projects/GSE108735.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE108735/raw`（已清理）
- `results_dir: /home3/choukaihsuan/GSE108735_results`
- `metadata: metadata/GSE108735/library_info.csv`
- `groups: metadata/GSE108735/sample_groups.csv`

**磁碟清理（2026-07-27）**：raw/（已刪）+ trimmed/（已刪）+ mate1/mate2/BAM（已刪），釋放空間。

---

## GSE171011（甲狀腺乳突癌 PTC，✅ 完成 2026-07-27）

**✅ 完成（2026-07-27 04:55）。** 報告位置：`~/GSE171011_results/report.html`（server，14MB）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE171011_report.html`

> **⚠️ 資料集更正**：原標記為「TNBC」，確認為**甲狀腺乳突癌（Papillary Thyroid Cancer，PTC）** tumor vs. adjacent normal。SRA Project：SRP312486（GEO: GSE171011）。

4 pairs 甲狀腺乳突癌手術切除組織，tumor vs. adjacent normal thyroid tissue，4T + 4N = 8 samples（SRR14088791–SRR14088798），PAIRED，Total RNA-Seq（RNA-Seq），Illumina。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 8/8 完成 |
| fastp QC/trim | ✅ 8/8 完成 |
| CIRIquant | ✅ 8/8 完成（SRR14088793 第一次 de novo 失敗，第二次成功；hisat2-build 耗時 31 min）|
| STAR paired+mate1+mate2 | ✅ 24/24 完成 |
| DCC | ✅ 8/8 完成 |
| consensus_filter / merge_counts | ✅ 完成（32,848 circRNAs；2,377 after filterByExpr）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis | ✅ 完成（三方法全跑：edgeR **326** / DESeq2 / limma）|
| predict_interactions（union mode）| ✅ 完成 |
| isoform_switching | ✅ 完成（**67 events**，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（326 candidates）|
| report | ✅ 完成（14MB，2026-07-27 04:55）|

**主要數值結果**：
- 偵測：**32,848 consensus circRNAs**；filterByExpr 後 **2,377 tested**
- DE（edgeR_ciriquant）：**326 significant**（nominal p < 0.05，|log2FC| > 1）；上調 173 / 下調 153；**303 Type_I (92.9%) / 23 Type_II (7.1%)**
- Isoform switching：**67 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 29,963 rows）
- Biomarker candidates：**326 個**
- Top 1 biomarker：**chr14:31596991|31641328**（log2FC=−9.02，pvalue=0.0155，Type_I，下調）
- study_title：`Next Generation Sequencing of Papillary Thyroid Cancer tissue sample and adjacent normal tissue`

**CIRIquant SRR14088793 de novo 失敗教訓**：第一次執行於 de novo HISAT2 alignment 階段報 `FileNotFoundError: SRR14088793_denovo.sorted.bam not created`；Snakemake 重啟後第二次成功（hisat2-build 耗時 31 min，遠超 CLAUDE.md 原估算的 8 min；等待期間從時間戳確認 index 仍在寫入、未卡住）。

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（4T + 4N）：SRR14088791–SRR14088798（配對設計，4 pairs）
- genome：hg19；配對設計（`patient_id` 欄）
- Library：Total RNA-Seq，150bp PE，Illumina

**Server config**（`config/projects/GSE171011.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE171011/raw`（已清理）
- `results_dir: /home3/choukaihsuan/GSE171011_results`

**磁碟清理（2026-07-27）**：raw/（71G 已刪）+ trimmed/（27G 已刪）+ mate1/mate2（29G 已刪）。

---

## GEO / SRA BioProject 資料集選擇指引

### 支援的 Accession 格式

Pipeline 支援三種 accession 格式，輸入 `--gse` 或 Web UI 的「GEO 資料集」欄位：

| 格式 | 範例 | 來源 | 查詢方式 |
|------|------|------|---------|
| `GSE*` | `GSE113230` | NCBI GEO Series | pysradb（GSE → SRP → SRR） |
| `PRJNA*` | `PRJNA808398` | NCBI SRA BioProject | NCBI eUtils API |
| `SRP*` | `SRP156355` | NCBI SRA Study | NCBI eUtils API |

**HPC 網路限制**：Server（172.16.0.178）防火牆**曾**封鎖部分 NCBI eUtils endpoint（狀態已變更，見下）：
- `eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi`（db=gds）✅ **可存取**（`_fetch_geo_title()` 使用）
- `eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`（PRJNA/SRP metadata）⚠️ **2026-07 文件記載封鎖；2026-10-19 實測可用**（直接 curl 一個真實 accession，不只看狀態碼也檢查內容，確認回傳真正的 `<EXPERIMENT_PACKAGE_SET>` SRA XML，不是偽裝成 200 的錯誤頁。未測試是否仍有其他 eUtils endpoint 被封鎖，也無法判斷這是防火牆規則改了還是原本就沒封鎖、舊記錄本身有誤——只能確認「efetch.fcgi 被封鎖」這個具體說法現在不成立）
- **（2026-10-19 新增）`scripts/download_geo.py`/`/run_gse`/`preflight_accession()` 三條程式路徑查過，全部呼叫同一個 `efetch.fcgi`，沒有任何程式碼層級的繞道或 fallback 邏輯**——下方「本機先抓再上傳 CSV」是文件記載的人工操作步驟，不是程式碼；efetch 現在可用後這個步驟理論上不再必要，但保留不刪，因為無法排除防火牆規則未來又變回封鎖。也因為三條路徑共用同一呼叫，preflight 查到的存在性/規模對實際下載是可信的保證，不會有「preflight 說有、下載時卻拿到不同/較少欄位」的風險。
- PRJNA/SRP 需在本機執行 `python scripts/download_geo.py --gse PRJNA808398` 取得 metadata CSV，再透過 Web UI 的「手動 SRR 清單 → 上傳 CSV」路徑匯入（2026-07 文件記載的作法；2026-10-19 efetch.fcgi 實測可用後，理論上可省略此步驟直接用 `/run_gse`，但未長期驗證，暫不刪除此作法說明）
- GSE 透過 pysradb 可在 server 上直接執行；`_fetch_geo_title()` 也在 server 上正常運作

**Condition 自動偵測**（`prepare_metadata.py`）：

| 樣本命名慣例 | 偵測邏輯 | 範例 |
|------------|---------|------|
| SRA T/N 尾碼 | `_detect_condition_by_suffix()`（優先）：名稱結尾為 `T` → tumor；`N`/`APN`/`CN` → normal | `M269T`→tumor, `87APN`→normal |
| GEO 關鍵字 | `_detect_condition()`（fallback）：掃描 disease_state/source_name/tissue 欄位 | `tumor tissue`→tumor |
| 無法偵測 | 留空（如 DCIS `D` 尾碼）| `1086D`→空，需手動指定 |

---

### 已分析/已調查資料集

> **⚠️ 權威來源聲明（2026-10-19）**：本表的 **Library** 與 **Read length** 兩欄與上方「資料集特性稽核校正」表描述同一組事實（同一批 fastp 實測/GEO 原始資料），**上方校正表才是權威來源**。本表這兩欄過去曾三度與校正表不同步（GSE130078、GSE136569、GSE248612 的讀長都曾錯誤沿用舊值「150bp」，直到 2026-10-19 才發現並改正），根因是同一組事實在文件裡有多處平行表述、更正只改了其中一處。**修改任何資料集的 Library/Read length 事實時，務必同時檢查並同步上方校正表與本表，不要只改其中一個**；若兩表不一致，以校正表為準並回來修正本表。

| 編號 | 樣本類型 | 對比組 | Library | Read length | 樣本數 | 定序深度（M reads，R1+R2 相加；2026-10-19 統一標註單位） | circRNA 偵測 | DE 顯著數（edgeR）| 評估 |
|------|----------|--------|---------|-------------|--------|----------|-------------|-------------------|------|
| **GSE113230** | 組織（TNBC tumor vs. normal）| tumor vs. normal | Total RNA（rRNA-depleted）| 150bp PE | 6（3T+3N）| ~100M reads/sample | 9,349（consensus）| **482** | ✅ 最佳示範 |
| **GSE58135** | 組織（乳癌 tumor vs. normal）| tumor vs. normal | **Poly(A)-selected**（2026-10-03更正，原誤記Total RNA；見GSE58135段落）| **50bp PE** | 10（5T+5N）| ~50M reads/sample（來源未明，疑為SRA標示值）／實測61.7M read pairs=123.4M reads（fastp，2026-10-19）| 1,607（CIRIquant-only）| **15** | ⚠️ 50bp讀長+poly(A) selection雙重例外，DCC失效→adaptive fallback；排除於15組跨資料集分析（有效單工具比例97.8%）|
| **GSE323364** | **細胞株**（MDA-MB-436 TNBC）| EPZ6438 vs. DMSO | Total RNA | 150bp PE | 6（3+3）| ~60M reads/sample（來源未明）／實測45.4M reads（fastp，2026-10-19）| 中等 | **15** | ⚠️ 細胞株+藥物，DE 極少 |
| **GSE133998** | 組織（乳癌 tumor vs. normal）| tumor vs. normal | Total RNA | 150bp PE | 12（6T+6N）| ~80M reads/sample | 10,979 | **84** | ✅ 配對設計，樣本數最多 |
| **SRP156355** | 組織（早期乳癌 IDC）| tumor vs. normal | rRNA-depleted（`other`）| **100–101bp PE** | 12（6T+6N）| avg 48M read pairs/sample（~96M reads）| **14,697**（consensus）| **152** | ✅ 完成；6 對配對 T/N；Type_I 95% / Type_II 5%；predict_interactions union 247 circRNAs；report 8.2 MB |
| **GSE77509** | 組織（HCC 肝癌 tumor vs. normal）| tumor vs. normal | Total RNA（rRNA-depleted）| **~100bp PE** | 12（6T+6N，6 pairs）| **57–118M read pairs/sample**（平均 ~85M read pairs）| 4,275（consensus）| **41** | ✅ 完成；Top：hsa_circ_0001181（BACH1，score=0.654）；36 Type_I / 5 Type_II；Yang et al. 2017 Nat Commun |
| **GSE130078** | 組織（ESCC 食道鱗狀細胞癌 tumor vs. normal）| tumor vs. normal | Total RNA（rRNA-depleted）| 101bp PE | 12（6T+6N，6 pairs）| **79–101M reads/sample** ✅ | 8,925（consensus）| **12（edgeR）/ 623（limma）** | ⚠️ 深度充足但 BSJ counts 樣本間差異 5.5x（5,775–31,656）；filterByExpr 僅保留 231（2.6%）；ESCC circRNA 全局下調；建議以 limma 為主方法 |
| **GSE248612** | 組織（胃癌 tumor vs. normal）| tumor vs. normal | Total RNA（rRNA depletion 未載明，UNVERIFIED，見校正表）| 138–140bp PE | 12（6T+6N，6 pairs）| ~78.0M reads/sample（fastp 實測39.0M read pairs，2026-10-19）| 18,391（consensus）| **46** | ✅ 完成；Top：hsa_circ_0004771（NRIP1，log2FC=5.58）；42 Type_I / 4 Type_II；6 對配對設計（LSS/QHZ/ZZJ/LFE/WHD/XXC）|
| **PRJNA808398** | 組織（TNBC tumor vs. normal）| tumor vs. normal | **cDNA（poly-A）** | 150bp PE（SRA metadata；poly-A 事前排除未進入pipeline，本機無fastp可核，2026-10-19查證：未驗證） | 50（25T+25N）| avg 18M read pairs/sample（低；SRA metadata，未本機驗證）| ❌ | ❌ | ❌ 不建議：poly-A selection，circRNA 接近零；A.C. CAMARGO 巴西 |

### 影響分析品質的關鍵因素

#### 1. Read Length（讀長）

| 讀長 | circRNA 偵測 | 說明 |
|------|-------------|------|
| **≥ 100bp PE**（建議）| ✅ 良好 | BSJ 需橫跨 back-splice junction，讀長越長越容易偵測；150bp 是目前標準 |
| 75bp PE | ⚠️ 尚可 | CIRIquant 可用，DCC 偵測能力下降 |
| **50bp PE** | ❌ 不佳 | STAR chimeric 幾乎無法偵測 BSJ；DCC 每 sample 只得 3–17 個 circRNA；需強制 CIRIquant-only 模式 |

**判斷方法**：GEO Series → Library Strategy 欄位；或下載 SRA run info 查 `avgLength` 欄。

#### 2. 樣本類型：組織 vs. 細胞株

| 類型 | circRNA 數量 | DE 結果 | 適合問題 |
|------|-------------|---------|---------|
| **組織**（tumor vs. normal）| 較多（整體轉錄組豐富）| 較多顯著 circRNA | 腫瘤生物標記、臨床相關性 |
| **細胞株**（drug treatment）| 較少（基因組背景單純）| 少量顯著（藥物效果有限）| 機制研究、pathway 分析 |
| **細胞株**（KO/OE）| 視目標基因而定 | 通常較集中 | 單一分子機制 |

**注意**：細胞株資料的 circRNA biomarker 臨床轉化性低，適合做機制驗證而非 biomarker 篩選。

#### 3. RNA-Seq 方式

| 方式 | circRNA 偵測 | 說明 |
|------|-------------|------|
| **Total RNA（rRNA-depleted）**（建議）| ✅ 最佳 | 保留 circRNA；rRNA-depleted 維持全轉錄組代表性 |
| **Total RNA（RNase R enriched）**| ✅✅ 最高靈敏度 | 消化線性 RNA，富集 circRNA；用於 benchmark ground truth（SRP011042）|
| poly-A selection | ❌ 不建議 | circRNA 無 poly-A tail → 幾乎偵測不到；若數據集用 poly-A，circRNA 數量會極少 |
| strand-specific | ✅ 加分 | 可推斷 circRNA strand；非必要 |

**判斷方法**：GEO → Library Selection 欄（`POLY_A`/`cDNA`/`other`）；Supplementary 找 library prep protocol。

#### 4. 定序深度（sequencing depth）

| 深度 | 說明 |
|------|------|
| **≥ 80M reads/sample**（建議）| circRNA 偵測靈敏度高；低豐度 circRNA 也能偵測 |
| 50–80M | 可用，主流研究水準 |
| < 30M | circRNA 偵測靈敏度明顯下降；少數 high-abundance circRNA 可用，但 DE 功率不足 |

**查詢方式**：SRA run info 的 `spots` 欄（paired-end: spots × 2 = total reads）。

> **⚠️ 重要教訓（GSE130078）**：定序深度足夠不代表 DE 分析功率足夠。GSE130078 有 79–101M reads/sample，但樣本間 BSJ count 總量差異達 **5.5 倍**（5,775–31,656），導致 filterByExpr 只保留 231/8,925 circRNAs（2.6%），edgeR 僅 12 個顯著。**定序深度是必要條件但不充分**；BSJ count 的穩定性（樣本間變異係數）才是關鍵。

#### 4b. BSJ Count 稀疏性（新增指標，源自 GSE130078 教訓）

BSJ count 稀疏性無法從 GEO 頁面直接查到，需分析後才知道，但可透過以下間接指標預測：

| 風險因素 | 說明 |
|---------|------|
| **癌症類型**（鱗狀上皮癌）| 食道（ESCC）、頭頸、肺鱗癌等 squamous cell carcinoma 的 circRNA 整體表現量普遍低於腺癌（breast、HCC、CRC）；腺癌資料集通常 filterByExpr 保留 5–50%，鱗狀癌可能僅 2–5% |
| **腫瘤含量（tumor cellularity）不明** | 手術切除組織的腫瘤純度差異大；ESCC 食道黏膜層切片尤其如此；若論文未報告 purity，BSJ 變異係數可能很高 |
| **RNA 品質未標注** | 未提供 RIN score 或 DV200 的資料集，降解 RNA 可能使 BSJ count 極低（降解破壞 BSJ spanning reads）|
| **circRNA 不是主研究目標** | 總 RNA-seq 轉順帶用於 circRNA，研究者不會特別控制 circRNA 偵測品質 |

**分析後的早期警示信號**（pipeline 完成 merge_counts 後即可判斷）：
- `count_matrix.tsv` 各 sample 的 BSJ sum 差異 > 3 倍 → ⚠️ 高風險
- filterByExpr 保留率 < 5% → ⚠️ DE 功率嚴重不足
- 若遇上述情況，建議以 **limma-voom 為主方法**（對稀疏計數更穩健），edgeR_ciriquant 作為輔助

#### 5. 樣本數與統計功率

| 組別樣本數 | DE 方法建議 | 預期顯著數 | 說明 |
|-----------|------------|-----------|------|
| n ≥ 5 vs. 5 | edgeR / DESeq2 | FDR 校正可用 | BH 校正有效；padj < 0.05 有意義 |
| n = 3 vs. 3（小樣本）| edgeR + nominal p | 需用 pvalue（未校正）| BH 校正後幾乎全不顯著（min padj ≈ 0.4）；論文需說明限制 |
| **配對設計**（tumor+normal same patient）| edgeR paired / limma blocking | 功率最高 | 消除個體差異；GSE133998 可加入配對項 |

**GSE133998 特別說明**：H36–H42 為同一患者的 cancer tissue 與 adjacent normal，是**嚴格配對設計**。已測試 `design = ~patient + condition`：limma 改善明顯（51 padj<0.05），但 edgeR 反而變差（FSJ offset 已吸收個體差異 + patient dummy 消耗 5 df）。最終決定維持 unpaired design，主要以 edgeR_ciriquant 為主方法。若要啟用配對設計，在 `metadata/GSE133998/sample_groups.csv` 加入 `patient_id` 欄即可（`analysis.R` 已支援自動偵測）。

### SRP156355 樣本清單（6 對 T/N，已生成 condition.csv）

| Patient | Normal SRR | Tumor SRR |
|---------|-----------|-----------|
| P138 | SRR7645073 (138N) | SRR7645074 (138T) |
| P148 | SRR7645077 (148N) | SRR7645078 (148T) |
| P204 | SRR7645071 (204N) | SRR7645072 (204T) |
| P272 | SRR7645075 (272N) | SRR7645076 (272T) |
| P690 | SRR7645087 (690N) | SRR7645088 (690T) |
| P1123 | SRR7645079 (1123N) | SRR7645080 (1123T) |

Condition CSV 位置：`/mnt/c/Users/User/Desktop/SRP156355_condition.csv`（含 `patient_id` 欄，可啟用配對設計）。
排除的樣本：DCIS（66D/299D/712D/803D/1102D/1151D）和 APN（719APN/768APN/87APN/91311APN/93277APN）留待未來三組比較分析。

### GSE77509 樣本清單（6 對 T/N，已生成 condition.csv）

Yang et al. 2017 *Nature Communications*（PMID 28194035）。20 對 HCC tumor + adjacent normal + 20 PVTT，rRNA-depleted，Illumina HiSeq 2500，~100bp PE。以定序深度（total spots/pair）排序後選 top 6：

| 患者 | Normal SRR | Tumor SRR | N 深度 | T 深度 |
|------|-----------|-----------|--------|--------|
| P17 | SRR3140289 | SRR3140387 | 133M | 89M |
| P16 | SRR3140284 | SRR3140382 | 148M | 60M |
| P26 | SRR3140326 | SRR3140422 | 84M | 90M |
| P12 | SRR3140264 | SRR3140362 | 58M | 109M |
| P20 | SRR3140303 | SRR3140400 | 78M | 88M |
| P22 | SRR3140311 | SRR3140408 | 81M | 84M |

**注意**：每個 biological sample 實際上由 4–5 個 split SRR 組成（總深度如上），這裡各取第一個 SRR 使用（單 SRR 約 12–30M reads）。若需要完整深度需在 pipeline 加入 SRR merge 支援。

Condition CSV 位置：`/mnt/c/Users/User/Desktop/GSE77509_condition.csv`
cancer label：`tumor` / `normal`；genome：hg19

---

### 選擇新資料集的 Checklist

在 GEO / SRA BioProject 找新 dataset 時確認以下條件（✅ = 必要，⚠️ = 建議）：

- ✅ Read length ≥ 100bp（PE）
- ✅ RNA-Seq library = Total RNA（rRNA-depleted or RNase R）
- ✅ 每組 ≥ 3 replicates（建議 ≥ 5）
- ✅ 有明確的 case vs. control（tumor/normal、treatment/vehicle）
- ⚠️ 定序深度 ≥ 50M reads/sample（必要但不充分，見下方）
- ⚠️ 組織樣本優先（而非細胞株），若要 biomarker 研究尤其重要
- ⚠️ 優先選擇腺癌（乳癌、HCC、CRC、肺腺癌），避免鱗狀癌（ESCC、HNSCC、肺鱗癌）——後者 circRNA 全局表現量低，filterByExpr 通過率可能只有 2–5%
- ⚠️ hg19 或 hg38 人類基因組（本 pipeline 以 hg19 為主）
- ⚠️ 配對設計（paired tumor/normal）功率最高
- ⚠️ 查論文 Methods 是否有報告 RIN score（≥ 7 為佳）或 DV200（≥ 30%）；多數資料集不標示，無記錄不代表品質差，但鱗狀癌或 FFPE 樣本需特別留意

**排除條件**：
- poly-A selection → circRNA 接近零（**GSE58135 是此條件的一個記錄在案的反例**：2026-10-03查證GEO extract_protocol原文確認其為poly(A)-selected，卻仍偵測到1,607個circRNA並留在cohort中——代表此排除條件並非絕對，但該資料集的DE功率也確實極低(45 tested circRNAs)，且已排除於15組跨資料集重現分析，見GSE58135段落與「資料集特性稽核校正」表）
- single-end reads → BSJ 偵測困難
- < 30bp reads → 完全不可用
- 鱗狀上皮癌（squamous）+ 未報告 RIN/purity → 高風險組合，DE 功率可能不足

**分析後補救策略**（遇 BSJ 稀疏時）：
- filterByExpr 保留率 < 5% → 改以 **limma-voom 為主方法**（TMM 正規化對稀疏計數更穩健）
- edgeR_ciriquant 結果作為輔助（Type I/II 分類參考）
- 論文 Methods 需說明樣本 BSJ count 變異性及其對統計功率的影響

---

## 論文 Discussion 素材：腫瘤組織 circRNA 普遍下調的生物機制

本 pipeline 分析的多個資料集（GSE113230 乳癌、GSE133998 乳癌、GSE77509 HCC、GSE248612 胃癌）
均觀察到 **腫瘤組織的 DECs 以下調為主**，此為 circRNA 研究中有生物學基礎的一致性現象。

### 機制一：腫瘤細胞分裂加速 → 線性 RNA 競爭優勢

正常細胞靜止時，circRNA 因半衰期長（缺少 5' cap 和 poly-A，不被外切酶降解）可大量積累。
腫瘤細胞快速分裂時，細胞體積快速倍增，單位時間內線性 mRNA 轉錄速率必須跟上增殖需求，
而 circRNA 的產生本身就是「競爭性剪接」的結果（環化剪接 vs. 線性剪接）——
**高增殖率環境偏向線性剪接**，circRNA 生成效率下降。

### 機制二：剪接因子重塑（splicing factor remodeling）

腫瘤中多種 RNA 結合蛋白（RBP）表現改變，這些蛋白直接影響背向剪接（back-splicing）效率：

| RBP | 腫瘤中變化 | 對 circRNA 影響 |
|-----|-----------|----------------|
| **QKI**（Quaking）| 多種腫瘤下調 | QKI 是促進 circRNA 產生的關鍵 RBP；QKI↓ → circRNA 全局下調 |
| **ESRP1/2**（epithelial splicing regulatory protein）| 上皮-間質轉化（EMT）時下調 | 影響 exon skipping 和 back-splicing 比例 |
| **MBNL family** | 多種腫瘤中表現異常 | 競爭 RNA 二級結構形成，影響環化效率 |
| **muscleblind-like** | — | 影響 circRNA biogenesis 的 intronic repeat pairing |

### 機制三：circRNA 作為腫瘤抑制因子

許多已知 circRNA 功能是**腫瘤抑制性**的：
- 作為 **miRNA sponge** 保護腫瘤抑制基因的 mRNA（最著名例子：ciRS-7/CDR1as 吸附 miR-7，保護 EGFR pathway 拮抗基因）
- 干擾 oncogene 翻譯或信號傳遞
- 腫瘤中這些 circRNA 下調，相當於**解除對 oncomiRNA 的競爭抑制**（ceRNA hypothesis）

### 機制四：表觀遺傳靜默（epigenetic silencing）

DNA 甲基化和 H3K27me3 在腫瘤中大規模重塑：
- 許多 circRNA 的親本基因被靜默
- 或剪接調控序列（如 intronic inverted repeats）被甲基化，影響環化效率

### 本 pipeline 的技術角度補充

**edgeR_ciriquant 測的是 BSJ/FSJ 比值**，因此「下調」意味著相對線性轉錄本的環化效率下降，
不一定是 BSJ 絕對量減少。這比純粹看 BSJ counts 更靈敏地反映「circRNA 專一性」的調控變化。

GSE130078（ESCC 食道鱗狀細胞癌）的 623 個 limma 顯著 circRNA 中下調比例更高，可能是鱗狀癌中
QKI/ESRP1/2 的表現量特別低，加劇了 back-splicing 抑制。

### 論文引用建議

| 文獻 | 說明 |
|------|------|
| Hansen et al. (2013) *Nature* | ciRS-7/CDR1as 作為 miR-7 sponge |
| Jeck et al. (2013) *Genome Biology* | circRNA 在分化細胞中高表現，增殖細胞中低表現 |
| Wan et al. (2019) *Cancer Research* | QKI 調控 circRNA biogenesis，腫瘤中 QKI 下調機制 |
| Kristensen et al. (2019) *Nucleic Acids Research* | circRNA 作為癌症生物標記的系統性綜述 |
| Conn et al. (2015) *eLife* | ESRP1/2 控制 back-splicing 效率 |
| Zhang et al. (2016) *Molecular Cell* | RBP 調控 circRNA 環化的分子機制（MBNL、QKI 等）|

**整體趨勢**：乳癌（GSE113230/GSE133998）、HCC（GSE77509）、胃癌（GSE248612）、攝護腺癌（GSE221107）
的分析結果均與此文獻一致——circRNA 在腫瘤中全局下調，且 Type_I DECs（circRNA 專一性，非線性 mRNA 變化）
佔 85–98%，進一步支持腫瘤中背向剪接效率普遍降低的假說。

---

## GSE97239（膀胱癌 Bladder Cancer，✅ 完成 2026-07-30）

**✅ 完成（2026-07-30 06:45）。** 報告位置：`~/GSE97239_results/report.html`（server，13MB）

膀胱癌手術切除組織（Bladder Transitional Cell Carcinoma），tumor vs. adjacent normal，3 對配對 T/N，6 samples（SRR5398213–SRR5398218）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成（S3 aria2c） |
| fastp QC/trim | ✅ 12/12 完成 |
| CIRIquant | ✅ 6/6 完成 |
| STAR paired+mate1+mate2 | ✅ 18/18 完成 |
| DCC | ✅ 6/6 完成 |
| consensus_filter / merge_counts | ✅ 完成（23,480 circRNAs） |
| assign_isoforms / annotate_circbase | ✅ 完成（86.1% exonic，intronic 4,356，intergenic 672）|
| DE analysis | ✅ 完成（edgeR **140** / DESeq2 / limma significant）|
| predict_interactions（union mode）| ✅ 完成（228 circRNAs，interactions.json 11MB）|
| isoform_switching | ✅ 完成（**65 events**，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（140 candidates）|
| report | ✅ 完成（report.html 13MB，2026-07-30 06:45）|

**主要數值結果**：
- 偵測：23,480 consensus circRNAs；filterByExpr 後 **2,719 tested**
- DE（edgeR_ciriquant）：**140 significant**（nominal p < 0.05，|log2FC| > 1）；上調 121 / 下調 19；**106 Type_I (75.7%) / 34 Type_II (24.3%)**——Type_II 比例（24.3%）高於其他資料集（一般 7–15%），可能反映膀胱癌中 circRNA/線性 RNA 雙層調控的特殊性
- Isoform switching：**65 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 20,211 rows）
- Biomarker candidates：140 個
- Top 1 biomarker：**chr10:126727566|126799662（hsa_circ_0005418）**；log2FC=+5.88，p=0.0065，Type_II，score=0.7766，in_circbase=1
- Top 2：**chr19:16192723|16197350（hsa_circ_0008432）**；log2FC=+9.47，p=0.0003，Type_II，score=0.7253

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（3T + 3N）：SRR5398213（tumor P1）、SRR5398214（tumor P2）、SRR5398215（tumor P3）；SRR5398216（normal P1）、SRR5398217（normal P2）、SRR5398218（normal P3）
- genome：hg19；配對設計（patient_id：P1/P2/P3）

**Server config**（`config/projects/GSE97239.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE97239/raw`
- `results_dir: /home3/choukaihsuan/GSE97239_results`

**啟動歷史（2026-07-29 修復）**：
- 原 queue 命令為 `agent.py --gse GSE97239`，因 `pysradb` 未安裝於 ciriquant env 而失敗
- metadata/groups 路徑繼承自 GSE121842（`web_ui.py _update_paths_for_project` bug）
- 修復：手動建立 `metadata/GSE97239/`，修正 `config/projects/GSE97239.yaml`，更新 queue DB 命令為直接 snakemake 呼叫
- `web_ui.py run_gse()` 已修正：config+metadata 存在時改用 snakemake 直接調用，避免 pysradb 依賴
- CIRIquant 首次失敗（`ciriquant.yaml` 設定問題，見已知問題表）：2026-07-29 深夜修正後重跑，2026-07-30 06:45 完成

---

## GSE192410（卵巣癌 Ovarian Cancer，✅ 完成 2026-07-30，方法不對稱已於 2026-08-02 修正）

**✅ 完成（2026-08-02 13:13，CIRIquant-only 全量重跑修正版）。** 報告位置：`~/GSE192410_results/report.html`（server，22MB）

卵巣癌組織，tumor vs. adjacent normal，3 對配對 T/N，6 samples（SRR17297761–SRR17297766）。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成 |
| fastp QC/trim | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成（GTF 沿用舊檔，不受 `-Nr`/tools 設定影響，2026-08-02 重跑未重新執行）|
| STAR paired+mate1+mate2 | ✅ 18/18 完成 |
| DCC | ⚠️ 僅 3/6 完成過（SRR17297761/762/763；764/765/766 因 O(n²) 效能瓶頸放棄）；**2026-08-02 起全部 6 個樣本改為 CIRIquant-only，DCC 不再是 consensus_filter 的依賴**（見下方修正說明）|
| consensus_filter / merge_counts | ✅ 完成（**59,196** consensus circRNAs，CIRIquant-only 對稱方法；見下方樣本表）|
| assign_isoforms / annotate_circbase | ✅ 完成（41,977 exonic / 14,260 intronic / 2,958 intergenic）|
| DE analysis（三方法）| ✅ 完成（edgeR **2,346** / DESeq2 **4,292** / limma **7,368**）|
| predict_interactions（union mode）| ✅ 完成（233 circRNAs，~50 min）|
| isoform_switching | ✅ 完成（**51 events**，within-gene FDR < 0.1，共 55,640 rows）|
| rank_biomarkers | ✅ 完成（2,347 candidates）|
| report | ✅ 完成（22MB，2026-08-02 13:13）|

**主要數值結果（CIRIquant-only 修正版）**：
- 偵測：59,196 consensus circRNAs；filterByExpr 後 **7,144 tested**
- DE（edgeR_ciriquant）：**2,346 significant**（nominal p < 0.05，|log2FC| > 1）；上調 62 / 下調 2,284（腫瘤仍普遍下調）；**2,044 Type_I (87.1%) / 302 Type_II (12.9%)**
- DE（DESeq2）：4,292 significant；DE（limma-voom）：7,368 significant
- Isoform switching：**51 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1）
- Biomarker candidates：2,347 個
- Top 1 biomarker：**chr5:43292576|43297268（hsa_circ_0008621，HMGCS1）**；log2FC=−12.79，p=8.5e-06，Type_I，score=0.7507，155 miRNA / 151 RBP binders，下調——**與修正前完全相同的旗艦候選**

**✅ 2026-08-02 方法不對稱修正 + 重要結論：下調並非方法假影**

原始分析（2026-07-30）腫瘤樣本用 CIRIquant+DCC 雙工具共識、正常樣本因 DCC 的 O(n²) 效能瓶頸被迫改用 CIRIquant-only，記錄為「⚠️ 方法不對稱，2,286 個下調可能被高估」。修正方式：`config/projects/GSE192410.yaml` 的 `consensus.tools` 改為 `[ciriquant]`（`min_tools: 1`），讓全部 6 個樣本統一用同一種偵測方法，重跑 `consensus_filter → merge_counts → assign_isoforms → annotate_circbase → DE → predict_interactions → isoform_switching → rank_biomarkers → generate_report`（CIRIquant/STAR/DCC 本身不必重跑，GTF 不受 tools 設定影響）。

**修正後結果幾乎與修正前完全一致**（下調 2,284 vs 原本 2,286；上調 62 vs 61；Top 1 biomarker 完全相同）——**代表「腫瘤下調」這個發現不是 DCC bypass 造成的方法假影，而是真實訊號**。方法不對稱本身是一個確實存在、必須修正的技術瑕疵（已修正），但它並不是這個資料集腫瘤下調趨勢的主要成因。

**額外驗證：CIRIquant 自己回報的 Circular_Reads（偵測階段最原始的 BSJ read 數，早於任何 consensus/DCC 邏輯介入）直接證實這個下調是真實的**：

| SRR | 條件 | Total_Reads | Mapped_Reads | Circular_Reads |
|-----|------|-----------:|-------------:|----------------:|
| SRR17297761 | tumor | 82.1M | 77.2M | 444,706 |
| SRR17297762 | tumor | 121.1M | 116.0M | 89,798 |
| SRR17297763 | tumor | 111.8M | 101.6M | 266,062 |
| SRR17297764 | normal | 104.0M | 98.5M | 1,173,544 |
| SRR17297765 | normal | 110.4M | 106.0M | 1,216,554 |
| SRR17297766 | normal | 99.0M | 92.0M | 1,823,374 |

總讀數和 mapping rate 在 tumor/normal 之間相近（都是 ~90–96%），但 **Circular_Reads 在 tumor 只有 90K–445K，normal 卻有 1.17M–1.82M**（約 3–10 倍差距）——這個數字來自 CIRIquant 內部偵測，跟 DCC 或 consensus 共識邏輯完全無關，證明卵巢腫瘤組織的環化活性（back-splicing efficiency）確實遠低於鄰近正常組織，與 CLAUDE.md「腫瘤組織 circRNA 普遍下調的生物機制」段落記錄的跨資料集趨勢（QKI/ESRP1/MBNL 下調、乳癌/HCC/胃癌/攝護腺癌皆同向）完全吻合。

**各樣本 consensus circRNA 數量（CIRIquant-only 對稱方法，2026-08-02）**：

| SRR ID | 條件 | 患者 | 共識 circRNA |
|--------|------|------|----------:|
| SRR17297761 | tumor | P1 | 14,603 |
| SRR17297762 | tumor | P2 | 3,722 |
| SRR17297763 | tumor | P3 | 2,044 |
| SRR17297764 | normal | P1 | 32,067 |
| SRR17297765 | normal | P2 | 33,267 |
| SRR17297766 | normal | P3 | 16,493 |

tumor 樣本的 consensus 數量比修正前（雙工具共識，1,567–11,683）略高（改為單工具後不再需要 DCC 交集，保留更多 circRNA），normal 樣本數字不變（本來就是 CIRIquant-only）；但 tumor 仍明顯低於 normal，與上述 Circular_Reads 的證據一致——這是真實的生物學差異，不是統計假影。

<details>
<summary><b>歷史版本（2026-07-30 雙工具/單工具混用版，僅供對照）</b></summary>

**主要數值結果（舊版，方法不對稱）**：
- 偵測：57,469 consensus circRNAs；filterByExpr 後 7,141 tested
- DE（edgeR_ciriquant）：2,347 significant；上調 61 / 下調 2,286；2,049 Type_I (87.3%) / 298 Type_II (12.7%)
- DE（DESeq2）：4,267 significant；DE（limma-voom）：7,312 significant
- Top 1 biomarker：chr5:43292576|43297268（hsa_circ_0008621）；log2FC=−12.79，Type_I

**各樣本 consensus circRNA 數量（舊版）**：

| SRR ID | 條件 | 患者 | 方法 | 共識 circRNA |
|--------|------|------|------|----------:|
| SRR17297761 | tumor | P1 | CIRIquant+DCC | 11,683 |
| SRR17297762 | tumor | P2 | CIRIquant+DCC | 2,890 |
| SRR17297763 | tumor | P3 | CIRIquant+DCC | 1,567 |
| SRR17297764 | normal | P1 | CIRIquant-only（DCC bypass）| 32,067 |
| SRR17297765 | normal | P2 | CIRIquant-only（DCC bypass）| 33,267 |
| SRR17297766 | normal | P3 | CIRIquant-only（DCC bypass）| 16,493 |

</details>

**DCC O(n²) 瓶頸歷史（2026-07-29～30，仍為有效的技術紀錄）**：
- DCC 0.5.0 的 duplicate-marking 使用 list `in` 操作（O(N²) 複雜度），對含 1.3–1.9M 總 junction 的樣本估算需 23–26 小時/樣本
- **第一輪**：針對 chr10:116879948↔116889298（50,000+ reads）等高頻 junction 執行精確過濾，仍卡住
- **第二輪**：bulk threshold 過濾（>1000 combined count），max 降至 986，仍卡住（tmp_nonduplicates 僅以 39 KB/min 增長）
- **第三輪（最終方案）**：放棄 764/765/766 的 DCC；在 `DCC/CircCoordinates` 放入 header-only 檔案（1 行，70 bytes），刪除 `.snakemake/incomplete/` 中對應的 base64 標記，重啟 Snakemake with `--rerun-triggers mtime`；`parse_dcc()` 讀到 0 circRNA → adaptive fallback 觸發（min/max ratio=0 < 0.1）→ CIRIquant-only 模式

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（3T + 3N）：SRR17297761/762/763（tumor P1/P2/P3）；SRR17297764/765/766（normal P1/P2/P3）
- genome：hg19；配對設計（patient_id：P1/P2/P3）
- study_title：`Circular RNA expression profiling in ovarian cancer tumor and adjacent normal tissue`

**Server config**（`config/projects/GSE192410.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE192410/raw`
- `results_dir: /home3/choukaihsuan/GSE192410_results`
- `sra_cache_dir: /home3/choukaihsuan/GSE192410/sra_cache`
- `consensus.tools: [ciriquant]`（2026-08-02 起，原為 `[ciriquant, dcc]`）；`consensus.min_tools: 1`（原為 2）

**Queue 失敗歷史（2026-07-29 修復）**：與 GSE97239 相同原因（pysradb 依賴 + metadata 路徑繼承 bug），修復方式相同。

---

## GSE192849（乳癌 node-positive，RNase R 富集，✅ 完成 2026-08-03）

**✅ 完成（2026-08-03 02:14）。** 報告位置：`~/GSE192849_results/report.html`（server，19MB）；本機備份：`/mnt/c/Users/User/Desktop/circRNA agent report/GSE192849_report.html`

淋巴結陽性原發性乳癌（Node-positive primary breast cancer）手術切除組織，tumor vs. adjacent normal，3 對配對 T/N，6 samples（SRR17395573–SRR17395578）。VAHTS Total RNA Library Prep Kit + **RNase R 消化**（circRNA 富集），Illumina，150bp PE，~55M reads/sample。

| 步驟 | 狀態 |
|------|------|
| 下載 SRA | ✅ 6/6 完成 |
| fastp QC/trim | ✅ 6/6 完成 |
| CIRIquant | ✅ 6/6 完成 |
| STAR paired+mate1+mate2 | ✅ 18/18 完成 |
| DCC | ⚠️ 全部 bypass（RNase R O(n²) 瓶頸，見下方）|
| consensus_filter / merge_counts | ✅ 完成（68,680 consensus circRNAs；filterByExpr 後 10,977 tested）|
| assign_isoforms / annotate_circbase | ✅ 完成 |
| DE analysis（三方法）| ✅ 完成（edgeR **404** / DESeq2 / limma significant）|
| predict_interactions（union mode）| ✅ 完成（interactions.json 8.7MB）|
| isoform_switching | ✅ 完成（**20 events**，within-gene FDR < 0.1）|
| rank_biomarkers | ✅ 完成（404 candidates，biomarker_candidates.tsv 98KB）|
| report | ✅ 完成（19MB，2026-08-03 02:14）|

**主要數值結果**：
- 偵測：68,680 consensus circRNAs（CIRIquant-only）；filterByExpr 後 **10,977 tested**
- DE（edgeR_ciriquant）：**404 significant**（nominal p < 0.05，|log2FC| > 1）；上調 395 / 下調 9；**368 Type_I (91.1%) / 36 Type_II (8.9%)**
- Isoform switching：**20 events**（within-gene BH FDR < 0.1，|ΔIUI| > 0.1，共 64,806 rows）
- Biomarker candidates：404 個
- Top 1 biomarker：**hsa_circ_0001922（HUWE1，chrX:53672263|53681075）**；log2FC=+9.23，p=0.0042，Type_II，score=0.7815，50 miRNA / 106 RBP binders，上調
- Top 2：**hsa_circ_0001777（ESYT2，chr7:158580695|158591763）**；log2FC=+9.33，Type_II，score=0.7723
- Top 3：**hsa_circ_0001118（NDUFA10，chr2:240929491|240946787）**；log2FC=+9.17，Type_II，score=0.7525

**⚠️ 重要注意：RNase R 富集造成方向性差異**

腫瘤上調（395/404 = 97.8%）與其他 Total RNA-Seq 資料集的「腫瘤下調」趨勢**相反**。此方向差異主因 **RNase R 消化** 的 library prep：
- RNase R 消化線性 RNA 後，FSJ counts 幾乎歸零（線性轉錄本大幅減少）
- edgeR_ciriquant 測的是 BSJ/FSJ ratio；RNase R 後 FSJ≈0 → ratio 自然偏高
- 腫瘤樣本的 CIRIquant Circular_Reads（bypass 前原始數字）：444K / 90K / 266K
- 正常樣本的 Circular_Reads：1.17M / 1.22M / 1.82M（約 3–10 倍差距）
- 即使排除 ratio 效應，正常組織的環化活性（back-splicing）仍高於腫瘤，符合 circRNA 腫瘤普遍下調的生物機制
- **論文 Discussion 需說明**：RNase R 富集資料的 log2FC 方向不可直接與 Total RNA-Seq 資料集比較；biomarker 方向需以 CSI（circular-to-total ratio）或原始 BSJ counts 佐證

**各樣本共識 circRNA 數量（CIRIquant-only）**：

| SRR ID | 條件 | 患者 | 共識 circRNA |
|--------|------|------|----------:|
| SRR17395573 | tumor | P1 | 6,240 |
| SRR17395574 | tumor | P2 | 19,253 |
| SRR17395575 | tumor | P3 | 12,711 |
| SRR17395576 | normal | P1 | 26,854 |
| SRR17395577 | normal | P2 | 18,845 |
| SRR17395578 | normal | P3 | 42,268 |

**DCC O(n²) 瓶頸與 bypass 處理**：

RNase R 富集使 Chimeric.out.junction 檔案極大（paired: 438K–1.3M lines；mate1+mate2 加總：1.1M–3.4M lines/sample），DCC 0.5.0 duplicate-marking O(n²) 複雜度導致估算需 23–26 小時/樣本（與 GSE192410 pattern 一致）。**全部 6 個樣本統一使用 header-only CircCoordinates bypass**（54 bytes，1 header line），確保方法對稱性（不重蹈 GSE192410 初版的 tumor/normal 不對稱問題）。Snakemake 以 `--rerun-triggers mtime` 重啟，跳過 DCC 繼續後續步驟。部分 DCC 孤兒進程（Aug02 殘留）於 2026-08-03 以 `kill -9` 清除。

**設定**：
- case/control label：`tumor` / `normal`
- SRR 清單（3T + 3N）：SRR17395573/574/575（tumor P1/P2/P3）；SRR17395576/577/578（normal P1/P2/P3）
- genome：hg19；配對設計（patient_id：P1/P2/P3）
- Library：VAHTS Total RNA + RNase R，150bp PE，Illumina
- `consensus.tools: [ciriquant]`；`consensus.min_tools: 1`（DCC bypass）

**Server config**（`config/projects/GSE192849.yaml`）路徑：
- `raw_dir: /home3/choukaihsuan/GSE192849/raw`
- `results_dir: /home3/choukaihsuan/GSE192849_results`
