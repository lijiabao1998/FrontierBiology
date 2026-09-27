# dataset_map｜BIO-001 r1｜2026-09-28（UTC 09-27 檢索）

## 已下載並審計（真實資料，存 experiments/lit_data/）

| 檔案 | 大小（實測） | 內容 | 用途 |
|---|---|---|---|
| `GSE133344_filtered_cell_identities.csv.gz` | 1,956,905 B | 每細胞：barcode、guide_identity、read/UMI count、coverage、gemgroup | **本輪已解析**：111,445 cells；81,542 single / 29,273 control / 630 dual；103 個單擾動；每擾動細胞數 min 126 / median 463 / max 3,824 |
| `GSE133344_filtered_genes.tsv.gz` | 264,791 B | ENSG↔symbol 對照 | 基因索引（下輪表達矩陣對照用） |

guide_identity 慣例（由實資料確認）：`<GENE>_<guide>__<GENE>_<guide>`；`NegCtrl*` 為對照 guide；兩半皆 NegCtrl = control cell。

## 可得但本輪未下載（DATA_TRACTABILITY_BLOCKED）

| 檔案 | 大小（GEO 目錄實測） | 阻碍 | 解除條件 |
|---|---|---|---|
| `GSE133344_filtered_matrix.mtx.gz` | 1.1 GB（壓縮） | 純 stdlib 串流解析 ~數 GB 文本/50M+ 行，超出本輪 90 min 時間盒；下游運算無 numpy | 下輪：`pip install --target .deps numpy scipy anndata`（裝於 repo 內）或實作串流取樣；MTX 為 gene-major MatrixMarket |
| `GSE133344_RAW.tar` | 10 GB | 不需要 | — |

## 資料授權與使用條件

- Norman et al. 2019（Science 2019）Perturb-seq（K562 細胞系，非人體/臨床樣本）：GEO 公開發佈；論文屬细胞系 in vitro 資料。使用條件：引用原始出處；無個資（細胞 barcode 非人類識別）。SAFETY.md 合規：公開、低風險、非臨床。
- PerturbVAE（brhanufen/perturbvae）：MIT；其 split 方法論（BLAKE2b 決定性、function-grouped holdout）被本輪**獨立重實作**（非複製程式碼）。
- Wei et al. Nat Methods 2025（題卡 primary）：cookie 牆未取得全文；引用條件待讀全文後補審。

## 下一步資料行動（精確）

1. `pip install --target FrontierBiology/problems/BIO-001/experiments/.deps numpy scipy`（僅在 repo 內）。
2. 串流解析 mtx.gz 的目標基因列（103 個單擾動 + control 的 ~9k 基因面板），或直接全量 sparse 載入。
3. 用本輪已交付的 `split_validator`/evaluator（`bio001_split_evaluator.py` 介面不變）跑真實 Norman：frozen（gene-function-grouped，需 GO/pathway 資源清單——下輪另查授權）vs random split 的 leakage gap 實測。
