# LITERATURE_MAP｜BIO-001｜20260927T190001038729Z-glm-BIO-001

檢索窗：2026-09-27T19:00Z–19:06Z。引擎：zcode-websearch／zcode-webfetch／curl(GEO)。

## 檢索紀錄（6 query）
| # | category | query | outcome |
|---|---|---|---|
| 1 | general | `"single-cell perturbation response prediction" benchmark generalization Nature Methods evaluation` | 定位 Wei et al. Nat Methods 2025（題卡 primary）；發現 PerturbVAE leakage-aware benchmark |
| 2 | discipline | `Norman 2019 perturb-seq GSE133344 data download size format counts csv` | GEO 定位；dense CSV 說法待實測 |
| 3 | discipline | `https://ftp.ncbi.nlm.nih.gov/geo/.../GSE133344/suppl/`（curl 實測） | identities 1.9MB／genes 259K／matrix 1.1GB；前兩檔下載審計 |
| 4 | solution | `https://github.com/brhanufen/perturbvae`（webfetch） | split 方法論 A/B/C（BLAKE2b）；Norman A→B Pearson 0.522→0.257；MIT |
| 5 | solution | `perturbation response prediction "simple baselines" hard to beat mean baseline ...` | 「baseline gap」2026 preprint；mean baseline 難擊敗為公認議題 |
| 6 | criticism | `"perturbation response" single-cell prediction benchmark criticism flawed leakage 2025` | **Systema**（Nat Biotechnol 44:1050–1059, 2025/26）：benchmarks 被 systematic variation 灌水 |

## 圖譜
- **rigorously/computationally known**：Norman 2019 資料結構（本輪實測：111,445 cells、105 singles [corrected]、41,759 duals [corrected]、126–3,824 cells/perturbation）。
- **conjectured/open**：跨細胞情境與未見擾動泛化——open（Wei 2025、PerturbVAE A→B gap、Systema、baseline-gap preprint 一致指向）。
- **disputed**：多數已發表模型的「泛化」宣稱受 leakage/systematic-variation 質疑（Systema；Camillo 對 virtual-cell 的 mean-baseline 批評）。
- **superseded**：無（題卡 known_result 仍準確）。
- **unresolved**：BIO-001 主問題本身。

## 限制
Nature primary 被 cookie 牆擋（僅 secondary 確認存在與卷期）；PerturbVAE 數值取自其 repo 報告層級未重算；6 query 有界；「無結果」非未解證明。
