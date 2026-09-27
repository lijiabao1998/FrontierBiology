# ROUND REPORT｜BIO-001｜20260927T190001038729Z-glm-BIO-001

- branch：`glm/BIO-001-baseline-r1`（base `e0c162b2705760eee8cc22c1400da2370984fac4`）；pin `07d2b130`
- verdict：**NO_RESOLUTION_FOUND**（leakage/泛化為公認未解：Systema、PerturbVAE、baseline-gap preprint 一致）
- 本輪判定：**LEAKAGE_DEMO_PASS**（4/4 預登記 checks）

## 交付（全部 stdlib、確定性）
1. **真實資料審計**（dataset_map.md＋results/r1/bio001_r1_results.json）：Norman 2019 GSE133344 metadata 實測——111,445 cells（81,542 single／29,273 control／630 dual）、103 個單擾動、每擾動 126–3,824 cells（median 463）；BLAKE2b 5-fold 分配已生成。
2. **split validator＋evaluator**（`bio001_split_evaluator.py`）：frozen（function-grouped）/ leaky（同組鄰居留在 train）/ donorleak 三模式，group-mean 基線、Pearson/MSE。
3. **洩漏負控制（合成 ground truth，seed 42/12345）**：frozen Pearson **0.000**（留出組無任何同組訓練訊號——正確的零）vs leaky Pearson **0.878**（gap **0.878** ≥ 0.2 ✓）；MSE 0.794→0.115；donor 洩漏 gap **0.187** ≥ 0.1 ✓。這以可控實驗量化了 Systema/PerturbVAE 在真實 benchmark 報告的膨脹機制。
4. **DATA_TRACTABILITY_BLOCKED（誠實記錄）**：表達矩陣 1.1GB mtx.gz 未下載（stdlib 時間盒）；解除路徑寫明（numpy/scipy in-repo 或串流取樣）。**未偽造任何真實資料結果**。
5. 失敗日誌：合成生成器初版偏斜權重使留出組 0 細胞（ZeroDivision）→ 改均勻指派後通過；Nature abstract 遭 cookie 牆（記 limitation）。

## 沒做成／限制
- 真實表達矩陣未處理 → 未對真實 Norman 宣稱任何泛化評估。
- function-grouped split 之 pathway 資源（GO/Reactome 授權版本）未整合——下輪資源清單。
- PerturbVAE 數值未獨立重算（僅方法論引用）。

## 下一輪最小下一步
1. in-repo 安裝 numpy/scipy → 串流解析 mtx → 真實 Norman 上跑同一 evaluator（frozen vs random）。
2. 取得 pathway 分組資源（授權審計）→ 真 function-grouped split。
3. 讀 Wei 2025 全文（proxy/校園管道）對齊其 split 定義。
