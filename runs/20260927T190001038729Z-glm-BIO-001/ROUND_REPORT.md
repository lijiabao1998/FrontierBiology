# ROUND REPORT｜BIO-001｜20260927T190001038729Z-glm-BIO-001

- branch：`glm/BIO-001-baseline-r1`（base `e0c162b2705760eee8cc22c1400da2370984fac4`）；pin `07d2b130`
- verdict：**NO_RESOLUTION_FOUND**（leakage/泛化為公認未解：Systema、PerturbVAE、baseline-gap preprint 一致）
- **REMEDIATION 2026-09-28（Codex review）**：下方初版數字 81,542/630/29,273/103 已被修正取代——正確值 57,831 single / 41,759 dual / 11,855 control、105 singles；初版宣稱僅適用於錯誤 parser，不應再引用。
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

## Remediation（2026-09-28，Codex review 回應）
- **P1 parser**：guide_identity 命名實為「token 過濾」語意（gene 名可含底線、guide 為 NegCtrl*/數字）；初版 rsplit 解析把 SET_KLF1__SET_KLF1 誤判單擾動、NegCtrl0_BAK1__... 誤判 control。修正後 class counts 大幅翻正（dual 630→41,759），並與 Norman 2019 雙擾動設計之文獻認知一致。回歸測試＝independent_check.py 內 4 條 spot assertions。
- **P1 donor 對照**：新增 donor_holdout（donor 3 僅出現在 test）之 donor-clean split；donor_leak_gap 改為 vs donor-clean（0.1869 ≥ 0.1 ✓）——原比較只反映特徵開關。
- **P1 獨立覆核**：`independent_check.py`（regex+memoisation 結構相異實作）6/6 checks 吻合（cells_total、class_counts、105 singles、median、fold 決定性、spot identities）；round.json 補 independent_check 紀錄。
- **P2 hashes**：重算並通過 `sha256sum -c`。
- 收窄宣稱：初版「103 個單擾動／630 dual cells」等數字撤回；dataset_map.md 已註記修正。

## Remediation v2（Convergence Wave，Codex 二審回應）
- P1 same-T donor 對照：CLEAN（donor-blind）vs LEAKY（donor-label 通道開）共用完全相同之 test rows T 與訓練列；same_test_rows_invariant=true；gap 0.1869 於相同列上定義；donor-3-only 子集兩 predictor 皆常數（Pearson 0）之事實如實記錄。
- P1 正名：僅 E2（group gap ≥0.2）為 preregistered；D1/D3/D4 改標 post_hoc_diagnostics，不回寫歷史。
- P1 獨立洩漏驗證：independent_leakage_verifier.py（不 import 主 evaluator，自行構造 split/預測/Pearson）5/5 吻合主 verdict——V1 group gap、V2 donor gap、V3 same-rows、V4/V5 frozen/leaky 重現。
- P2 manifest：verify_manifest.py + 全輸入（.gz）與輸出 hash，最後生成，全綠。
- P2 round.json summary 103→105 singles 已更正。
