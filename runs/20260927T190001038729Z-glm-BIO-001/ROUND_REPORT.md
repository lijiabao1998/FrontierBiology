# BIO-001 r1：metadata 審計與未完成的洩漏驗證

目前科學判定：**E2_INCONCLUSIVE_SIBLING_FREE_CLEAN_NOT_CONSTRUCTIBLE**。
BIO-001 保持 OPEN。原始 acceptance 與檢索紀錄保留；本次只修正已有證據，
沒有新增研究輪次或補造歷史檢索。

## 已有資料

Norman GSE133344 metadata 共 111,445 cells：57,831 single、41,759 dual、
11,855 control；105 個單擾動，每擾動 113–1,960 cells（median 495）。
表達矩陣 1.1GB 未下載或處理，沒有真實資料泛化評估。

## E2 負結果

同一 held-gene test cohort 為770列。LEAKY Pearson=0.8776；真正移除functional
siblings後，CLEAN train=0，沒有可用的sibling-free估計器。因此CLEAN Pearson
與group gap都是null，原先洩漏PASS撤回。r2改成模型消融也不能補救這個負結果。

## donor post-hoc 診斷與覆核修正

主實作的test固定為groups 6/7、donors 0/1/2，共1,105列。LEAKY train為其他
群組所有donors的4,504列；CLEAN train再限donor 3，共1,100列。兩者測試列
相同且無train/test重疊。因test groups未出現在train，LEAKY的估計公式
`donor_mean[d] + group_mean.get(group, global_mean) - global_mean`簡化成
`donor_mean[d]`；CLEAN為其train global mean。Pearson為0.1616與0.0。
這同時改變exposure與估計器，只是post-hoc描述，不能當作隔離出的洩漏效應。

舊v4 verifier用了另一個切分：non-held genes為train（5,230列），所有groups的
前三個donors為test（4,509列），其中3,929列同時在train/test。0.905因此不是同一
estimand，也不是已核實的獨立重現。舊MISMATCH JSON保存在
`results/r1/history/independent_leakage_verification.v4.json`。修正版不import主程式，
按上述row選擇與簡化公式獨立重建，並由空CLEAN train直接推導E2 verdict，
不從待驗證JSON的PASS/status字段推導答案。驗證結果與有界回歸紀錄見修復報告。

## r2 狀態及限制

r2過往未凍結0.2門檻、未證明本輪四路新檢索、沒有r2獨立verifier。已撤回
admission、完成與PASS標籤，僅保留探索性未驗證的0.8776對0.0數值。原始紀錄
另存historical_record_at_9590a99.json，沒有補做檢索冒充當時已完成。

## 重播與歷史

從repo根目錄執行：

```text
python -B problems/BIO-001/experiments/bio001_r1/bio001_split_evaluator.py
python -B problems/BIO-001/experiments/bio001_r1/independent_check.py
python -B problems/BIO-001/experiments/bio001_r1/independent_leakage_verifier.py
python -B problems/BIO-001/experiments/bio001_r2_ablation/eval_model_ablation.py
python -B problems/BIO-001/experiments/bio001_r1/tests_convergence.py
python -B problems/BIO-001/experiments/bio001_r1/verify_manifest.py --git-ref HEAD
```

目前結果、原始輸入與校驗器均受manifest涵蓋；產物不含即時時間戳，統一LF。
退出0可表示成功保存負結果，不表示科學假說通過。舊report另存
HISTORICAL_REPORT_at_9590a99.md；只有那份歷史檔內保留撤回的PASS與舊數字。
尚待獨立review及owner合併判斷；新研究必須另開有效預登記輪次。
