# 狀態｜2026-09-29

BIO-001仍OPEN。已有r1研究紀錄與r2探索產物；10張題卡中只有BIO-001實作了
metadata/evaluator工具。本次是已有證據修復，不新增研究輪次。

- r1：Norman metadata 111,445 cells／105 singles；E2 CLEAN train為0，
  **INCONCLUSIVE**，沒有洩漏PASS。donor 0.1616僅為post-hoc exposure與估計器比較。
- 舊verifier的0.905來自不同cohort及train/test重疊；歷史MISMATCH保存，
  修正版按主實作原有estimand獨立重建。具體檢查及review狀態見修復報告。
- r2：admission與PASS已撤回，**EXPLORATORY_UNVERIFIED**；缺有效事前門檻、
  本輪四路新檢索與獨立驗證。保留描述數值，不冒充完成研究。
- 表達矩陣1.1GB未處理；無真實資料泛化、機制、濕實驗或臨床結果。

權威紀錄：runs/20260927T190001038729Z-glm-BIO-001/ROUND_REPORT.md，以及
runs/20260929-gpt-BIO-001-convergence/REPAIR_REPORT.md。records CI只驗格式與紀錄，
不替科研結果背書。新研究須重新檢索、凍結、admit；本修復不補造過往搜尋。
