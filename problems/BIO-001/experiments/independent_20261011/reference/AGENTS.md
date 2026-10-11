# FrontierBiology agent 入口

先讀 README、STATUS、VALIDATION 與固定治理：
https://github.com/lijiabao1998/FrontierLab-Governance/tree/07d2b13051b83215182e411e1612f92f1912d8fb
其中 AGENTS、RESEARCH_PROTOCOL、EVIDENCE_POLICY、SAFETY 全部適用。

所有 agent 本次初始化後走 `<agent>/BIO-xxx-<topic>`，不直接寫 main、不自合；owner 明確批准後才合。fetch 記 base SHA，一工作目錄一寫入者；領題先查 PR／runs，避免重做。

每輪 start → 本輪四路聯網檢索並讀原文 → 凍結資料切分、指標、成功條件及預算 → admit → 重現 → 探索 → Verifier/Skeptic獨立檢查 → 記失敗 → PR。無網路不能假造檢索。確認外部同範圍完成就標 COMPLETED_EXTERNAL，停止新發現搜尋；未核實聲稱先 PAUSED。

禁止把相關性寫成因果、預測寫成干預有效、simulation寫成experiment、同資料重新切分寫成獨立replication。方法作者自己重跑不是獨立證據；學術指標與生物解釋分開。

只用公開合法去識別及低風險資料；病原／毒素功能優化、未審核濕實驗、個資、重新識別、臨床應用不在自動研究權限內。原文含高風險例子不代表授權本庫研究該例子。

研究檔案放 problems/<ID>/{experiments,proofs,results}/；每輪 runs/<round>/round.json。保存版本、seed、資料來源license/hash、指標及負結果；未跑寫 NOT_RUN。預設30分鐘／0美元／100次，到限就停。首輪只做 BIO-001 基線。
