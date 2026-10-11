# 每輪開始前檢索與研究協議

## 0. 先界定正在解什麼
讀題卡的 statement、known_result、open_gap、completion_criterion。記錄量詞、參數範圍、資料／物種／條件、假設及成功標準。論文解了特例、改善下界、更新榜單或只在某資料集成功，不能直接關閉更廣的父問題。

## 1. 四路新鮮檢索（每一輪、本輪開始後）
至少四個不同的實際查詢：
- general：問題全名、別名、作者及核心描述的廣域搜尋。
- discipline：數學 arXiv/math、作者頁與正式期刊；物理 arXiv/APS/合作組；生物 PubMed/bioRxiv/正式期刊；化學 ChemRxiv/ACS/RSC/IUCr/正式期刊。這是來源類型，不要求被封鎖平台必須可用。
- solution：問題名加 proof/solved/resolved/counterexample/exact/improved bound 或同義詞，搜尋最新年份與不限日期兩種視角。
- criticism：相關新論文加 correction/erratum/retraction/comment/failed replication；確認版本及適用範圍。

每個 query 記原字串、引擎／索引、查到什麼；至少閱讀兩個不同來源，其中至少一個原始論文、作者正式證書或原始合作組資料。優先原文，不用 AI 摘要互相背書。只看到摘要就標 abstract-only；無全文／無資料／未查到更正不能填成已核實。

時間記 ISO 8601 含時區。`admit` 要求搜尋在本輪開始後、24 小時內；新輪次即使在同一天也重新查。沒有網路或原文不足以判讀，記 BLOCKED；不得用模型記憶填搜尋紀錄。

**全網檢索是盡可能廣的、可重做的有界搜尋，不是宣稱完整掃描整個網路。沒有搜尋結果不是未解證明。** 寫明檢索遺漏、付費牆、未開啟來源及不確定性。

## 2. 檢索結論與立即處置
| 結果 | 操作 |
|---|---|
| NO_RESOLUTION_FOUND | 只允許按既定預算進行本輪，保留未解狀態的不確定性 |
| PARTIAL_PROGRESS | 更新已知結果；只研究剩餘明確缺口 |
| CLAIMED_RESOLVED | 停止把原題當全新發現目標；轉核查解答，不能直接標完成 |
| RESOLVED_EXTERNAL | 精確 scope 一致且獨立核實後，立即在工作分支標 COMPLETED_EXTERNAL／已由外部完成，附作者、原始來源、獨立核查、日期；停止重複探索並提狀態更新 PR |
| BLOCKED | 缺網路／來源／資料／授權，停止探索，說明下一步 |

resolved 紀錄必填 `problem_id, scope_statement, scope_match, basis, confirmed_by, independent_check, independently_verified, confirmed_at, sources`。至少一個 primary 加另一個可核對來源或證書；兩篇轉載不等於獨立核查。腳本檢查欄位存在，**scope 和解答真實性必須人工／獨立 reviewer 審核**。

外部完成與本庫完成分開：`COMPLETED_INTERNAL` 只在本庫精確命題經合適驗證、獨立覆核及最新文獻排重後由批准的合併流程設定。新結果是否新穎與是否正確分開記。預印本有爭議則 CLAIMED_RESOLVED；已解問題仍可開明確標示 REPLICATION 的子任務，不冒充新解。

## 3. 有界施工
先執行 `start`，填輪次紀錄再 `admit`。每輪上限 120 分鐘，預設 30 分鐘／0 美元／100 次；付費 API、GPU、長任務先取得業主授權。搜尋、重現、探索、審核各記耗用；不要讓協調吃掉所有研究時間。

每輪的 acceptance 在實驗前凍結。保存原始輸入、程式、環境、seed、超參數、資料切分、evaluator版本、結果與 hashes。紀錄可以提交失敗，但不可把失敗記成成功。先跑已知正例和必須失敗的反例；測試不會紅則 verifier 不可信。

## 4. 結束與合併
結果說清：做成什麼、沒做成什麼、是否只是重現、失敗可否重用、下一個最小動作。把輪次 FINISHED 與問題 COMPLETED 分開。PR 的本輪來源檢索在合併前過期或新解答出現，重新查後再合。驗證器修改單獨審；原問題不能因不會解而偷偷改成較弱命題。

## 5. 生物、化學與物理的完成界線
廣泛機制問題不能只用 benchmark 分數關閉。第一輪先把 benchmark 子任務、切分、誤差界和適用條件凍結；完成子任務只關子任務。需要新實驗才能辨識的假說標 EXPERIMENT_REQUIRED。計算研究仍可繼續，只是不能把計算結果冒充新實驗或普遍機制。
