# FrontierBiology
前沿生物

## 主線：公開資料 → 基線重現 → 跨條件驗證 → 可檢驗假說
第一輪選 **BIO-001：單細胞擾動預測的跨情境泛化**，先核對資料授權與切分、重現簡單基線，再考慮新模型。不是同時重訓十個大模型，也不把預測相關性寫成已發現生物機制。

本批10題是依原始研究整理的、尚未完整解決的研究缺口；**問題的精確研究範圍由本庫題卡定義，不冒稱是來源作者提出的同一條正式猜想**。`OPEN` 是初始有界篩查狀態，每輪必須重查。已有特定方法／資料集的成功不代表所有情境都已解。

| ID | 問題 | 優先級 | 起步任務 |
|---|---|---|---|
| [BIO-001](problems/BIO-001/problem.json) | 單細胞擾動反應跨細胞情境泛化 | A | 留出細胞類型／donor 的簡單基線 |
| [BIO-002](problems/BIO-002/problem.json) | 基因調控網路的因果可識別性 | A | 合成真值與公開擾動資料分開驗證 |
| [BIO-003](problems/BIO-003/problem.json) | 從單細胞快照推斷細胞命運 | B | 動力學假設與獨立 lineage 對照 |
| [BIO-004](problems/BIO-004/problem.json) | 非編碼變異跨組織調控效應 | B | 染色體／細胞情境留出及校準 |
| [BIO-005](problems/BIO-005/problem.json) | 良性蛋白多重突變的高階上位性 | B | 加成基線對未知多突變組合 |
| [BIO-006](problems/BIO-006/problem.json) | 蛋白構象系綜與環境條件 | B | 系綜可觀測量與留出條件比較 |
| [BIO-007](problems/BIO-007/problem.json) | 未知 RNA 摺疊與複合體結構 | B | 盲測模板洩漏與幾何合法性 |
| [BIO-008](problems/BIO-008/problem.json) | 老化標記的跨族群效度與因果解釋 | C | 公開去識別資料的跨研究校準 |
| [BIO-009](problems/BIO-009/problem.json) | 生態／微生物動力系統的可識別性 | B | 非病原合成系統的參數與預測退化 |
| [BIO-010](problems/BIO-010/problem.json) | 物種分布的時空外推與觀測偏差 | A | 空間區塊留出、偏差與校準基線 |

A/B/C 是起步順序而非成功機率。每題記 statement、已知結果、缺口、最小任務、evaluator、限制、完成條件與来源。

## 開始一輪
讀 AGENTS、STATUS、VALIDATION 和治理 pin。使用獨立分支／工作目錄，先看已有 PR 與失敗紀錄。

```bash
# 治理庫放旁邊，checkout GOVERNANCE.lock.json 的精確 commit
python3 ../FrontierLab-Governance/tools/frontier.py validate .
python3 ../FrontierLab-Governance/tools/frontier.py start . BIO-001 --agent gpt
# 真正查本輪最新文獻：general / discipline / solution / criticism，填 round.json
python3 ../FrontierLab-Governance/tools/frontier.py admit . runs/<round-id>/round.json
```

發現同範圍問題已被外部可靠解決，立即記 `COMPLETED_EXTERNAL`，附原作者、來源與核查；只見預印本解答聲稱則 `CLAIMED_RESOLVED`，先停下核實。某benchmark達標只關該子任務，不把父問題整體關閉。

## 證據與安全
只做合法公開、低風險、計算型研究。蛋白與核酸任務限良性基準的預測／驗證；不設計或增強病原、毒素或危害性功能。不自動跑濕實驗、人體／動物介入，不作臨床診斷或治療建議。

資料有自己的使用條件。論文可讀不代表個體資料可下載；無授權則記 BLOCKED。donor、研究批次、同源群與時間依任務切分，避免同一受試者／近親序列同時進訓練及測試。一次高分不是新機制。

## 本次交付
題卡與研究流程已建立；原創研究輪次0，各題分析／evaluator尚未實作，資料未批量下載。共用CI驗記錄，不會替 agent 搜尋或替生物實驗背書。沒有常駐 agent、付費服務或研究自動排程。
