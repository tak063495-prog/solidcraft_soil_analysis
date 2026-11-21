# REPRODUCEの補助監査コマンド

`REPRODUCE.ja.md`とSPECは新規求解開始時の束を保持する。
ここではその後追加した集計/生産CLI照合のコマンドを示す。新しい物理・許容値ではない。
同REPRODUCEの環境設定後、リポジトリルートで実行する。

```powershell
.venv\Scripts\python "$B/cli_controls.py" runs/NEW-b01-fresh runs/NEW-b01-cli-controls --cli target/release/soil-fem-3d.exe
.venv\Scripts\python "$B/cli_binding_audit.py" runs/NEW-b01-fresh runs/NEW-b01-cli-controls runs/NEW-b01-cli-binding.json
.venv\Scripts\python "$B/work_audit_b01.py" runs/NEW-b01-fresh runs/NEW-b01-work-audit.json
.venv\Scripts\python "$B/verify_contract.py" runs/NEW-b01-fresh runs/NEW-b01-probes/scope.json runs/NEW-b01-mutations/mutations.json runs/NEW-b01-contract.json
.venv\Scripts\python interim/2026-12-02-tri6-offset-kink-b03/report_validation.py runs/NEW-b01-rust runs/NEW-b01-python runs/NEW-b01-validation.json
```

追加5〜20分、RAM 4〜8 GB（未実測）。通常CLIで5ケースを新たに求解する。
CLI自身のcheckpoint-10を別PIDで再開し、時間以外のRun全体を厳密比較する。
観測版と通常CLIは別実行ファイルなので、最初のRun完全一致はfalseのまま保存する。
補助監査では各実行版SHAを実ファイルへ結び、差が`provenance.solver_binary_sha256`だけであることを検証する。
他の全場/監査/履歴の差を許可しない。同じ実行版の再開比較からSHAを除外しない。

`work-audit`は各受理端点の台形観測であり、設置原点を各履歴自身から求める。
`contract`は全受理ID/土GP/板評価点の網羅と実多面体幾何を検証する。
`validation`は全発見ID・実行IDを照合する。準備エラー/未実行/skipを成功へ変えない。
この監査の成功をB-01全体の数値受入へ転用しない。
