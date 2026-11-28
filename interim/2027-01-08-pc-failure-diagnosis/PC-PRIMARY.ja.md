# PC一次資料：依頼本文の埋め込みログ

出典は本依頼で提示されたログ要旨。PCのファイルを取得・再読込した記録ではない。
PCのPID 18836へ接続・操作していない。以下の旧SHA・保存数値を変更しない。
日付は2026-11-18、時刻はJST (Asia/Tokyo)。2027-01-08は指定保存先ラベル。
11-15の`/workspace/pc_runs/pc_queue_small_20261115-results.zip`は旧フルcheckoutが2GB capで
status=memoryとなった別記録で、今回の6件を含まない。この旧メモリ失敗を今回へ転用しない。

## b02-n4

07:16:22–07:16:48、cwd `D:\Work\NEW-B02`、SHA `e1e49b78f3f5138f662e94dd4125cc07afae1b9c`。
timeout 7500s、memcap 8,589,934,592 bytes、elapsed 26.6s、peak_ws 131,485,696、rc 1、
timed_out=false、memory_exceeded=false。

```text
.venv\Scripts\python.exe interim/2026-12-24-b02-offset-composite-cause-diagnosis\b02cause_driver.py --output D:\Work\NEW-B02\runs\NEW-b02cause-PC-n4 --observer D:\Work\NEW-B02\target\release\b02-observer.exe --cli D:\Work\NEW-B02\target\release\soil-fem-3d.exe --steps 4
b02cause_driver.py:97: assert process['pid'] != native['pid']
AssertionError
```

command()はPopen後にwaitし終了済みp.pidを記録。全exit_code=0：native=1156 (3.710s)、
independent=22044 (0.707s)、plain=18148 (0.281s)、resumed-6=17284 (3.161s)、
resumed-10=1156 (2.704s)。終了済みnativeとresumed-10が同じPID。
tet12-n4の上記出力は存在する。後続のexact比較（98–99行）は未到達。
tet96-n4、steps8/16は未到達。

## b01-evidence

07:19:54–07:19:59、cwd `D:\Work\NEW-B01-curved`、SHA `3747630ab615b4be2117ca1aaa13e798cd436374`。
timeout1800s、memcap4,294,967,296、elapsed2.1s、peak_ws0、rc1。

```text
.venv\Scripts\python.exe interim/2026-12-25-b01-curved-tri6-cause-diagnosis\b01cause_evidence.py --verify interim/2026-12-25-b01-curved-tri6-cause-diagnosis --extract D:\Work\NEW-B01-curved\work\b01cause-pc-unsealed
b01cause_evidence.py:23: if p.stat().st_size != part['size'] or sha(p) != part['sha256']:
FileNotFoundError: [WinError 2] ... interim\2026-12-25-b01-curved-tri6-cause-diagnosis\evidence.zip.part000
```

説明は「missing evidence is recorded, later steps still run」だが例外停止。
unittestは11件中10件ok、同じ欠落によるtest_sealed_evidence_and_false_status_are_separateのみ失敗。
remeasureは未実行。当該partはsparse除外。

## b01-cli

同じSHA/cwd。07:20:05–07:20:07、timeout1800s、memcap4,294,967,296、elapsed2.1s、peak_ws0、rc1。

```text
.venv\Scripts\python.exe interim/2026-12-25-b01-curved-tri6-cause-diagnosis\b01cause_cli_controls.py --solver D:\Work\NEW-B01-curved\target\release\soil-fem-3d.exe --output D:\Work\NEW-B01-curved\work\b01cause-pc-cli-controls
b01cause_cli_controls.py:40: raw = (PRIOR/part['path']).read_bytes()
FileNotFoundError: [Errno 2] ... D:\Work\NEW-B01-curved\interim\2026-12-03-tri6-curved-b01\evidence.zip.part000
```

sparse除外のpartが不在、求解前に停止。

## b01-roof-n8

同じSHA/cwd。07:20:16–07:20:22。cases生成rc0 (2.1s)。
求解：`target\release\soil-fem-3d.exe solve-stages ...\b01cause-pc-faceted-roof-n8\{mesh,ground,stages,run}.json`。
timeout1800s、memcap8,589,934,592、elapsed4.2s、peak_ws113,528,832、rc1。

```text
numerical solve failed: load_plate_beam: cutback exhausted: numerical solve failed: beam interior true residual exceeded tolerance
```

run.json：n_frames=2、n_attempts=12。last_committed：stage_index=1、stage_progress=1.0、boundary_factor=0.0。
last_attempt：stage=load_plate_beam、progress_from=0.0、progress_to=0.001953125、accepted=False、
reasonは上記、committed_state_unchanged_on_rejection=True。
n=2/4は各5frames、16accepted/0rejected。真の曲面殻でなく平面多面体屋根。

## a04-matrix

07:41:09–07:41:27、cwd `D:\Work\NEW-A04`、SHA `33b185b060c26cf2d0957271391b7bec86e491da`。
timeout7500s、memcap8,589,934,592、elapsed18.0s、peak_ws170,000,384、rc1。

```text
.venv\Scripts\python.exe interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis\a04cause_driver.py --output D:\Work\NEW-A04\runs\NEW-a04-cause-PC
study.py --steps 4 8 16: exit_code=1 (7.24s)
native: construction complete: 20 frames, 1 removals
2026-11-29-richards-wet-joint-a04/soil_reference.py:175:
    result=m.step(v['h'],v['trace'],v['storage'],period,step)
2026-11-22-wet-joint-transient-water/hydraulic_reference.py:160:
    raise ArithmeticError(('hydraulic line search', merit))
ArithmeticError: ('hydraulic line search', 0.011860214796293166)
```

soil/VG-n4。外側RuntimeErrorは子の非ゼロ終了の包装。soil.logのzipfile Duplicate name警告
（study.py等）は終了原因ではない。32/64/128未実行。

## a01-evidence

08:45:37–08:45:42、cwd `D:\Work\NEW-A01`、SHA `7b254b1cdbe9c288989631944cd698aa861f83ac`。
verify timeout1800s、extract-old timeout3600s、各memcap8,589,934,592、elapsed2.1s、peak_ws0、rc1。

```text
a01cause_evidence.py --verify interim/2026-12-20-a01-joint-combo-cause-diagnosis
FileNotFoundError: interim\2026-12-20-a01-joint-combo-cause-diagnosis\evidence.zip.part000
a01cause_evidence.py --extract-old interim/2026-11-23-joint-combo-extension --output D:\Work\NEW-A01\work\NEW-A01-old-verified
FileNotFoundError: interim\2026-11-23-joint-combo-extension\evidence.zip.part000
```

a01-g32/7route/stacked/tet10はdone。Tet4求解は止まっていない。

## 共通の取得条件・未実行

sparseはinterim下のzip/zip.part*/tar/tar.gzとそのpart、tgz/gz/7z/npz、全*.part*を除外。
フルcheckoutは先端付近約29GB、証拠アーカイブ約27.66GB、source-ish約312–461MBという依頼上の規模。
今回の6件をtimeout/memcap起因とはしない。peak_ws=0は計測値の記録であり、必要メモリ0の主張ではない。
b02-n8/n16、a04-steps-32/64/128は先行失敗によるskipped・未実行。今回の自動再実行対象にしない。
