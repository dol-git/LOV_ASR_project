# gREST 실행 및 수렴 점검 기록

검토일: 2026-10-21. 기존 입력·production log·분석 코드·저장된 보고서와
진단 CSV를 읽어 정리했다. 분석 코드를 실행하거나 DCD를 다시 분석하지 않았다.
아래 실행 길이는 입력의 목표값뿐 아니라 각 replica 로그의 마지막 완료 step을
확인한 값이다. 구조 지표의 수렴 평가는 기존 보고서의 증거 검토이며 독립 재현이 아니다.

**해당 계산 길이와 평가 지표에서 수렴을 확인하지 못해,
정량적인 ensemble 비교의 주근거로 사용하지 않았다.**

## 실행 범위와 조건

| 항목 | WT | Anc |
|---|---|---|
| 검토한 50 ns production 세트 | 1개 | 1개 |
| 교환 execution replica 수/세트 | 4 | 4 |
| 각 replica 실제 마지막 step | 12,500,000 | 12,500,000 |
| timestep / 명목 완료 길이 | 0.004 ps / 50 ns | 0.004 ps / 50 ns |
| 저장 프레임/replica | 5,000, 0.01–50.00 ns | 동일 |
| bath / thermostat | NVT, 298.15 K, BUSSI | 동일 |
| solute parameter state 1–4 | 298.15 / 331 / 365 / 400 K | 동일 |
| REST 설정 | 1 dimension, `param_type1=D L` | 동일 |
| 실제 production tempering selection | PROA residue 120–142 | PROA residue 120–142 |
| 대응 Jα 분석 구간 | 120–142 | 119–141 |

입력은 CHARMM, PME, vdW force switch, switch/cutoff/pairlist 10/12/13.5 Å,
VVER, rigid bonds, PBC를 지정한다. HMR는 `hydrogen_mr=YES`, `hmr_ratio=3.0`,
`hmr_ratio_xh1=2.0`, `hmr_target=solute`다. 이는 읽은 입력 설정이며 실제 사용된
모든 force-field 파일·생산 바이너리의 출처나 물리적 정확성 검증을 대신하지 않는다.
좌표 저장과 교환 주기는 각각 2,500 step=10 ps다. 모든 8개 50 ns log의
마지막 INFO step은 12,500,000, 출력 시간은 50,000.0024 ps다. 명목 시간축은
step×0.004 ps로 정의하며 소수 오차를 추가 sampling으로 세지 않는다.

이전 10 ns pilot도 WT/Anc 각각 4 replica가 존재한다. 8개 pilot log에서
2,500,000 step(출력 10,000.0005 ps) 완료를 확인했고, 기존 manifest/검증 기록은
각 1,000 frame, 10 ps 간격을 기록한다. 10 ns pilot과 50 ns 첫 10 ns의
거리 trace는 동일하지 않다. 그러나 restart/초기조건의 독립성을 확립하지 않았으므로
pilot을 추가 독립 반복으로 세거나 10+50=60 ns 연속 run으로 연결하지 않는다.
pilot의 모든 생산 조건이 50 ns와 같다고 가정하지 않는다. 짧은 시험·평형화·archive
복사본 역시 추가 독립 run으로 집계하지 않는다.

## replica, walker, physical state의 구분

- **Execution replica / walker:** 해당 실행 replica의 원본 DCD와 REM history를
  함께 추적한다. 이 구현의 기록은 좌표 출력 뒤 parameter assignment를 갱신한다.
  같은 walker의 Hamiltonian state는 교환에 따라 달라진다. 4개 walker는 교환으로
  결합된 한 세트이며 독립 production run 4개가 아니다.
- **Parameter state:** 고정된 REST 조건이다. 상위 state의 K 표시는 선택 항의
  유효 solute tempering parameter이며 bath 전체를 그 온도로 가열한 별도 NVT가
  아니다. `D L`은 선택된 dihedral/LJ 항을 tempering한다.
- **Physical state:** `remd_convert`, `convert_type=PARAMETER`, `convert_ids=1`로
  매 저장 시점 state 1을 가진 replica를 선택한다. state 1은 298.15 K baseline다.
  GENESIS 소스의 `coeff_full=T_bath/T_rest`, `coeff_half=sqrt(coeff_full)`와
  입력에 따르면 state 1 scaling은 unity다. 이는 소스/입력 유도이며 별도의 runtime
  scale 출력값 검증은 아니다.

저장 frame i(1-based)의 step=i×2500에는 REM step=(i−1)×2500의 상태를 적용한다.
기존 QC는 REM/log 전이 일치 및 변환 좌표·box와 원본 state-1 소유 replica의 일치를
전 프레임에서 확인했다. 후속 walker audit도 state-1 owner mapping을 대조했다.
변환 출력 DCD의 기본 시간 header는 실제 시간을 나타내지 않으므로
`physical_frame_mapping.csv`의 10 ps/frame을 사용해야 한다.

Physical-state 궤적은 서로 다른 walker의 좌표를 이은 parameter-sorted sample다.
연속 단일 분자의 kinetics·전이율·수명으로 해석하지 않는다. 같은 walker의
high-state와 이후 state-1 구조 유사성도 기하학적 provenance 점검이며 물리적
복귀 시간이나 구조적 평형의 증거가 아니다.

## 확인된 sampling 진단

작은 공개 진단표: [exchange_diagnostics.csv](../analysis/grest/exchange_diagnostics.csv).
기존 `qc/{WT,ANC}/exchange_acceptance.csv`의 6행을 그대로 옮겼다.
수용률을 새로 계산하거나 교환 이벤트를 다시 분석하지 않았다.

| 시스템 | state 1↔2 | state 2↔3 | state 3↔4 | 각 쌍 시도 수 |
|---|---:|---:|---:|---:|
| WT | 24.68% | 25.84% | 28.16% | 2,500 |
| Anc | 21.68% | 25.36% | 28.56% | 2,500 |

원래 코드는 `(step, state pair)`로 양쪽 replica의 교환 기록을 중복 제거하고
결과 일치를 확인한다. 기존 summary는 REM/log 불일치 0건, 양 끝 state 간
비중복 A→B→A 왕복 WT 212회·Anc 211회를 기록한다. 각 walker가 온도 사다리를
이동한 occupancy/trace도 저장돼 있다. 마지막 불완전 왕복은 별도 censored 기록이며
마지막 교환 뒤 상태 점에 추가 체류 시간을 주지 않는다. 이 지표들은 **parameter-space
mixing**이 일어났음을 보여 주지만 느린 구조 전이의 충분한 sampling이나 평형을
보장하지 않는다. 이 문서는 구조 효과크기·서열·좌표 표를 공개하지 않는다.

## 수렴 평가의 지표와 비교 대상

수렴의 주된 점검 대상은 **같은 시스템의 같은 physical state 1**에서 시간 구간별
분포가 얼마나 재현되는지다. 서로 다른 유효 온도 state의 분포 차이 자체는
수렴 실패 근거가 아니다. 목표 Hamiltonian이 다르므로 차이가 날 수 있으며,
재가중 없이 전체 state를 섞어 baseline ensemble로 쓰지 않는다.

| 기존 점검 | 비교·정의 | 증거가 허용하는 해석 |
|---|---|---|
| Jα RMSD | core backbone 비가중 fitting 후 Jα N/CA/C/O, 자기 reference; Jα 재정렬 없음 | 시간 block 변화와 fitting 민감도 점검 |
| RMSF | 전체·block별 fluctuation/profile 비교 | 구간 의존성; 평균 하나만으로 수렴 판정 불가 |
| 접촉 거리 | 7개 대응 residue pair 중원자 최소거리; 다섯 10 ns block의 median/IQR/ECDF | 일부 pair의 block 의존성, 모든 pair가 실패했다는 뜻은 아님 |
| pilot 대 50 ns | 같은 pair/원자 선택의 거리 분포 | 별도 관측 집합 비교; 독립 표본 검정 아님 |
| 고정 cluster 분율 | 전체 state-1에서 정의한 label을 모든 block에 유지; 1.00/1.25 Å cutoff | block별 population 안정성의 기술적 점검 |
| 확장 descriptor audit | 접촉 수, COM, Jα 방향, DSSP H fraction, PCA 축 품질; block/cumulative 분포 | descriptor별 drift/block sensitivity; 안정으로 분류돼도 평형 입증 아님 |

50 ns 거리 보고서는 (0,10]…(40,50] ns의 5개 block을 사용한다. ECDF 차이는
sup|F−G|, JSD는 log2 divergence(bits; 제곱근 distance 아님), TV는 ½Σ|p−q|다.
WT/Anc cluster label은 독립 정의라 두 시스템 사이 cluster 벡터 JSD는 계산하지
않는다. 전체 분포는 비교 block을 포함하고 cumulative 표본도 겹친다. 따라서
전체와의 차이가 작아지는 것만으로 수렴을 주장하지 않고 block 간 증거를 우선한다.
연속 프레임 독립성을 가정한 p값·SEM이나 사전 합격 threshold는 제공되지 않았다.

기존 보고서에는 동일 state-1 내 일부 거리 분포의 위치/폭 변화와 고정 label
분율의 block 의존성이 남아 있다. 일부 접촉·구조 지표는 비교적 재현되므로 모든
관측량의 수렴 실패로 일반화하지 않는다. 이러한 증거와 독립 반복 부재는 안정된
정량 ensemble 비교를 주근거로 삼지 않는 **보수적 판단**을 지지한다. 단순히
분포가 다르다는 이유로 엄밀한 통계적 비수렴을 입증한 것은 아니다.

## 기존 판단의 타당성과 미해결 사항

1. 이전 README의 “convergence problems”는 범위가 모호했다. 이 갱신에서는
   “해당 계산 길이와 평가 지표에서 수렴을 확인하지 못해, 정량적인 ensemble
   비교의 주근거로 사용하지 않았다.”로 한정한다. HMR의 수렴이 입증됐다는 뜻도 아니다.
2. WT 최초 fit core는 1–119, Anc는 1–118이다. N말단 영향 점검 후 stable core는
   WT11–119/Anc11–118, 대응 core sensitivity는 WT12–119/Anc11–118을 썼다.
   이들은 최근 HMR의 core100 Cα 정의와 다르다. 최초 core RMSD의 high tail만으로
   전역 core 붕괴·변성을 주장할 수 없고 정의가 다른 통계를 섞으면 안 된다.
3. Anc 생산 tempering 120–142와 대응 Jα119–141의 불일치는 실제 입력 및 stdout에서
   확인된다. 후처리 구간 변경은 생산 Hamiltonian을 고치지 않는다. 수렴 여부와 별개의
   비교 조건 문제이며 영향의 크기는 여기서 확정하지 않는다.
4. 초기 gREST pilot 문서의 기본 Jα120–142 표기와 후속 대응 분석119–141을 구분한다.
   최초 50 ns RMSF 보고서에는 DSSP 미산출이 남아 있지만 이후 별도 reanalysis에는
   MDTraj DSSP가 있다. 서로 다른 단계의 보고서를 같은 분석 버전으로 취급하지 않는다.
5. 유효 독립 표본수·불확실도·독립 run 재현성·rare-state 안정성은 확인되지 않았다.
   교환률·왕복·같은 state에서의 walker 차이나 소유 비율 변화만으로 성공/실패를
   확정하지 않는다. 고정 adduct는 암회복 반응이나 transition state를 직접 sampling하지
   않으며, clustering Shannon H는 구조 열역학 entropy나 Eyring ΔS‡의 직접 증거가 아니다.
6. 후속 검증에는 대응 tempering 영역·입력/힘장 출처 확인, 독립 초기조건의 교환 세트,
   사전 정의한 동일-state 지표와 상관을 고려한 불확실도 평가가 필요하다. 여기서는
   새 계산·연장·재가중·FES 생성·분석 재실행을 하지 않았다.

## 코드, 외부 자료와 실행 순서

[코드 provenance](../analysis/grest/code_provenance.csv)는 원본 SHA256을 기록한다.
세 Python 파일은 byte-identical 스냅샷이다. 실행하지 않고 구문만 검사했다.
[launch 예제](../analysis/grest/submit_50ns.example.sh)는 기존 WT/Anc scheduler
script에서 개인 경로·allocation account·site module을 제거하고 입력과 binary를
환경변수로 받도록 만든 예제다. 과거 실행 파일의 동일 사본이 아니며 job 제출로 검증하지
않았다. 생산 입력·좌표를 제공하지 않으므로 clone만으로 계산을 재현할 수 없다.

스냅샷은 `analysis/grest/Heo/analysis/grest50/`에 있다. 이 sibling layout을
별도 새 외부 작업공간에 복사하고 필요한 원본 자료를 배치하는 것이 전제다.
원본 코드의 실행 동작은 다음과 같으며 이번 작업에서 실행한 명령이 아니다.

- `analyze.py prepare`: snapshot inventory, DCD record/frame 점검, log/REM 교환
  QC와 trace·왕복 집계, 변환 input 생성. 출력 파일을 덮어쓸 수 있다.
- `analyze.py convert WT` / `convert ANC`: 외부 GENESIS `remd_convert` 실행 및
  변환 좌표 검증. `analyze.py immutability`: 기존 inventory와 원자료 hash 점검.
- `contact_distances/convergence_10ns_50ns/analyze_convergence.py`: 기존 50 ns
  거리 CSV와 pilot DCD/REM/PSF를 사용한다. pilot 거리도 재계산하고 그림을 생성하므로
  단순 CSV 요약 명령이 아니다. 출력은 대체로 exclusive create를 사용한다.
- `clustering/block_population_stability/analyze.py`: 기존 membership·population·
  physical frame mapping으로 분율/JSD/TV를 계산한다. 새 clustering은 하지 않으나
  보고서와 그림을 생성한다. module top-level에서 작업하므로 import도 실행에 해당한다.

필요한 외부 자료: `Fugaku_Heo_20260928/HMR_{WT,ANC}/`의 생산 DCD/REM/log/restart,
production input·stdout, 같은 snapshot의 `input/` PSF/PDB 및 force-field,
`tools/genesis-2.1.6/bin/remd_convert`, physical-state 변환 결과·mapping,
기존 거리 CSV·원자 선택, cluster membership·population, pilot 분석 폴더와 원본 파일.
완전한 dependency bundle이나 재현용 환경 export는 포함하지 않았다. 코드는 NumPy,
SciPy, Matplotlib 및 Python 표준 라이브러리를 사용하며 원 production binary/build의
정확한 검증과 전체 패키지 환경 고정은 외부 기록으로 남는다.

[외부 증거 index](../analysis/grest/external_evidence_index.csv)는 읽은 보고서·metadata의
작업공간 상대경로와 SHA256만 제공한다. `Heo/analysis/grest50/REPORT.ko.md`,
`FINAL_EVIDENCE_AUDIT.ko.md`, `walker_provenance/REPORT.ko.md`,
`contact_distances/convergence_10ns_50ns/REPORT.ko.md`,
`clustering/block_population_stability/REPORT.ko.md` 및 별도
`reanalysis_20261006/04_phase3A_ensemble/` 문서가 핵심 근거다.
production 입력·8개 log의 완료 step도 직접 읽어 확인했지만 파일 자체는 제외했다.
공개 허용이 불명확한 서열·좌표·생산 입력·구조 수치 결과, 원문 보고서, 전체 로그,
DCD/restart/checkpoint, 환경 폴더, 개인 경로·account 정보는 추가하지 않았다.
공개 수치는 이 요청 범위의 실행 조건 및 작은 교환 sampling 진단으로 한정했다.
