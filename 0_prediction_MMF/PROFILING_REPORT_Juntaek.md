# MATLAB 시뮬레이션 리소스 프로파일링

> 대상: `main_minrate.m` (계측본: `main_minrate_profiled.m`)
> 프로젝트 루트: `MMF_in_the_Dark_Ver_final`

## 1. 진행 프로젝트 개요
* **프로젝트명/주제:** MMF (Max-Min Fairness) beamforming under imperfect CSI — RSMA(Rate-Splitting Multiple Access) 기반 하향링크 빔포밍 ("MMF in the Dark")
* **주요 목적:**
  * 채널을 완전히 알지 못하는(imperfect/reconstructed CSI) 상황에서 **최소 사용자 전송률(min spectral efficiency)을 최대화**하는 빔포머 설계 성능 비교.
  * NOMP(`trivariate_NOMP`) 기반 채널 추정 + Fisher/CRLB 기반 **오차 공분산 행렬(ECM, `H_hat_err_cov`)** 을 반영한 **제안 기법 MMF-GPI-RS** 를, 기준 기법들(WMMSE, RZF, MRT, RS 미적용 GPI, ECM 미반영 버전)과 SNR(25~40 dB)에 따라 비교.
  * 결과를 `figure(1)` (SNR vs. minimum spectral efficiency) 및 `result.mat` 로 저장.

## 2. 주력 데이터 스케일 (Matrix Size)
설정값: `N = 12` (안테나), `K = 4` (사용자), `Ns = 128` (subcarrier), `L = 5` (path), `iter = 15`, SNR 4점(25/30/35/40 dB) → **외부 루프 총 60회 반복**.

| 배열 / 연산 위치 | 차원 | 타입 | 대략 크기 |
|---|---|---|---|
| `S = eye(M)`, `M = Ns*N` (`trivariate_NOMP.m`) | **1536 × 1536** | real double | **≈ 18.9 MB** (최대 단일 행렬, 채널당 K회 임시 생성) |
| 재구성 채널 벡터 `h_reconst_dl` / `y_full` (kron 결과) | 1536 × 1 | complex double | ≈ 24 KB (다수) |
| 빔포밍 KKT 행렬 `A,B,C,D{k}` (`MMF_GPI_RS`, `func_rate_compute`) | **60 × 60** = N(K+1) | complex double | ≈ 57.6 KB × 16개 ≈ 0.9 MB |
| `H_ch`, `H_hat` | 12 × 4 | complex double | < 1 KB |
| `H_hat_err_cov` (ECM) | 12 × 12 × 4 | complex double | ≈ 9 KB |
| Fisher 정보행렬 (`get_fisher`) | 20 × 20 (4L × 4L) | real double | < 4 KB |
| Jacobian `G` (`get_G`) | 20 × 12 | complex double | < 4 KB |

* **데이터 타입:** 대부분 **complex double** (채널·빔포머·KKT 행렬). `eye(...)`·rate 값·Fisher 등 일부만 real double. **single 정밀도는 미사용.**
* **예상 메모리 점유량 (알고리즘 워킹셋):** 순수 알고리즘이 동시에 잡는 배열은 **약 20~30 MB 수준** (지배 요인 = `eye(1536)` ≈ 19 MB, 순차적으로 생성·해제됨). 실제 측정되는 Peak RSS 는 여기에 **MATLAB 런타임 + CVX 솔버 기반 메모리(수백 MB~수 GB)** 가 더해진 값 → 5절 측정치 참조.
  * 참고(최적화 힌트): `trivariate_NOMP.m`의 `S = eye(M); y = S*y_noisy;` 는 사실상 `y = y_noisy` 와 동일하므로, 19 MB짜리 항등행렬 생성·곱을 제거하면 채널 생성 단계 메모리/시간을 줄일 수 있음.

## 3. 핵심 연산 및 기술
* **가장 부하가 큰 내장 연산:**
  1. **`pinv` (SVD 기반 유사역행렬)** — 반복 while 루프 안에서 다수 호출: `MMF_GPI_RS`(블록별 12×12, 매 반복 K+1회), `without_RS`(N·K×N·K = 48×48 전체 pinv 매 반복), `zero_forcing`, `WMMSE`, `channel_gen_est`(Fisher 20×20 pinv). → **가장 반복 횟수가 많은 무거운 연산.**
  2. **CVX 볼록 최적화 (`cvx_begin ... cvx_end`)** — `WMMSE_ver6.m`. SDP/원뿔 프로그램을 반복적으로 풀며 **단일 호출당 시간·메모리 비용이 가장 큼**.
  3. **`kron` (Kronecker product)** — 1536-length 채널/신호 벡터 합성, `channel_gen_est`·`trivariate_NOMP` 에서 L회 × 반복.
  4. **`get_fisher` 의 4중 중첩 루프** — `L²·N·Ns ≈ 25×12×128 ≈ 3.8만` 회 `get_der` 평가 (내장 연산은 아니지만 순수 MATLAB 루프 CPU 부하 큼).
  5. 조밀 행렬 곱 `f'*A*f` (60×60) 를 GPI while 루프에서 반복.
  * ※ `svd`, `eig`, `fft`, `inv` 직접 호출은 없음 (`pinv` 내부에서 SVD 사용).
* **주요 사용 Toolbox / 외부 패키지:**
  * **CVX** (외부 볼록최적화 패키지, `./cvx/`) — `WMMSE_ver6.m` 에서 사용. *(MathWorks 제품 아님)*
  * **Statistics and Machine Learning Toolbox** — `unifrnd` (`channel_gen_est.m`, `trivariate_NOMP.m`).
  * 그 외 `pinv`, `kron`, `blkdiag`, `nchoosek`, `randn`, `reshape` 등은 **Base MATLAB**.

## 4. 병렬 처리 (Parallel Computing)
* **parfor 사용 여부:** **X (미사용)** — 모든 루프가 `for` 직렬 실행.
* **Worker(코어) 수:** 해당 없음 (Parallel Computing Toolbox·`parpool` 미사용). 단, `pinv`·행렬곱 등은 **MATLAB 내장 BLAS/LAPACK 의 암묵적 멀티스레딩** 으로 자동 병렬화됨.
* **Worker당 다루는 데이터 크기:** 해당 없음 (단일 프로세스·단일 워커).

## 5. 실행 결과 (실제 측정치)
> `main_minrate_profiled.m` 실행 완료 시 콘솔 하단에 자동 출력됩니다. 아래 값을 채워 넣으세요.

* **실행 시간:** `2022.72 초 ( 33.71 분 )`  ← 콘솔 `총 실행 시간 (Elapsed Time)`
* **Peak RAM 사용량:** `1310.3 MB ( 1.280 GB )`  ← 콘솔 `Peak RAM (Peak RSS)`

### 측정 방식 요약
* `tic/toc` 로 전체 실행 시간(wall-clock) 측정.
* 0.5초 주기 백그라운드 `timer` 가 macOS 에서 `ps -o rss= -p <PID>` (Windows: `tasklist`) 로 **프로세스 실제 점유 메모리(RSS)** 를 샘플링, 최댓값을 Peak RAM 으로 기록.
* Peak RSS 에는 MATLAB 엔진 + CVX 솔버 상주 메모리가 모두 포함됨 (프로세스 전체 관점). 순수 알고리즘 워킹셋은 2절 참조.
