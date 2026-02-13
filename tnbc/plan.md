# QSP Model Reproduction Plan
## Metastatic TNBC - Pembrolizumab Monotherapy
### Based on: Arulraj et al., Sci. Adv. 9, eadg0289 (2023)

---

## 1. Project Overview

**Goal:** Arulraj et al. (2023)의 전이성 삼중음성유방암(mTNBC) QSP 모델을 R로 재현한다.
- 원 논문은 MATLAB SimBiology로 구현됨
- R의 `deSolve`, `mrgsolve` 등의 ODE solver 패키지를 활용하여 재현

**Model Scale:**
- 746 species (상태변수)
- 737 parameters
- 1,428 reactions (ODE 정의)
- 314 algebraic rules
- 65 discrete events
- 127 virtual patient parameter distributions

---

## 2. Model Architecture (Compartments)

### 2.1 Compartment Structure (9개 주요 구획)

| Compartment | Description | Initial Volume |
|-------------|-------------|----------------|
| V_C | Central (Blood) | 5 L |
| V_P | Peripheral | 60 L |
| V_T | Primary Tumor | Dynamic (starts ~0.000001 mL) |
| V_T_Ln1 | Lung Metastatic Tumor 1 | Dynamic |
| V_T_Ln2 | Lung Metastatic Tumor 2 | Dynamic |
| V_T_other | Other Metastatic Tumor | Dynamic |
| V_LN | Primary Tumor-draining LN | 1 LN unit |
| V_LN_Ln | Lung Tumor-draining LN | 1 LN unit |
| V_LN_other | Other Tumor-draining LN | 1 LN unit |

### 2.2 Compartment Connectivity
```
Central (Blood) <---> Peripheral
    |
    |---> Primary Tumor <---> Primary LN
    |---> Lung Met Tumor 1 <---> Lung LN
    |---> Lung Met Tumor 2 <---> Lung LN
    |---> Other Met Tumor <---> Other LN
```

---

## 3. Model Species (주요 세포 및 분자)

### 3.1 Cell Types
| Category | Species | Description |
|----------|---------|-------------|
| Cancer | C1-C5 | 5개 cancer clones (각기 다른 성장률, neo-epitope 발현) |
| T cells | nT0 (naive CD4), nT1 (naive CD8) | Naive T cells |
| | T0 (Treg), Th (Helper T) | CD4+ T cell subsets |
| | T1-T8 | Neo-antigen specific cytotoxic T cells (8 specificities) |
| | T1_exh-T8_exh | Exhausted cytotoxic T cells |
| APCs | APC, mAPC | Immature/Mature antigen-presenting cells |
| Macrophages | M1, M2 | M1 (antitumor) / M2 (protumor) macrophages |
| MDSCs | MDSC | Myeloid-derived suppressor cells |

### 3.2 Cytokines & Soluble Factors
| Species | Role |
|---------|------|
| IL-2 | T cell proliferation |
| IL-10 | Immunosuppressive |
| IL-12 | APC maturation, M2->M1 polarization |
| IFN-gamma | Antitumor, PD-L1 upregulation |
| TGF-beta | Immunosuppressive, Treg induction |
| CCL2 | MDSC/Macrophage recruitment |
| Angiogenic factors | Tumor vasculature growth |

### 3.3 Checkpoint Molecules
| Molecule | Expression |
|----------|-----------|
| PD-1 | Cytotoxic T cells, Macrophages |
| PD-L1 | Cancer cells (upregulated by IFN-gamma) |
| PD-L2 | Cancer cells |
| CD28 | T cells |
| CTLA-4 | T cells |
| CD80/CD86 | APCs |
| CD47/SIRPalpha | Cancer cells / Macrophages |

---

## 4. Core Model Equations

### 4.1 Cancer Cell Dynamics
- **Modified Gompertzian Growth:**
  ```
  dC_i/dt = k_C_growth_i * C_i * log(K / C_total)
            - k_C_death * C_i
            - k_C_T * f(T_cyt, C) * C_i
            - k_C_M1 * f(M1, C) * C_i
  ```
  - C_i: i번째 cancer clone
  - K: carrying capacity (종양 혈관에 의해 동적으로 변화)
  - k_C_T: cytotoxic T cell에 의한 cancer cell killing rate
  - k_C_M1: M1 macrophage에 의한 phagocytosis rate

### 4.2 Tumor Carrying Capacity (Angiogenesis)
- **Vasculature-dependent carrying capacity:**
  ```
  dK/dt = k_K_g * c_vas / (c_vas + c_vas_50) - k_K_d * K
  ```
  - c_vas: angiogenic factor concentration
  - Angiogenic factors secreted by cancer cells and M2 macrophages

### 4.3 Naive T Cell Dynamics
- **Central compartment:**
  ```
  dnT/dt = k_nT_source - k_nT_death * nT
           - k_nT_mig * (nT - nT_peripheral)
           - k_nT_LN * nT (trafficking to LN)
  ```

### 4.4 T Cell Activation (in Lymph Nodes)
- **Activation depends on TCR-pMHC ligation:**
  ```
  activation_rate = f_Hill(pMHC_on_mAPC)
  ```
- **Number of divisions (linear sum):**
  ```
  n_div = n_div_TCR + n_div_CD28 + n_div_IL2
  ```
  - CD28 engagement competitively inhibited by CTLA-4
  - IL-2 contribution from helper T cells

### 4.5 Cytotoxic T Cell Dynamics (in Tumor)
- **Infiltration:**
  ```
  k_T_infiltration * V_T * f_vas * T_central
  ```
- **Cancer cell killing:**
  ```
  killing = k_C_T * T_cyt / (T_cyt + C_total * k_T_C_ratio) * C_i
  ```
  - Inhibited by: TGF-beta, Arg-I, NO, PD-1/PD-L1 ligation
- **Exhaustion:**
  ```
  exhaustion_rate = k_T1 * f(PD1-PDL1) + k_T_IL10 * f(IL10)
  ```

### 4.6 Macrophage Dynamics
- **Recruitment:** CCL2-dependent
  ```
  recruitment = k_Mac_recruit * CCL2 / (CCL2 + CCL2_50)
  ```
- **Polarization:**
  - M1 -> M2: promoted by IL-10, TGF-beta
  - M2 -> M1: promoted by IL-12, IFN-gamma
- **Phagocytosis:**
  - Inhibited by IL-10 and CD47-SIRPalpha interaction

### 4.7 MDSC Dynamics
- **Recruitment:** CCL2-dependent
- **Effects:**
  - Release Arg-I and NO
  - Inhibit cytotoxic T cell activity

### 4.8 Cytokine Dynamics
- Each cytokine: secretion - degradation - consumption
  ```
  dCytokine/dt = secretion_terms - k_deg * Cytokine
  ```

### 4.9 PD-1/PD-L1 Checkpoint Module
- **PD-L1 expression on cancer cells:**
  ```
  PD-L1 = PD-L1_base * (1 + r_PDL1_IFNg * IFNg / (IFNg + IFNg_50))
  ```
- **PD-1/PD-L1 binding inhibits:**
  - T cell-mediated cancer cell killing (Hill function)
  - M1 macrophage phagocytosis (Hill function)

### 4.10 Pembrolizumab PK/PD
- **PK (2-compartment model):**
  ```
  dA_central/dt = -CL * A_central/V_C - Q * (A_central/V_C - A_peripheral/V_P)
  dA_peripheral/dt = Q * (A_central/V_C - A_peripheral/V_P)
  ```
  - Dose: 200 mg IV every 3 weeks
- **PD:**
  - Anti-PD-1 binds PD-1 on T cells and macrophages
  - Blocks PD-1/PD-L1 and PD-1/PD-L2 interactions
  - Reduces inhibitory signaling -> enhanced killing/phagocytosis

### 4.11 Diversity Indices (for Biomarker Analysis)
- **Richness (S):** Number of unique cancer clones or T cell specificities
- **Shannon Index:**
  ```
  H = -sum(p_k * ln(p_k))    (k = 1..S)
  ```
- **Evenness:**
  ```
  J = H / ln(S)
  ```

---

## 5. Key Parameters Summary

### 5.1 Physical/Physiological Parameters
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| V_C | 5 | L | Central compartment volume |
| V_P | 60 | L | Peripheral compartment volume |
| vol_cell | 2572 | um^3 | Cancer cell volume |
| vol_Tcell | 175 | um^3 | T cell volume |
| Ve_T | 0.37 | - | Tumor void fraction |
| k_cell_clear | 0.1 | 1/day | Dead cell clearance rate |

### 5.2 Virtual Patient Parameters (127 varied parameters)
| Parameter | Distribution | Median (SD) | Description |
|-----------|-------------|-------------|-------------|
| Initial met diameter | Log-normal | 1.65 cm (0.3) | Initial tumor size |
| Growth rate constants | Log-normal | Clone-specific | Cancer clone growth |
| Neo-antigen concentrations | Log-normal | ~6.75e-14 mol/cell | Antigen expression |
| Seeding time | Sampled | Variable | Metastasis timing |
| PD-L1 expression | Log-normal | Variable | Baseline PD-L1 |

### 5.3 Sensitive Parameters (Response Rate에 영향이 큰 파라미터들)
- k_C_T: Cancer cell killing rate by cytotoxic T cells
- k_Mac_recruit: Macrophage recruitment rate
- k_vas_growth: Tumor vasculature growth rate
- CCL2_50: Half-maximal CCL2 for MDSC recruitment
- div_T0: T cell diversity
- k_Treg: T cell death rate by Tregs
- k_CCL2_deg: CCL2 degradation rate
- agconc_C_0: Self-antigen concentration

---

## 6. Implementation Plan (R Code)

### Phase 1: Infrastructure Setup
**목표:** R 프로젝트 구조 설정 및 필수 패키지 설치

```
TNBC/
├── plan.md                    # 이 파일
├── R/
│   ├── 01_parameters.R        # 파라미터 정의
│   ├── 02_model_ode.R         # ODE 시스템 정의
│   ├── 03_events.R            # Discrete events (dosing, seeding)
│   ├── 04_rules.R             # Algebraic rules
│   ├── 05_initial_conditions.R # 초기조건 설정
│   ├── 06_simulation.R        # 시뮬레이션 실행
│   ├── 07_virtual_patients.R  # Virtual patient 생성
│   ├── 08_biomarker_analysis.R # 바이오마커 분석
│   ├── 09_visualization.R     # 결과 시각화
│   └── utils.R                # 유틸리티 함수
├── data/
│   ├── adg0289_Data_S1.xlsx   # 모델 파라미터 원본
│   └── adg0289_Data_S2.xlsx   # 바이오마커 분석 원본
├── output/
│   ├── figures/
│   └── results/
└── main.R                     # 메인 실행 스크립트
```

**Required R Packages:**
- `deSolve` or `mrgsolve`: ODE solver
- `readxl`: Excel 파일 읽기
- `lhs`: Latin Hypercube Sampling (virtual patients)
- `ggplot2`: 시각화
- `dplyr`, `tidyr`: 데이터 처리
- `parallel`: 병렬 연산 (VP 시뮬레이션)

### Phase 2: Parameter & Species Definition (01_parameters.R, 05_initial_conditions.R)
1. Data_S1.xlsx의 Parameters sheet에서 737개 파라미터 로드
2. Data_S1.xlsx의 Species sheet에서 746개 초기 상태변수 정의
3. Data_S1.xlsx의 Compartments sheet에서 구획 정보 로드
4. 단위 변환 확인 (mole, cell, L, mL 등)

### Phase 3: ODE System Construction (02_model_ode.R)
1. Data_S1.xlsx의 Reactions sheet (1,428 reactions)를 R ODE 함수로 변환
2. **모듈별 구현 순서:**
   - (a) Cancer cell dynamics (Gompertzian growth + killing)
   - (b) Naive T cell trafficking
   - (c) APC dynamics & antigen processing
   - (d) T cell activation in LN
   - (e) T cell infiltration & effector functions in tumor
   - (f) T cell exhaustion
   - (g) Macrophage dynamics (recruitment, polarization, phagocytosis)
   - (h) MDSC dynamics
   - (i) Cytokine dynamics (IL-2, IL-10, IL-12, IFN-g, TGF-b, CCL2)
   - (j) Checkpoint molecules (PD-1, PD-L1, CD28, CTLA-4, CD47)
   - (k) Angiogenesis / carrying capacity
   - (l) Pembrolizumab PK/PD

### Phase 4: Algebraic Rules & Events (03_events.R, 04_rules.R)
1. 314개 algebraic rules 구현
   - Tumor volume 계산 (세포 수 -> volume)
   - Total cell counts (C_total, T_total 등)
   - Derived quantities (ratios, fractions)
2. 65개 discrete events 구현
   - Metastatic seeding events (time-triggered)
   - Pembrolizumab dosing schedule (200 mg Q3W)

### Phase 5: Single Patient Simulation (06_simulation.R)
1. 단일 환자에 대해 ODE 시스템 풀기
2. 치료 전 (tumor growth) 시뮬레이션
3. 치료 (pembrolizumab) 시뮬레이션
4. RECIST v1.1 기준 response 평가
   - CR/PR: tumor diameter 30% 이상 감소
   - PD: 20% 이상 증가
   - SD: 그 외

### Phase 6: Virtual Patient Generation (07_virtual_patients.R)
1. Latin Hypercube Sampling으로 127개 파라미터 샘플링
2. 파라미터 분포: 대부분 log-normal
3. 1000명 VP 생성 -> 816명 필터링 (target tumor diameter 도달)
4. 각 VP별 시뮬레이션 수행 (병렬 처리)

### Phase 7: Clinical Trial Simulation
1. Virtual clinical trial 수행
   - Pembrolizumab 200 mg Q3W
   - Response 평가 (RECIST v1.1)
2. **Calibration targets (KEYNOTE-119):**
   - Overall response rate (ORR)
   - Duration of response
   - Time to response
   - % patients with lung metastases (~65%)

### Phase 8: Biomarker Analysis (08_biomarker_analysis.R)
1. 94개 individual biomarker 후보 테스트
2. Response probability 및 RIS 계산
3. Biomarker combination (pairwise) 분석
4. 결과 ranking 및 Data_S2.xlsx와 비교

### Phase 9: Visualization (09_visualization.R)
1. Fig. 1 재현: Model schematic (diagram)
2. Fig. 2 재현: Calibration (cell fractions comparison)
3. Fig. 3 재현: Virtual clinical trial results
   - Waterfall plot, Spider plot
4. Fig. 4 재현: Treatment effects (responders vs non-responders)
5. Fig. 5 재현: Individual biomarkers
6. Fig. 6 재현: Biomarker combinations
7. Fig. 7 재현: Parameter sensitivity
8. Fig. 8 재현: Primary vs metastatic response comparison

---

## 7. Technical Considerations

### 7.1 ODE Solver Selection
- **mrgsolve 권장:** C++ 기반으로 대규모 ODE에 적합, event handling 우수
- **대안:** deSolve의 `lsoda` (stiff system 대응 가능)
- 746개 상태변수 -> stiff system일 가능성 높음

### 7.2 Computational Considerations
- 1000 VP x (pre-treatment + treatment) 시뮬레이션 = 높은 연산 비용
- `parallel` 또는 `future` 패키지로 병렬 처리 필수
- mrgsolve 사용 시 C++ 컴파일로 속도 향상

### 7.3 Validation Strategy
1. **단위 검증:** 모든 수식의 단위 일관성 확인
2. **단일 환자 검증:** 논문 Fig. S3, S4의 time series와 비교
3. **VP 검증:** 논문 Fig. 3의 response rates, waterfall plot과 비교
4. **바이오마커 검증:** Data_S2.xlsx의 결과와 비교

### 7.4 Key Challenges
1. **모델 규모:** 746 ODE + 314 rules -> 디버깅 난이도 높음
2. **Stiffness:** 세포 역학과 분자 역학의 시간 스케일 차이
3. **Event handling:** Metastatic seeding, drug dosing
4. **Parameter uncertainty:** 86개 파라미터에 강한 실험적 근거 부족
5. **MATLAB -> R 변환:** SimBiology reaction format을 R ODE로 변환

---

## 8. Implementation Timeline (Suggested)

| Week | Task | Deliverable |
|------|------|-------------|
| 1-2 | Phase 1-2: Setup + Parameters | 프로젝트 구조, 파라미터 로딩 |
| 3-5 | Phase 3: ODE System (모듈별) | 핵심 ODE 함수 구현 |
| 6 | Phase 4: Rules + Events | Algebraic rules, events 구현 |
| 7-8 | Phase 5: Single Patient | 단일 환자 시뮬레이션 검증 |
| 9-10 | Phase 6-7: VP + Clinical Trial | 가상 임상시험 시뮬레이션 |
| 11-12 | Phase 8-9: Biomarker + Viz | 바이오마커 분석 및 시각화 |

---

## 9. Reference Resources

- **논문:** Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
- **MATLAB Code:** http://dx.doi.org/10.17632/r46rk4vwdv.1
- **Base Model (Wang et al.):** iScience 25, 104702 (2022)
- **Data S1:** Model reactions and parameters (adg0289_Data_S1.xlsx)
- **Data S2:** Biomarker analysis results (adg0289_Data_S2.xlsx)
- **Clinical Trial:** KEYNOTE-119 (Winer et al., Lancet Oncol. 2021)

---

## 10. Data S1 Excel Structure (참조)

| Sheet | Rows | Description |
|-------|------|-------------|
| Compartments | 131 | 구획 정의 (이름, 용량, 단위) |
| Species | 747 | 상태변수 정의 (초기값, 단위, 위치) |
| Parameters | 738 | 파라미터 (값, 단위, 설명) |
| Uncertain parameters | 87 | 불확실 파라미터 |
| Reactions | 1429 | ODE 반응식 |
| Rules | 315 | 대수 규칙 |
| Events | 66 | 이벤트 (seeding, dosing) |
| Parameter distributions | 128 | VP 생성용 분포 |
