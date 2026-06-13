
# 🚀 RENE DAQ: RAW to Flat PRD ROOT Converter

![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)
![Uproot](https://img.shields.io/badge/Uproot-5.x-yellow.svg)
![Awkward](https://img.shields.io/badge/Awkward_Array-2.x-orange.svg)
![NumPy](https://img.shields.io/badge/NumPy-Optimized-brightgreen.svg)

본 프로젝트는 고에너지 물리 실험(HEP)의 FADC 및 SADC 장비로부터 획득한 커스텀 C++ 객체 기반의 원시 데이터(RAW ROOT)를 **순수 파이썬(Pure Python) 환경에서 초고속으로 파싱**하여, 범용적이고 분석하기 쉬운 **Flat ROOT 형식(PRD)**으로 변환하는 초고속 데이터 전처리 파이프라인입니다.

## 🌟 주요 기능 및 특징 (Key Features)

1. 🛡️ **C++ Dictionary 종속성 완전 탈출 (Zero C++ Dependency)**
   - `libRawObjs.so`와 같은 전용 C++ 라이브러리(Dictionary) 없이 구동됩니다.
   - Uproot의 해석 엔진이 뻗어버리는 복잡한 커스텀 객체(`ArrayS`)를 **바스켓(Basket) 단위 바이트(Byte) 레벨에서 직접 역설계(Reverse Engineering)하여 100% 에러 없이 안전하게 디코딩**하는 지능형 자체 파서를 탑재하고 있습니다.
2. ⚡ **압도적인 고속 벡터 연산 (Ultra-fast Vectorization)**
   - 파이썬의 느린 `for` 루프 병목을 완전히 배제하고, `Numpy`와 `Awkward Array`의 C-백엔드 벡터 연산을 극한으로 활용합니다.
   - 단일 코어 환경에서 **544만 개 이상의 이벤트를 단 8분(초당 ~11,000 이벤트 처리)** 만에 매칭 및 파싱합니다.
3. 🎯 **정밀한 FADC-SADC 동기화 (High-precision Matching)**
   - 이벤트 트리거 번호(`TrgNum`)를 기준으로 글로벌 맵(Global Map)을 구축하고, 이진 탐색(Binary Search) 알고리즘을 적용하여 서브런(Subrun) 경계를 넘나드는 데이터를 완벽하게 페어링합니다.
4. 📦 **분석 친화적 Flat 구조 (Universal Flat Schema)**
   - 변환된 `PRD` 파일은 복잡한 중첩 레코드 대신 고정 길이 다차원 배열(`int32_t[4]`) 및 `std::vector` 호환 가변 배열 형식으로 평탄화(Flatten)되어 저장됩니다.
   - Numpy, Pandas, RDataFrame, Machine Learning 툴 등 어떤 환경에서도 즉시 로드 가능합니다.
5. 🗜️ **초고속 LZ4 압축 알고리즘 (I/O Optimization)**
   - 디스크 용량을 쥐어짜는 구형 Zlib 대신, 데이터 로딩(Read) 속도를 수십 배 향상시키는 `LZ4` 압축 알고리즘을 채택하여 분석 단계에서의 I/O 병목을 제거했습니다.

---

## 🛠️ 요구 환경 (Prerequisites)

순수 파이썬 패키지만으로 동작하며 PyROOT 설치가 필요하지 않습니다.
```bash
# Python 3.8 이상 권장
pip install numpy awkward uproot tqdm

```

---

## 📂 디렉토리 구조 (Directory Structure)

스크립트가 정상적으로 작동하기 위해 아래의 디렉토리 구조를 권장합니다. (경로는 소스 코드 내 `BASE_DIR` 변수에서 수정 가능합니다.)

```text
Cern_root/
 ├── RAW/                 # FADC 및 SADC 원본 ROOT 파일 모음
 │    └── 004183/
 │         ├── FADC_004183.root.00000 ~ 
 │         └── SADC_004183.root.00000 ~ 
 ├── TCBLOG/              # 장비 세팅값 텍스트 로그 파일
 │    └── TCB_004183.log (또는 .log.gz)
 ├── PRD/                 # 스크립트 실행 시 자동 생성되는 출력 디렉토리
 │    └── 004183/
 │         └── PRD_004183.00000.root ~
 └── rene_processing.py          # 메인 변환 실행 스크립트

```

---

## 🚀 사용 방법 (Usage)

1. 스크립트 상단의 사용자 환경 설정 변수를 세팅합니다.

```python
RUN_NUM = 4183
BASE_DIR = "/home/kds/mywork/Cern_root"

```

2. 스크립트를 실행합니다.

```bash
python target.py

```

3. 프로그레스 바를 통해 진행률이 실시간으로 표시되며, 완료 시 요약 리포트가 출력됩니다.

---

## 📊 출력 데이터 스키마 (PRD ROOT Schema)

생성된 `PRD` 파일에는 메타데이터를 담은 **`Run` TTree**와 메인 데이터인 **`Event` TTree**가 분리되어 저장됩니다.

### 1. `Run` Tree (파일당 1개 Entry)

해당 런(Run) 전체에 공통으로 적용된 하드웨어 및 소프트웨어 셋업 상수입니다.

| Branch Name | Type | Description |
| --- | --- | --- |
| `RunNumber` | `uint32_t` | 런 번호 |
| `nF` / `nS` | `int32_t` | FADC / SADC 활성 채널 개수 |
| `F_PmtID` / `S_PmtID` | `int32_t[]` | 각 채널에 매핑된 PMT ID 배열 |
| `F_THR`, `F_DLY` | `int32_t[]` | FADC 임계값(Threshold) 및 딜레이(Delay) |
| `S_THR`, `S_DLY` | `int32_t[]` | SADC 임계값 및 딜레이 |

### 2. `Event` Tree (이벤트별 1개 Entry)

트리거 단위로 결합된 채널별 측정 데이터 및 원시 파형(Waveform)입니다.

| Branch Name | Type | Description |
| --- | --- | --- |
| `TrgNum` | `uint32_t` | 글로벌 트리거 번호 (매칭 기준) |
| `EventType` | `uint32_t` | 임계값 초과 이벤트 여부 플래그 (1=FADC, 2=SADC, 3=Both) |
| `TCBTRGTime` | `double` | TCB 기준 트리거 발생 시간 |
| `F_Pedestal` | `int16_t[4]` | FADC 채널별 베이스라인 잡음(Pedestal) |
| `S_ADC` | `int32_t[4]` | SADC 채널별 적분 전하량 (ADC) |
| `S_PeakTime` | `double[4]` | SADC 채널별 피크 발생 시간 |
| `nF_Waveform_X` | `int32_t` | X번 채널 파형의 샘플 길이 |
| **`F_Waveform_X`** | **`uint16_t[]`** | **[가변 배열] X번 채널 FADC 파형 시계열 데이터 (`std::vector` 완벽 호환)** |

---

## 📈 벤치마크 및 성능 요약 (Performance Benchmark)

* 🖥️ **테스트 환경**: Single-core Python (`awkward`, `uproot`)
* 📦 **데이터 규모**: Run 4183 (12개 Subrun 파일, 원본 `7.03 GB`)
* ⏱️ **처리 시간**: **`8.33 분 (500.0 초)`**
* 🎯 **처리량**: **`5,440,123 개 이벤트`** 무결성 매칭 및 역설계 변환 성공 (유실률 0%)
* 💾 **디스크 최적화**: 분석 단계의 읽기 속도를 극한으로 끌어올리기 위한 **LZ4 고속 압축** 및 다차원 배열 정렬 뼈대(Padding) 생성으로 인해 파일 용량은 약 `9.73 GB`로 증가함 (I/O 속도를 위한 의도된 최적화).

---

## 🧩 트러블슈팅 및 기술적 배경 (Technical Notes)

**Q: Uproot의 `UnknownInterpretation` 에러를 어떻게 해결했나요?**

> RAW 파일에 기록된 `fColl.ArrayS` 객체는 C++ Dictionary가 없어 Uproot의 자동 파싱이 불가능합니다. 본 코드의 `parse_ArrayS_branch` 함수는 ROOT Basket 바이트 구조를 직접 열어 Object-wise 및 Member-wise 직렬화 마스킹 헤더(`0x40000000`)를 지능적으로 동적 판별하는 독자적인 바이트 디코더입니다. 버퍼 경계 잘림(Buffer cutoff) 시에도 데이터를 짝수(`uint16`) 바이트로 자동 보정하여 에러를 원천 봉쇄합니다.

**Q: Awkward Array 변환 중 `TypeError: option[var * uint16]`는 왜 발생했나요?**

> 채널 간 병합 시 빈공간 패딩(`pad_none`)을 거친 배열은 결측치(None) 이력이 남아 ROOT의 `std::vector`로 직접 저장이 불가능합니다. 이를 해결하기 위해 배열을 1차원 Numpy 평면으로 완전히 으깬(`flatten`) 뒤 본래 길이대로 재조립(`unflatten`)하여, 옵션 타입 꼬리표를 강제로 떼어내는 '메모리 세탁(Memory Cleansing)' 기법을 적용했습니다.

---

*Developed and optimized for High-Energy Physics Data Processing Pipeline.*

```
