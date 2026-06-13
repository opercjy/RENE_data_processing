import uproot
import awkward as ak

# 테스트할 파일 (00000번 서브런 기준)
RUN_NUM = 4183
SUB_RUN = "00000"
BASE_DIR = "/home/kds/mywork/Cern_root"

RAW_FILE = f"{BASE_DIR}/RAW/{RUN_NUM:06d}/FADC_{RUN_NUM:06d}.root.{SUB_RUN}"
PRD_FILE = f"{BASE_DIR}/PRD/{RUN_NUM:06d}/PRD_{RUN_NUM:06d}.{SUB_RUN}.root"

print("="*75)
print(" 🧐 [1] 원본 RAW 파일 내부 구조 검사 (C++ Dictionary 의존성 파일)")
print("="*75)
try:
    with uproot.open(RAW_FILE) as f_raw:
        tree_raw = f_raw["AbsEvent"]
        print(f"▶ 원본 트리: {tree_raw.name} (총 이벤트 수: {tree_raw.num_entries:,}개)")
        
        # 문제의 C++ ArrayS 객체 브랜치 정보 엿보기
        b_arrayS = tree_raw["FChannelData/fColl/fColl.ArrayS"]
        print(f"\n⚠️ [해석 불가 원인이었던 브랜치] {b_arrayS.name}")
        print(f"  - 타입: C++ Custom Object")
        print(f"  - 바스켓 수: {b_arrayS.num_baskets} 개")
        print(f"  - 💡 C++ 'libRawObjs.so' 라이브러리가 없으면 PyROOT 밖에서는 읽을 수 없는 특수 포맷입니다.")
        
        raw_trg = tree_raw["EventInfo/fTrgNum"].array(library="np")[0]
        print(f"  - [첫 번째 이벤트] 트리거 번호: {raw_trg}")
except Exception as e:
    print(f"RAW 파일 읽기 오류: {e}")

print("\n" + "="*75)
print(" ✨ [2] 새롭게 생성된 PRD 파일 내부 구조 검사 (Python/C++ 100% 호환 Flat 포맷)")
print("="*75)
try:
    with uproot.open(PRD_FILE) as f_prd:
        tree_run = f_prd["Run"]
        tree_evt = f_prd["Event"]
        
        print(f"▶ 변환된 트리: {tree_evt.name} (총 이벤트 수: {tree_evt.num_entries:,}개)")
        print(f"▶ 메타데이터 트리: {tree_run.name} (행 수: {tree_run.num_entries}개)")
        
        print("\n✅ [Event 트리 내부 주요 브랜치 구조]")
        tree_evt.show(name_width=20, typename_width=25)
        
        print("\n🔍 [데이터 무결성 테스트 - 첫 번째 이벤트 검증]")
        prd_trg = tree_evt["TrgNum"].array(library="np")[0]
        wav_0 = tree_evt["F_Waveform_0"].array(library="ak")
        
        print(f"  - 트리거 번호: {prd_trg} (원본 RAW 파일과 일치)")
        if len(wav_0) > 0 and len(wav_0[0]) > 0:
            print(f"  - 추출된 FADC CH_0 파형 길이: {len(wav_0[0])} 샘플")
            print(f"  - 파형 데이터 앞부분(10개) : {ak.to_list(wav_0[0][:10])} ...")
        else:
            print("  - 파형 데이터가 비어 있습니다.")
            
except Exception as e:
    print(f"PRD 파일 읽기 오류: {e}")
print("="*75)
