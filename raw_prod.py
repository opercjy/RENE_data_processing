import os
import glob
import gzip
import struct
import time
import numpy as np
import awkward as ak
import uproot
from tqdm.notebook import tqdm

# =========================================================================
# ⚙️ 1. 사용자 설정 (Run Number 및 기본 경로)
# =========================================================================
RUN_NUM = 4183
BASE_DIR = "/home/kds/mywork/Cern_root"
RAW_DIR = f"{BASE_DIR}/RAW/{RUN_NUM:06d}"
PRD_DIR = f"{BASE_DIR}/PRD/{RUN_NUM:06d}"
TCBLOG_DIR = f"{BASE_DIR}/TCBLOG" 
os.makedirs(PRD_DIR, exist_ok=True)

print("="*60)
print(f"🚀 본격적인 ROOT 데이터 초고속 플랫 변환 시작 (Run {RUN_NUM})")
print("="*60)
print(f"📁 입력(RAW) 경로: {RAW_DIR}")
print(f"📁 출력(PRD) 경로: {PRD_DIR}")

# 전체 실행 시간 및 용량 측정용 변수
start_time = time.time()
raw_total_bytes = 0
prd_total_bytes = 0
total_events_processed = 0

# =========================================================================
# 🛠️ 2. 파서 클래스 및 헬퍼 함수
# =========================================================================
class RunInfo:
    def __init__(self, runNumber, dictFADC, dictSADC):
        nF, nS = 0, 0
        for key, value in dictFADC.items():
            if key == 'NADC': continue
            setattr(self, f"F_{key}", np.array(value, dtype=np.int32))
            nF = len(value)
        for key, value in dictSADC.items():
            if key == 'NADC': continue
            setattr(self, f"S_{key}", np.array(value, dtype=np.int32))
            nS = len(value)
        self.runNumber = np.array([runNumber], dtype=np.uint32)
        self.nF = np.array([nF], dtype=np.int32)
        self.nS = np.array([nS], dtype=np.int32)

class TCBLogReader:
    def __init__(self, runNumber, log_dir):
        fName = f"{log_dir}/TCB_{runNumber:06d}.log"
        self.lines = []
        if os.path.exists(fName):
            with open(fName, 'rt') as f: self.lines = f.readlines()
        elif os.path.exists(fName+".gz"):
            with gzip.open(fName+".gz", 'rt') as f: self.lines = f.readlines()
        else:
            print(f"⚠️ 경고: TCB 로그 파일({fName})을 찾을 수 없어 기본값으로 초기화합니다.")

    def ExtractWJ(self):
        infoFADC, infoSADC = {}, {}
        for line in self.lines:
            line = line.strip().split(' ', 2)
            if len(line) < 3 or line[0] != "WJ" or '=' not in line[-1]: continue
            varName, values = line[-1].split('=', 1)
            varName = varName.strip().split()[0]
            values = [int(x) for x in values.split()]
            if line[1] == 'FADC': infoFADC[varName] = values
            elif line[1] == 'SADC': infoSADC[varName] = values
        return infoFADC, infoSADC

def get_b(tree, keyword, lib="ak"):
    for k in tree.keys():
        if k.endswith(keyword) or k == keyword:
            arr = tree[k].array(library=lib)
            if lib == "ak" and hasattr(arr, "fields") and "fArray" in arr.fields: return arr["fArray"]
            return arr
    raise KeyError(f"'{keyword}' 브랜치를 찾을 수 없습니다.")

def get_branch_obj(tree, keyword):
    for k in tree.keys():
        if k.endswith(keyword) or k == keyword: return tree[k]
    raise KeyError(f"'{keyword}' 객체를 찾을 수 없습니다.")

def to_rectangular(arr, length, fill_val=0):
    length = int(length)
    padded = ak.pad_none(arr, length, axis=1, clip=True)
    return ak.to_numpy(ak.fill_none(padded, fill_val))

def format_size(bytes_size):
    """바이트 크기를 보기 좋은 단위로 변환"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0: return f"{bytes_size:.2f} {unit}"
        bytes_size /= 1024.0

# -------------------------------------------------------------------------
# 🛡️ ArrayS 무결점 역설계 바이트 디코더
# -------------------------------------------------------------------------
def parse_ArrayS_branch(branch, ch_counts):
    all_waves, lengths = [], []
    ch_counts_np = ak.to_numpy(ch_counts).astype(np.int32)
    global_event_idx, total_events = 0, len(ch_counts_np)
    
    for i in range(branch.num_baskets):
        basket = branch.basket(i)
        data = basket.data
        if hasattr(data, 'tobytes'): data = data.tobytes()
        else: data = bytes(data)
            
        border = getattr(basket, 'border', len(data))
        data = data[:border]
        offsets = basket.byte_offsets
        if offsets is None:
            step = border // basket.num_entries if basket.num_entries > 0 else border
            offsets = [j * step for j in range(basket.num_entries + 1)]
            
        for evt in range(basket.num_entries):
            if global_event_idx >= total_events: break
            
            start, stop = offsets[evt], offsets[evt+1] if evt + 1 < len(offsets) else border
            start, stop = min(start, border), min(stop, border)
            evt_data, evt_data_len, cursor = data[start:stop], stop - start, 0
            
            M_channels = ch_counts_np[global_event_idx]
            global_event_idx += 1
            if M_channels == 0: continue
            
            if evt_data_len < 4: 
                for _ in range(M_channels):
                    all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)
                continue
                
            val = struct.unpack_from(">I", evt_data, cursor)[0]
            if val & 0x40000000:
                for _ in range(M_channels):
                    if cursor + 4 <= evt_data_len:
                        val = struct.unpack_from(">I", evt_data, cursor)[0]
                        cursor += 6 if val & 0x40000000 else 4
                    if cursor + 4 <= evt_data_len:
                        fN = struct.unpack_from(">i", evt_data, cursor)[0]; cursor += 4
                    else: fN = 0
                        
                    flag = evt_data[cursor] if cursor < evt_data_len else 0
                    cursor += 1 if cursor < evt_data_len else 0
                        
                    if flag == 1 and 0 < fN < 1000000:
                        byte_len = fN * 2
                        if cursor + byte_len <= evt_data_len:
                            arr = np.frombuffer(evt_data[cursor:cursor+byte_len], dtype='>u2').astype(np.uint16)
                            cursor += byte_len
                            all_waves.append(arr); lengths.append(fN)
                        else:
                            rem = ((evt_data_len - cursor) // 2) * 2
                            if rem > 0:
                                arr = np.frombuffer(evt_data[cursor:cursor+rem], dtype='>u2').astype(np.uint16)
                                cursor += rem
                                all_waves.append(arr); lengths.append(rem // 2)
                            else:
                                all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)
                    else:
                        all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)
            else:
                fN_list = []
                for _ in range(M_channels):
                    if cursor + 4 <= evt_data_len:
                        fN = struct.unpack_from(">i", evt_data, cursor)[0]; cursor += 4
                    else: fN = 0
                    fN_list.append(fN)
                    
                for fN in fN_list:
                    flag = evt_data[cursor] if cursor < evt_data_len else 0
                    cursor += 1 if cursor < evt_data_len else 0
                        
                    if flag == 1 and 0 < fN < 1000000:
                        byte_len = fN * 2
                        if cursor + byte_len <= evt_data_len:
                            arr = np.frombuffer(evt_data[cursor:cursor+byte_len], dtype='>u2').astype(np.uint16)
                            cursor += byte_len
                            all_waves.append(arr); lengths.append(fN)
                        else:
                            rem = ((evt_data_len - cursor) // 2) * 2
                            if rem > 0:
                                arr = np.frombuffer(evt_data[cursor:cursor+rem], dtype='>u2').astype(np.uint16)
                                cursor += rem
                                all_waves.append(arr); lengths.append(rem // 2)
                            else:
                                all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)
                    else:
                        all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)

    while global_event_idx < total_events:
        for _ in range(ch_counts_np[global_event_idx]):
            all_waves.append(np.array([], dtype=np.uint16)); lengths.append(0)
        global_event_idx += 1

    if len(all_waves) > 0:
        flat_waves = np.concatenate(all_waves)
        jagged_waves = ak.unflatten(flat_waves, lengths)
    else:
        jagged_waves = ak.Array([])
    return ak.unflatten(jagged_waves, ch_counts_np)

# =========================================================================
# 🚀 3. 메인 파이프라인 (인덱싱 및 플랫 직렬화)
# =========================================================================
subruns = []
fNamesFADC, fNamesSADC = {}, {}
for fNameFADC in sorted(glob.glob(f"{RAW_DIR}/FADC_{RUN_NUM:06d}.root.*")):
    subrun = fNameFADC.rsplit('.', 1)[-1]
    fNameSADC = f"{RAW_DIR}/SADC_{RUN_NUM:06d}.root.{subrun}"
    if os.path.exists(fNameSADC):
        subruns.append(subrun)
        fNamesFADC[subrun] = fNameFADC
        fNamesSADC[subrun] = fNameSADC
        raw_total_bytes += os.path.getsize(fNameFADC) + os.path.getsize(fNameSADC)

if not subruns:
    raise RuntimeError("에러: 처리 가능한 ROOT 파일이 존재하지 않습니다!")

infoFADC, infoSADC = TCBLogReader(RUN_NUM, TCBLOG_DIR).ExtractWJ()
if not infoFADC: infoFADC = {'PmtID': [0]*4, 'DLY': [0]*4, 'THR': [0]*4, 'RL': [0]*4}
if not infoSADC: infoSADC = {'PmtID': [0]*4, 'DLY': [0]*4, 'THR': [0]*4, 'GW': [0]*4}

runInfo = RunInfo(RUN_NUM, infoFADC, infoSADC)
nF_int, nS_int = int(runInfo.nF[0]), int(runInfo.nS[0])
F_PID = np.asarray(getattr(runInfo, 'F_PID', getattr(runInfo, 'F_PmtID', np.zeros(nF_int))), dtype=np.int32)
S_PID = np.asarray(getattr(runInfo, 'S_PID', getattr(runInfo, 'S_PmtID', np.zeros(nS_int))), dtype=np.int32)

sadc_trg_global, sadc_file_idx, sadc_entry_idx = [], [], []

print(f"\n✅ 총 {len(subruns)}개의 Subrun 파일 짝을 로드했습니다. (원본 누적 용량: {format_size(raw_total_bytes)})")
print("⏳ SADC 트리거 글로벌 맵 구축 중...")

for i_sub, sub in enumerate(tqdm(subruns, desc="Indexing SADC")):
    with uproot.open(fNamesSADC[sub]) as f:
        trgs = get_b(f["AbsEvent"], "EventInfo/fTrgNum", "np")
        sadc_trg_global.append(trgs)
        sadc_file_idx.append(np.full(len(trgs), i_sub, dtype=np.int32))
        sadc_entry_idx.append(np.arange(len(trgs), dtype=np.int32))

sadc_trg_global = np.concatenate(sadc_trg_global)
sadc_file_idx = np.concatenate(sadc_file_idx)
sadc_entry_idx = np.concatenate(sadc_entry_idx)

sort_mask = np.argsort(sadc_trg_global)
sadc_trg_global = sadc_trg_global[sort_mask]
sadc_file_idx = sadc_file_idx[sort_mask]
sadc_entry_idx = sadc_entry_idx[sort_mask]

# [미션 2] 트리거 인덱싱 요약 정보
print(f" └─ 📊 매칭 대기 중인 전체 SADC 이벤트 수 : {len(sadc_trg_global):,} 개")
print(f" └─ 🔢 트리거 번호 탐색 범위 : {sadc_trg_global[0]} ~ {sadc_trg_global[-1]}")

first_file_saved = False # 딕셔너리 구조 출력을 위한 스위치
print("\n🚀 본격적인 데이터 변환 및 플랫 파일 작성 시작!")

for subrun in tqdm(subruns, desc=f"Converting Run {RUN_NUM}"):
    fNameFADC = fNamesFADC[subrun]
    out_filename = f"{PRD_DIR}/PRD_{RUN_NUM:06d}.{subrun}.root"

    # [A] FADC 로드
    with uproot.open(fNameFADC) as f_fadc:
        tree_f = f_fadc["AbsEvent"]
        f_trg = get_b(tree_f, "EventInfo/fTrgNum", "np")
        f_tcb = get_b(tree_f, "EventInfo/fTCBTrgTime", "np")
        f_id_raw  = get_b(tree_f, "fColl.fID", "ak")
        f_bit_raw = get_b(tree_f, "fColl.fBit", "ak")
        f_ped_raw = get_b(tree_f, "fColl.fPedestal", "ak")
        
        ch_counts = ak.num(f_id_raw)
        f_wav_raw = parse_ArrayS_branch(get_branch_obj(tree_f, "fColl.ArrayS"), ch_counts)
        f_ndp_raw = ak.num(f_wav_raw, axis=-1)

    # [B] SADC 고속 이진 탐색 매칭
    indices = np.searchsorted(sadc_trg_global, f_trg)
    valid = indices < len(sadc_trg_global)
    matched = valid.copy()
    matched[valid] = (sadc_trg_global[indices[valid]] == f_trg[valid])
    
    m_fadc_idx = np.where(matched)[0]
    M = len(m_fadc_idx)
    total_events_processed += M
    if M == 0: continue

    m_sadc_idx_global = indices[matched]
    m_sadc_file_idx = sadc_file_idx[m_sadc_idx_global]
    m_sadc_ent_idx = sadc_entry_idx[m_sadc_idx_global]

    # [C] 매칭된 SADC 데이터 추출
    s_tcb_chunks, s_id_chunks, s_bit_chunks, s_adc_chunks, s_time_chunks = [], [], [], [], []
    for file_idx in sorted(np.unique(m_sadc_file_idx)):
        mask = (m_sadc_file_idx == file_idx)
        entries = m_sadc_ent_idx[mask]
        
        with uproot.open(fNamesSADC[subruns[file_idx]]) as fs:
            ts = fs["AbsEvent"]
            s_tcb_chunks.append(get_b(ts, "EventInfo/fTCBTrgTime", "np")[entries])
            s_id_chunks.append(get_b(ts, "fColl.fID", "ak")[entries])
            s_bit_chunks.append(get_b(ts, "fColl.fBit", "ak")[entries])
            s_adc_chunks.append(get_b(ts, "fColl.fADC", "ak")[entries])
            s_time_chunks.append(get_b(ts, "fColl.fTime", "ak")[entries])
            
    s_tcb = np.concatenate(s_tcb_chunks)
    s_id   = to_rectangular(ak.concatenate(s_id_chunks), nS_int)
    s_bit  = to_rectangular(ak.concatenate(s_bit_chunks), nS_int)
    s_adc  = to_rectangular(ak.concatenate(s_adc_chunks), nS_int)
    s_time = to_rectangular(ak.concatenate(s_time_chunks), nS_int)

    # [D] Numpy & Awkward 초고속 벡터 연산
    f_tcb_cut = f_tcb[m_fadc_idx]
    f_id  = to_rectangular(f_id_raw[m_fadc_idx], nF_int)
    f_bit = to_rectangular(f_bit_raw[m_fadc_idx], nF_int)
    f_ped = to_rectangular(f_ped_raw[m_fadc_idx], nF_int)
    f_ndp = to_rectangular(f_ndp_raw[m_fadc_idx], nF_int)
    f_wav = f_wav_raw[m_fadc_idx]
    
    f_thr_2d = np.tile(np.array(runInfo.F_THR), (M, 1))
    f_dly_2d = np.tile(np.array(runInfo.F_DLY), (M, 1))
    f_wavestart = f_tcb_cut[:, np.newaxis] - f_dly_2d
    f_threshold = f_ped + np.array(runInfo.F_THR)
    
    has_fadc_over = np.zeros(M, dtype=bool)
    clean_wavs = []
    
    for iCH in range(nF_int):
        # 파이썬 Numpy 평탄화 기법을 통한 완전한 배열 Option 세탁
        has_ch = ak.to_numpy(ak.num(f_wav) > iCH)
        valid_wavs = f_wav[has_ch, iCH]
        
        lengths = np.zeros(M, dtype=np.int32)
        lengths[has_ch] = ak.to_numpy(ak.num(valid_wavs))
        
        flat = ak.to_numpy(ak.flatten(valid_wavs)).astype(np.uint16) if len(valid_wavs) > 0 else np.array([], dtype=np.uint16)
        clean_ch = ak.unflatten(flat, lengths)
        clean_wavs.append(clean_ch)
        
        is_over = clean_ch > f_threshold[:, iCH]
        has_fadc_over |= ak.to_numpy(ak.any(is_over, axis=-1))

    s_thr_2d = np.tile(np.array(runInfo.S_THR), (M, 1))
    s_dly_2d = np.tile(np.array(runInfo.S_DLY), (M, 1))
    s_wavestart = s_tcb[:, np.newaxis] - s_dly_2d
    s_peaktime = np.where(s_time - s_wavestart < 0, -99.0, s_time - s_wavestart)
    
    is_s_over = s_adc > np.array(runInfo.S_THR)
    has_sadc_over = np.sum(is_s_over, axis=-1) > 0
    event_type = np.where(has_fadc_over, 1, 0) + np.where(has_sadc_over, 2, 0)

    # [E] LZ4 압축을 적용하여 플랫 ROOT 트리 쓰기
    with uproot.recreate(out_filename, compression=uproot.LZ4(4)) as f_out:
        run_dict = {
            "RunNumber": np.array([RUN_NUM], dtype=np.uint32),
            "nF": np.array([nF_int], dtype=np.int32), "F_PmtID": np.array([F_PID], dtype=np.int32),
            "F_DLY": np.array([runInfo.F_DLY], dtype=np.int32), "F_THR": np.array([runInfo.F_THR], dtype=np.int32), "F_RL": np.array([getattr(runInfo, 'F_RL', np.zeros(nF_int, dtype=np.int32))], dtype=np.int32),
            "nS": np.array([nS_int], dtype=np.int32), "S_PmtID": np.array([S_PID], dtype=np.int32),
            "S_DLY": np.array([runInfo.S_DLY], dtype=np.int32), "S_THR": np.array([runInfo.S_THR], dtype=np.int32), "S_GW": np.array([getattr(runInfo, 'S_GW', np.zeros(nS_int, dtype=np.int32))], dtype=np.int32)
        }

        event_dict = {
            "TrgNum": np.asarray(f_trg[m_fadc_idx], dtype=np.uint32),
            "EventType": np.asarray(event_type, dtype=np.uint32),
            "TCBTRGTime": np.asarray(f_tcb_cut, dtype=np.float64),
            "nCH_FADC": np.full(M, nF_int, dtype=np.int32),
            "F_PmtID": f_id.astype(np.int32), "F_THR": f_thr_2d.astype(np.uint16),
            "F_Triggered": f_bit.astype(np.int32), "F_WaveStartTime": f_wavestart.astype(np.float64),
            "F_Pedestal": f_ped.astype(np.int16), "F_NDP": f_ndp.astype(np.int32),
            "nCH_SADC": np.full(M, nS_int, dtype=np.int32),
            "S_PmtID": s_id.astype(np.int32), "S_THR": s_thr_2d.astype(np.uint16),
            "S_Triggered": s_bit.astype(np.int32), "S_WaveStartTime": s_wavestart.astype(np.float64),
            "S_PeakTime": s_peaktime.astype(np.float64), "S_ADC": s_adc.astype(np.int32)
        }
        
        # 완벽히 정제된 파형 데이터 결합
        for iCH in range(nF_int):
            event_dict[f"F_Waveform_{iCH}"] = clean_wavs[iCH]

        # [미션 2] 첫 번째 서브런일 때 딕셔너리 구조 출력
        if not first_file_saved:
            print("\n" + "-"*65)
            print("📋 [생성된 PRD 트리의 딕셔너리 스키마 구조 요약]")
            print(f" ▶ Run Tree : {list(run_dict.keys())}")
            print(f" ▶ Event Tree : {list(event_dict.keys())[:8]} ... (이하 파형 채널 생략)")
            print("-" * 65 + "\n")
            first_file_saved = True

        f_out["Run"] = run_dict
        f_out["Event"] = event_dict
    
    prd_total_bytes += os.path.getsize(out_filename)

# =========================================================================
# 📊 4. [미션 2] 최종 실행 시간 및 압축률 리포트
# =========================================================================
elapsed_time = time.time() - start_time
comp_ratio = (prd_total_bytes / raw_total_bytes) * 100 if raw_total_bytes > 0 else 0

print("\n" + "="*60)
print("🏆 [작업 요약 리포트] 완벽한 플랫 루트 데이터가 생성되었습니다!")
print("="*60)
print(f"⏱️ 총 소요 시간 : {elapsed_time/60:.2f} 분 ({elapsed_time:.1f} 초)")
print(f"🎯 처리된 이벤트: {total_events_processed:,} 개 매칭 및 변환 성공")
print(f"💾 원본 파일 용량 (RAW): {format_size(raw_total_bytes)}")
print(f"📦 생성된 플랫 용량 (PRD): {format_size(prd_total_bytes)}")
print(f"🔥 파일 압축 효율 : {comp_ratio:.1f}% 수준으로 다이어트 성공! (기존 대비 약 {100-comp_ratio:.1f}% 절약)")
print("="*60 + "\n")
