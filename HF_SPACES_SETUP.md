# Hugging Face Spaces 자동 동기화 설정

## 🔄 방법 1: HF Spaces에서 GitHub 직접 연결 (권장)

### 단계별 설정:

1. **HF Spaces 설정 페이지 접속**
   - https://huggingface.co/spaces/Flamehaven/arr-medic-cyp3a4-demo/settings

2. **Repository 설정**
   - "Repository" 섹션에서 "Connect to GitHub" 클릭
   - GitHub repository 선택: `Flamehaven/Arr-Medic-CYP3A4`
   - Branch: `master`
   - Path: `demo/` (데모 파일들이 있는 폴더)

3. **자동 동기화 활성화**
   - "Auto-sync" 옵션 체크
   - "Sync on push" 활성화

### 결과:
- ✅ GitHub `demo/` 폴더 수정 → 자동으로 HF Spaces 업데이트
- ✅ Real-time 동기화
- ✅ 별도 설정 불필요

---

## 🤖 방법 2: GitHub Actions 파이프라인

### 필요한 설정:

1. **GitHub Secrets 추가**
   ```
   Repository Settings > Secrets > Actions
   - HF_TOKEN: Hugging Face API 토큰
   ```

2. **HF 토큰 생성**
   - https://huggingface.co/settings/tokens
   - "Write" 권한으로 토큰 생성

3. **파이프라인 활성화**
   - 이미 `.github/workflows/deploy-hf-spaces.yml` 생성됨
   - `demo/` 폴더 수정시 자동 배포

### 트리거 조건:
- `demo/**` 파일 수정
- `backend/predictor.py` 수정
- `requirements.txt` 수정
- 수동 실행 가능 (workflow_dispatch)

---

## 🎯 권장사항

**방법 1 (HF Spaces 직접 연결)**을 먼저 시도하세요:
- 설정이 더 간단함
- HF에서 공식 지원
- 실시간 동기화
- GitHub Actions 리소스 불필요

**방법 2**는 더 고급 제어가 필요한 경우:
- 커스텀 배포 로직
- 다중 환경 배포
- 복잡한 빌드 과정

---

## 🔧 현재 파일 구조

```
ARR-medic-cyp3a4-opensource/
├── demo/
│   ├── app.py              # 멀티링구얼 Gradio 앱
│   ├── requirements.txt    # HF Spaces 호환 의존성
│   └── Dockerfile         # 컨테이너 설정
├── backend/
│   └── predictor.py       # 예측 엔진
└── .github/workflows/
    └── deploy-hf-spaces.yml # 자동 배포 파이프라인
```

## ✅ 테스트 방법

1. `demo/app.py`에 작은 수정 (예: 주석 추가)
2. GitHub에 push
3. HF Spaces가 자동으로 업데이트되는지 확인

---

**추천: 방법 1을 먼저 시도하고, 필요시 방법 2로 업그레이드하세요!**