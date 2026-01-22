# PubMed Review Automation

NCBI PubMed 메일을 읽어서 논문을 선별하고 LLM 요약을 만든 뒤 Google Sheets에 저장하는 자동화 스크립트입니다.

## 구성
- `config.yaml`: 모델/프롬프트, High IF 리스트, 시트 설정
- `pubmed_review/main.py`: 전체 실행 엔트리포인트
- GitHub Actions: 3일마다 실행

## 환경 변수
- `EMAIL_USER`: 메일 계정
- `EMAIL_PASS`: 메일 앱 비밀번호
- `OPENAI_API_KEY`: OpenAI API 키
- `GOOGLE_SERVICE_ACCOUNT_JSON`: 서비스 계정 JSON 문자열 (권장)
  - 또는 `GOOGLE_SERVICE_ACCOUNT_FILE` 경로 사용 가능
- `SPREADSHEET_ID`: (선택) Google Sheets 파일 ID. 설정하면 `config.yaml`의 `sheets.spreadsheet_id`보다 우선합니다.
- `CONFIG_PATH`: 설정 파일 경로 (기본 `config.yaml`)

## Google Sheets 설정 방법
1. Google Sheets에서 스프레드시트를 생성합니다.
2. 주소창 URL에서 `/d/<ID>/edit` 사이의 `<ID>`를 `config.yaml`의 `sheets.spreadsheet_id`에 넣거나 `SPREADSHEET_ID` 환경 변수로 설정합니다.
3. 서비스 계정 이메일을 스프레드시트에 **편집자**로 공유합니다.
4. 메일 제목에서 추출된 검색명이 없다면 `sheets.sheet_name` (없으면 기본 `PubMed`) 탭에 기록됩니다.

## 로컬 실행
```bash
python -m pubmed_review.main
```

## Google Sheets 컬럼
1. 날짜
2. PMID
3. 제목
4. 저널
5. 발행일
6. DOI
7. 선별 기준 (High IF, Novelty)
8. Novelty 근거
9. 요약
10. 우수성

## 시트 이름
메일 제목이 `What's new for 'Radiology NLP' in PubMed` 형태라면 따옴표 안의 이름을 추출해
해당 이름을 시트 탭 이름으로 그대로 사용합니다. (예: `Radiology NLP`)
