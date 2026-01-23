# PubMed Review Automation

NCBI PubMed API로 최신 논문을 검색하고 LLM 요약을 만든 뒤 Google Sheets에 저장하는 자동화 스크립트입니다.

## 구성
- `config.yaml`: PubMed 검색 설정, 모델/프롬프트, High IF 리스트, 시트 설정
- `pubmed_review/main.py`: 전체 실행 엔트리포인트
- GitHub Actions: 3일마다 실행

## 환경 변수
- `PUBMED_EMAIL`: NCBI 연락용 이메일 주소 (Entrez 필수). `config.yaml`의 `pubmed.email`보다 우선합니다.
- `PUBMED_SEARCH_QUERY`: PubMed 검색 쿼리. `config.yaml`의 `pubmed.search_query`보다 우선합니다.
- `OPENAI_API_KEY`: OpenAI API 키
- `GOOGLE_SERVICE_ACCOUNT_JSON`: 서비스 계정 JSON 문자열 (권장)
  - 또는 `GOOGLE_SERVICE_ACCOUNT_FILE` 경로 사용 가능
- `SPREADSHEET_ID`: (선택) Google Sheets 파일 ID. 설정하면 `config.yaml`의 `sheets.spreadsheet_id`보다 우선합니다.
- `CONFIG_PATH`: 설정 파일 경로 (기본 `config.yaml`)

## Google Sheets 설정 방법
1. Google Sheets에서 스프레드시트를 생성합니다.
2. 주소창 URL에서 `/d/<ID>/edit` 사이의 `<ID>`를 `config.yaml`의 `sheets.spreadsheet_id`에 넣거나 `SPREADSHEET_ID` 환경 변수로 설정합니다.
3. 서비스 계정 이메일을 스프레드시트에 **편집자**로 공유합니다.
4. `pubmed.search_name`이 없다면 `sheets.sheet_name` (없으면 기본 `PubMed`) 탭에 기록됩니다.

## 로컬 실행
```bash
python -m pubmed_review.main
```

## 스케줄과 검색 기간 동기화
GitHub Actions의 실행 주기(`workflow.schedule_days`)와 PubMed 검색 기간(`pubmed.reldate`)을 동일하게 맞춰야
누락 없이 해당 기간의 논문을 처리할 수 있습니다. `pubmed.reldate`를 비워두면 `workflow.schedule_days` 값을 사용합니다.

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
`pubmed.search_name`을 설정하면 시트 탭 이름으로 그대로 사용합니다. (예: `Radiology NLP`)

## LLM 출력 스키마
`llm.novelty_schema`, `llm.summary_schema`에 JSON Schema를 정의해 출력 포맷을 제어합니다.
