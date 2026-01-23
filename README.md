# PubMed Review Automation

NCBI PubMed API로 최신 논문을 검색하고 LLM 요약을 만든 뒤 Google Sheets에 저장하는 자동화 파이프라인입니다.

## 전체 흐름
1. PubMed 쿼리로 최근 등록 논문을 조회합니다.
2. 논문 메타데이터/초록을 가져옵니다.
3. LLM으로 Novelty/요약을 생성해 필터링합니다.
4. 결과를 Google Sheets에 기록합니다.

## 필수 환경 변수
- `PUBMED_EMAIL`: NCBI 연락용 이메일 주소 (Entrez 필수).
- `OPENAI_API_KEY`: OpenAI API 키.
- `GOOGLE_SERVICE_ACCOUNT_JSON`: 서비스 계정 JSON 문자열.
- `SPREADSHEET_ID`: (선택) Google Sheets 파일 ID. 설정 시 `config.yaml`의 `sheets.spreadsheet_id`보다 우선합니다.
- `CONFIG_PATH`: (선택) 설정 파일 경로 (기본 `config.yaml`).

## 설정 파일(`config.yaml`) 핵심
### 1) 단일 쿼리 + 단일 시트
```yaml
pubmed:
  email: "your_email@example.com"
  search_query: '("Radiology") AND ("large language model")'
  sheet_name: "Radiology NLP"
```

### 2) 복수 쿼리 + 복수 시트
```yaml
pubmed:
  searches:
    - query: '서치쿼리1'
      sheet_name: '시트네임1'
    - query: '서치쿼리2'
      sheet_name: '시트네임2'
```

### 3) 검색 기간 동기화
`pubmed.reldate`를 비워두면 `workflow.schedule_days` 값을 사용합니다.
GitHub Actions 실행 주기와 검색 기간을 동일하게 맞춰야 누락이 없습니다.

### 4) Google Sheets 출력 탭
`pubmed.sheet_name` → `sheets.sheet_name` → 기본 `PubMed` 순으로 사용됩니다.

## Google Sheets 준비
1. Google Sheets에서 스프레드시트를 생성합니다.
2. URL의 `/d/<ID>/edit` 구간에서 `<ID>`를 `sheets.spreadsheet_id`에 넣거나 `SPREADSHEET_ID`로 설정합니다.
3. 서비스 계정 이메일을 스프레드시트에 **편집자**로 공유합니다.

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

## LLM 출력 스키마
`llm.novelty_schema`, `llm.summary_schema`에 JSON Schema를 정의해 출력 포맷을 제어합니다.
