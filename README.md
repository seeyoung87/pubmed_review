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
- `CONFIG_PATH`: 설정 파일 경로 (기본 `config.yaml`)

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
`sheet_name_template`로 지정한 템플릿에 넣습니다. 예) `'Radiology NLP'의 정리 시트`.
