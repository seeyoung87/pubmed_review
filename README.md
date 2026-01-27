# PubMed Review Automation

## What It Does

```
PubMed 쿼리 검색 → Selected Journal? → Yes → AI 요약/강점 정리 → 시트 저장
                                   → No  → 참신함? → AI 요약/강점 정리 → 시트 저장
                                                                  → 아님  → 스킵
```


## How to Setup

### 1. Fork & `config.yaml` 수정

```yaml
pubmed:
  searches:
    - query: '("Radiology") AND ("AI")'
      sheet_name: "Radiology AI"
    - query: '("CT") AND ("deep learning")'
      sheet_name: "CT Deep Learning"

llm:
  model: "gpt-5-mini"  # gpt-5-nano | gpt-5-mini | gpt-5.2
```

### 2. GitHub Secrets 등록

**Settings → Secrets and variables → Actions**:

| Secret | 값 |
|--------|-----|
| `OPENAI_API_KEY` | OpenAI API 키 |
| `GOOGLE_SERVICE_ACCOUNT_JSON` | Google 서비스 계정 JSON 전체 |
| `SPREADSHEET_ID` | Google Sheets URL의 `/d/` 와 `/edit` 사이 문자열 (아래 참고) |
| `PUBMED_EMAIL` | PubMed API용 이메일  |

**SPREADSHEET_ID 찾는 법**

```
https://docs.google.com/spreadsheets/d/XXXXXXXX/edit
                                       ^^^^^^^^
                                       이 부분 
```

**Google Cloud 설정**

1. [Google Cloud Console](https://console.cloud.google.com)에서 프로젝트 생성
2. Service Account 생성 → JSON 키 다운로드
3. [Google Sheets API 활성화](https://console.developers.google.com/apis/api/sheets.googleapis.com/overview)
4. Google Sheets를 서비스 계정 이메일과 **편집자**로 공유

### 3. 실행

- **자동**: 3일마다 (GitHub Actions)
- **수동**: Actions 탭 → Run workflow

## Cost Estimation

| 모델 | Input ($/1M) | Output ($/1M) | 월간 예상* |
|------|-------------|--------------|-----------|
| **gpt-5-nano** | $0.05 | ~$0.40 | ~$0.02 |
| **gpt-5-mini** | $0.25 | $2.00 | ~$0.10 |
| **gpt-5.2** | $1.75 | $14.00 | ~$0.70 |

\* 3일마다 실행, 30개 논문/회 기준. PubMed API와 Google Sheets API는 무료.

## Miscellaneous

### Output 컬럼

Date | Title | Journal | Selection Criteria | Summary | Comment

### 저널 이름 매칭

`config.yaml`의 `high_if_journals`는 정확한 매칭 + 와일드카드를 지원합니다:

- `"Radiology"` → Radiology만 매칭 (European Radiology 제외)
- `"Nature*"` → Nature, Nature Medicine, Nature Biotechnology 등
- `"*Radiology"` → Radiology, European Radiology, Skeletal Radiology 등

저널 이름 확인 방법:

1. [PubMed](https://pubmed.ncbi.nlm.nih.gov/)에서 논문 검색
2. 논문 클릭 → 상세 페이지에서 **제목 바로 아래** 표시되는 저널명 확인
3. 예: `Radiology. 2025 Jan;314(1):e241536.` → 저널명은 `Radiology`
4. 이 저널명을 `config.yaml`의 `high_if_journals`에 그대로 사용

> PubMed에 표시되는 저널명과 정확히 일치해야 합니다 (대소문자는 무관).


### License

MIT — see [LICENSE](LICENSE)
