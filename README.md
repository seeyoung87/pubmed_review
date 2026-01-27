# PubMed Review Automation

PubMed 최신 논문을 자동 검색 → AI 평가 → Google Sheets 저장하는 도구.

## What It Does

```
PubMed 검색 → 중복 체크 → Selected Journal? → Yes → 요약 생성 → 시트 저장
                                             → No  → AI 참신성 평가 → 참신함? → 요약 생성 → 시트 저장
                                                                    → 아님  → 스킵
```

- 3일마다 GitHub Actions로 자동 실행
- Selected Journal 논문은 novelty 체크 생략 (비용 절감)
- 10개씩 배치 저장, 네트워크 오류 시 자동 재시도

## Setup

### 1. Fork & `config.yaml` 수정

```yaml
pubmed:
  email: "your_email@example.com"
  searches:
    - query: '("Radiology") AND ("AI")'
      sheet_name: "Radiology AI"
    - query: '("CT") AND ("deep learning")'
      sheet_name: "CT Deep Learning"

sheets:
  spreadsheet_id: "YOUR_SHEET_ID"

llm:
  model: "gpt-5-mini"  # gpt-5-nano | gpt-5-mini | gpt-5.2
```

### 2. GitHub Secrets 등록

**Settings → Secrets and variables → Actions**:

| Secret | 값 |
|--------|-----|
| `OPENAI_API_KEY` | OpenAI API 키 |
| `GOOGLE_SERVICE_ACCOUNT_JSON` | Google 서비스 계정 JSON 전체 |
| `SPREADSHEET_ID` | Google Sheets ID |

### 3. Google Cloud 설정

1. [Google Cloud Console](https://console.cloud.google.com)에서 프로젝트 생성
2. Service Account 생성 → JSON 키 다운로드
3. [Google Sheets API 활성화](https://console.developers.google.com/apis/api/sheets.googleapis.com/overview)
4. Google Sheets를 서비스 계정 이메일과 **편집자**로 공유

### 4. 실행

- **자동**: 3일마다 (GitHub Actions)
- **수동**: Actions 탭 → Run workflow

## Cost Estimation

| 모델 | Input ($/1M) | Output ($/1M) | 월간 예상* |
|------|-------------|--------------|-----------|
| **gpt-5-nano** | $0.05 | ~$0.40 | ~$0.02 |
| **gpt-5-mini** | $0.25 | $2.00 | ~$0.10 |
| **gpt-5.2** | $1.75 | $14.00 | ~$0.70 |

\* 3일마다 실행, 30개 논문/회 기준. PubMed API와 Google Sheets API는 무료.

**비용 줄이기**: 모델을 `gpt-5-nano`로 변경 / `retmax` 줄이기 / Selected Journal 리스트 확장 (novelty 호출 생략)

## Misc

### Output 컬럼

Date | Title | Journal | Selection Criteria | Summary | Comment

### 저널 이름 매칭

`config.yaml`의 `high_if_journals`는 정확한 매칭 + 와일드카드를 지원합니다:

- `"Radiology"` → Radiology만 매칭 (European Radiology 제외)
- `"Nature*"` → Nature, Nature Medicine, Nature Biotechnology 등
- `"*Radiology"` → Radiology, European Radiology, Skeletal Radiology 등

저널 이름은 [PubMed](https://pubmed.ncbi.nlm.nih.gov/) 논문 상세 페이지에서 확인하세요.

### Troubleshooting

- **Google Sheets API is disabled** → 에러 메시지의 링크 클릭해서 API 활성화
- **Permission denied** → Sheets를 서비스 계정 이메일과 편집자로 공유
- **비용 과다** → `gpt-5-nano` 모델 사용 / `retmax` 줄이기

### License

MIT — see [LICENSE](LICENSE)
