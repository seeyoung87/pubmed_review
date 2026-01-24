# PubMed Review Automation

PubMed 최신 논문을 자동으로 검색하고, AI로 평가한 뒤 Google Sheets에 정리하는 자동화 도구입니다.

**매일 실행** → 새 논문 발견 → AI 평가 → 스프레드시트에 자동 저장

## ✨ Features

- 🔍 **자동 검색**: PubMed에서 설정한 쿼리로 최신 논문 자동 수집
- 🤖 **AI 평가**: OpenAI GPT로 논문 참신성 평가 및 3줄 요약 생성
- 📊 **자동 저장**: Google Sheets에 결과 자동 저장 (컬럼 헤더 포함)
- 🎯 **스마트 필터링**: High IF 저널 또는 참신한 논문만 선별
- ♻️ **중복 방지**: 이미 처리한 논문은 자동으로 스킵
- 💾 **안전한 저장**: 10개씩 배치 저장으로 에러 시 데이터 손실 최소화
- 🔄 **자동 재시도**: 네트워크 오류 시 exponential backoff으로 최대 4회 재시도
- 💰 **비용 최적화**: High IF 논문은 novelty 체크 생략, 토큰 사용량 로깅

## 🚀 Quick Start

### 1. 저장소 Fork 및 설정

1. GitHub에서 이 저장소를 Fork
2. Fork한 저장소의 `config.yaml` 수정
   - `pubmed.email`: 본인 이메일
   - `pubmed.searches`: 검색 쿼리와 시트 이름
   - `sheets.spreadsheet_id`: Google Sheets ID

### 2. GitHub Secrets 설정

Fork한 저장소에서 **Settings → Secrets and variables → Actions** 이동 후 추가:

- `OPENAI_API_KEY`: OpenAI API 키
- `GOOGLE_SERVICE_ACCOUNT_JSON`: Google 서비스 계정 JSON (전체 내용)
- `SPREADSHEET_ID`: Google Sheets ID

<details>
<summary>📌 Google Service Account 생성</summary>

1. [Google Cloud Console](https://console.cloud.google.com) 접속
2. 프로젝트 생성 → APIs & Services → Credentials
3. Create Credentials → Service Account
4. Keys → Add Key → JSON 다운로드
5. JSON 파일 전체 내용을 `GOOGLE_SERVICE_ACCOUNT_JSON`에 복사
6. [Google Sheets API 활성화](https://console.developers.google.com/apis/api/sheets.googleapis.com/overview)

</details>

### 3. Google Sheets 준비

1. Google Sheets에서 새 스프레드시트 생성
2. URL에서 ID 복사: `https://docs.google.com/spreadsheets/d/`**`1AbC...xYz`**`/edit`
3. 서비스 계정 이메일(`xxx@xxx.iam.gserviceaccount.com`)을 **편집자**로 공유

### 4. 자동화 실행

**자동 실행:** 3일마다 자동으로 실행됩니다 (GitHub Actions)

**수동 실행:**
1. Fork한 저장소의 **Actions** 탭 이동
2. **PubMed Review Automation** 선택
3. **Run workflow** 클릭

**결과 확인:** Google Sheets에 자동으로 저장됩니다

<details>
<summary>📸 예상 결과물</summary>

| Date | PMID | Title | Journal | Selection | Summary |
|------|------|-------|---------|-----------|---------|
| 2026-01-24 | 38123456 | Deep learning for CT... | Radiology | High IF | This study uses... |
| 2026-01-24 | 38123457 | Novel AI for MRI | Other | Novelty | Introduces a... |

</details>

---

## 💰 Cost Estimation

**Monthly cost:** ~**$0.10 USD/month** (OpenAI API, gpt-4o-mini)
- 실행당: ~$0.01 USD
- 3일마다 자동 실행

**무료 서비스:**
- Google Sheets API: $0
- PubMed API: $0

---

## 📝 Configuration

<details>
<summary><b>config.yaml 상세 설정</b></summary>

### 단일 검색 설정

```yaml
pubmed:
  email: "your_email@example.com"
  search_query: '("Radiology") AND (("large language model") OR ("GPT"))'
  sheet_name: "Radiology NLP"  # 시트 탭 이름
  reldate: 3  # 최근 3일간 논문 (비우면 workflow.schedule_days 사용)
  retmax: 200  # 최대 논문 수

filters:
  high_if_journals: ["Nature", "Science", "Cell", "Nature Medicine", "The Lancet"]

llm:
  model: "gpt-4o-mini"
  temperature: 0.2

sheets:
  spreadsheet_id: "YOUR_SHEET_ID"

workflow:
  schedule_days: 3  # GitHub Actions 주기 (cron과 일치시켜야 함)
```

### 다중 검색 설정 (여러 주제를 각각 다른 시트에)

```yaml
pubmed:
  email: "your_email@example.com"
  searches:
    - query: '("Radiology") AND ("large language model")'
      sheet_name: "Radiology NLP"
    - query: '("CT") AND ("deep learning")'
      sheet_name: "CT Deep Learning"
    - query: '("MRI") AND ("artificial intelligence")'
      sheet_name: "MRI AI"
```

### LLM 프롬프트 커스터마이징

```yaml
llm:
  novelty_prompt: |
    다음 논문이 정말 새로운 방법론을 제시하는지 평가하세요.
    Title: {title}
    Journal: {journal}
    Abstract: {abstract}

  summary_prompt: |
    논문을 3줄 이내로 요약하고 강점을 설명하세요.
    Title: {title}
    Journal: {journal}
    Abstract: {abstract}
```

</details>

<details>
<summary><b>환경 변수 전체 목록</b></summary>

| 변수 | 필수 | 설명 | 기본값 |
|------|------|------|--------|
| `PUBMED_EMAIL` | ✅ | NCBI 연락용 이메일 | config.yaml |
| `OPENAI_API_KEY` | ✅ | OpenAI API 키 | - |
| `GOOGLE_SERVICE_ACCOUNT_JSON` | ✅ | Google 서비스 계정 JSON | - |
| `SPREADSHEET_ID` | ❌ | Google Sheets ID | config.yaml |
| `CONFIG_PATH` | ❌ | 설정 파일 경로 | `config.yaml` |
| `LOG_LEVEL` | ❌ | 로그 레벨 | `INFO` |
| `PUBMED_RELDATE` | ❌ | 검색 기간(일) | config.yaml |
| `PUBMED_RETMAX` | ❌ | 최대 논문 수 | `200` |

</details>

---

## 🔧 Troubleshooting

<details>
<summary><b>Google Sheets API is disabled</b></summary>

```
RuntimeError: Google Sheets API is disabled. Enable it at https://...
```

**해결**: 에러 메시지의 링크를 클릭해서 Google Sheets API 활성화

</details>

<details>
<summary><b>Permission denied</b></summary>

```
RuntimeError: Permission denied. Share the spreadsheet with: xxx@xxx.iam.gserviceaccount.com
```

**해결**: Google Sheets에서 서비스 계정 이메일을 편집자로 공유

</details>

<details>
<summary><b>Missing PUBMED_EMAIL</b></summary>

```
RuntimeError: Missing PUBMED_EMAIL or pubmed.email in config
```

**해결**: `config.yaml`에 `pubmed.email` 설정 또는 환경 변수 `PUBMED_EMAIL` 추가

</details>

<details>
<summary><b>중복 논문이 계속 저장됨</b></summary>

**원인**: 시트 이름이 잘못되었거나 PMID 컬럼(B열)이 없음

**해결**:
1. 시트 이름이 `config.yaml`의 `sheet_name`과 정확히 일치하는지 확인
2. 컬럼 헤더가 있는지 확인 (첫 실행 시 자동 생성됨)

</details>

<details>
<summary><b>LLM 비용이 너무 많이 나옴</b></summary>

**해결**:
1. `config.yaml`에서 `retmax` 줄이기 (예: 50)
2. `filters.high_if_journals`에 저널 추가 (novelty 체크 생략됨)
3. 로그에서 토큰 사용량 확인:
   ```
   LLM usage - prompt: 450, completion: 85, total: 535 tokens
   ```

</details>

---

## 📊 Output Format

Google Sheets에 다음과 같이 저장됩니다:

| 컬럼 | 설명 | 예시 |
|------|------|------|
| Date | 처리 날짜 | 2026-01-24 |
| PMID | PubMed ID | 38123456 |
| Title | 논문 제목 | Novel deep learning approach... |
| Journal | 저널명 | Radiology |
| Publication Date | 발행일 | 2026 Jan 15 |
| DOI | DOI | 10.1148/radiol.123456 |
| Selection Criteria | 선별 기준 | High IF, Novelty |
| Novelty Reason | 참신성 근거 | Introduces new architecture... |
| Summary | 3줄 요약 | This study presents a novel... |
| Strengths | 강점 | Strong validation on large dataset... |

---

## 📈 How It Works

```
PubMed 검색
    ↓
메타데이터 수집 (제목, 초록 등)
    ↓
중복 체크 (이미 처리한 PMID 스킵)
    ↓
각 논문마다:
  ├─ High IF 저널? → Yes → 요약 생성 → 저장
  └─ No → AI 참신성 평가 → 참신함? → Yes → 요약 생성 → 저장
                                    └─ No → 스킵
    ↓
10개씩 배치로 Google Sheets 저장
```

**비용 최적화**: High IF 논문은 novelty API 호출을 건너뛰어 비용 50% 절감

---

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details
