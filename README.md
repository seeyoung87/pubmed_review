# PubMed Review Automation

PubMed 최신 논문을 자동으로 검색하고, AI로 평가한 뒤 Google Sheets에 정리하는 자동화 도구입니다.

**매일 실행** → 새 논문 발견 → AI 평가 → 스프레드시트에 자동 저장

## ✨ Features

- 🔍 **자동 검색**: PubMed에서 설정한 쿼리로 최신 논문 자동 수집
- 🤖 **AI 평가**: OpenAI GPT로 논문 참신성 평가 및 요약 생성
- 📊 **자동 저장**: Google Sheets에 결과 자동 저장
- 🎯 **스마트 필터링**: High IF 저널 또는 참신한 논문만 선별
- ♻️ **중복 방지**: 이미 처리한 논문은 자동으로 스킵
- 💾 **배치 저장**: 10개씩 배치 저장으로 안정성 향상
- 🔄 **자동 재시도**: 네트워크 오류 시 exponential backoff으로 최대 4회 재시도
- 💰 **비용 최적화**: High IF 논문은 novelty 체크 생략

## 🚀 Quick Start

### 1. Fork & 설정

1. 이 저장소를 Fork
2. `config.yaml` 수정:
   ```yaml
   pubmed:
     email: "your_email@example.com"
     searches:
       - query: '("Radiology") AND ("AI")'
         sheet_name: "Radiology AI"

   sheets:
     spreadsheet_id: "YOUR_SHEET_ID"

   llm:
     model: "gpt-5-mini"  # or gpt-5-nano, gpt-5.2
   ```

### 2. GitHub Secrets 설정

**Settings → Secrets and variables → Actions**에서 추가:

- `OPENAI_API_KEY`: OpenAI API 키
- `GOOGLE_SERVICE_ACCOUNT_JSON`: Google 서비스 계정 JSON (전체 내용)
- `SPREADSHEET_ID`: Google Sheets ID

### 3. Google 설정

1. [Google Cloud Console](https://console.cloud.google.com)에서 프로젝트 생성
2. Service Account 생성 → JSON 키 다운로드
3. [Google Sheets API 활성화](https://console.developers.google.com/apis/api/sheets.googleapis.com/overview)
4. Google Sheets를 서비스 계정 이메일(`xxx@xxx.iam.gserviceaccount.com`)과 **편집자**로 공유

### 4. 실행

- **자동**: 3일마다 자동 실행 (GitHub Actions)
- **수동**: **Actions** 탭 → **Run workflow** 클릭

---

## 💰 Cost Estimation

### OpenAI GPT-5 시리즈 모델 비교

| 모델 | Input ($/1M) | Output ($/1M) | 월간 예상 비용* | 특징 |
|------|--------------|---------------|----------------|------|
| **gpt-5-nano** | $0.05 | ~$0.40** | ~$0.02 | 가장 빠르고 저렴, 분류/요약용 |
| **gpt-5-mini** | $0.25 | $2.00 | ~$0.10 | 균형잡힌 가격/성능 (권장) |
| **gpt-5.2** | $1.75 | $14.00 | ~$0.70 | 최고 성능, 코딩/복잡한 작업 |

\* 3일마다 실행, 논문 30개/회 기준
\** 정확한 output 가격 미공개, 일반적으로 input의 8-10배로 추정

### 비용 절감 팁

- **gpt-5-nano**: 단순 요약만 필요한 경우 최적
- **gpt-5-mini**: 대부분의 경우 권장 (현재 기본값)
- High IF 저널 리스트 확장 → novelty 체크 생략으로 API 호출 50% 감소
- `config.yaml`의 `retmax` 조정 (기본 200)

**무료 서비스:** Google Sheets API, PubMed API

---

## 📝 Configuration

### 다중 검색 설정 (여러 주제 → 각각 다른 시트)

```yaml
pubmed:
  email: "your_email@example.com"
  searches:
    - query: '("Radiology") AND ("large language model")'
      sheet_name: "Radiology NLP"
    - query: '("CT") AND ("deep learning")'
      sheet_name: "CT Deep Learning"
```

### LLM 프롬프트 커스터마이징

```yaml
llm:
  model: "gpt-5-mini"
  temperature: 0.2

  novelty_prompt: |
    이 논문이 혁신적인지 매우 엄격하게 평가하세요.
    Title: {title}
    Abstract: {abstract}

  summary_prompt: |
    논문을 2줄로 요약하고 강점을 1줄로 설명하세요.
    Title: {title}
    Abstract: {abstract}
```

### 주요 설정값

| 설정 | 기본값 | 설명 |
|------|--------|------|
| `llm.model` | gpt-5-mini | gpt-5-nano, gpt-5-mini, gpt-5.2 |
| `llm.temperature` | 0.2 | 창의성 수준 (0-2) |
| `filters.high_if_journals` | [Nature, Science, ...] | High IF 저널 리스트 |
| `pubmed.reldate` | 3 | 검색 기간(일) |
| `pubmed.retmax` | 200 | 최대 논문 수 |

---

## 📊 Output Format

| 컬럼 | 설명 |
|------|------|
| Date | 처리 날짜 |
| PMID | PubMed ID |
| Title | 논문 제목 |
| Journal | 저널명 |
| Publication Date | 발행일 |
| DOI | DOI |
| Selection Criteria | 선별 기준 (High IF / Novelty) |
| Novelty Reason | 참신성 근거 |
| Summary | 2줄 요약 |
| Strengths | 강점 1줄 |

---

## 🔧 Troubleshooting

### Google Sheets API is disabled
```
RuntimeError: Google Sheets API is disabled...
```
→ 에러 메시지의 링크를 클릭해서 API 활성화

### Permission denied
```
RuntimeError: Permission denied. Share the spreadsheet with: xxx@...
```
→ Google Sheets를 서비스 계정 이메일과 편집자로 공유

### 중복 논문이 계속 저장됨
→ 시트 이름이 `config.yaml`의 `sheet_name`과 정확히 일치하는지 확인

### LLM 비용이 너무 많이 나옴
1. `gpt-5-nano` 모델로 변경
2. `retmax` 줄이기 (예: 50)
3. High IF 저널 리스트 확장

---

## 📈 How It Works

```
PubMed 검색 → 메타데이터 수집 → 중복 체크
    ↓
각 논문마다:
  High IF 저널? → Yes → 요약 생성 → 저장
              → No → AI 참신성 평가 → 참신함? → Yes → 요약 생성 → 저장
                                            → No → 스킵
    ↓
10개씩 배치로 Google Sheets 저장
```

**비용 최적화**: High IF 논문은 novelty API 호출 생략 → 비용 50% 절감

---

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details

---

## 🔗 Sources

- [GPT-5 mini Model](https://platform.openai.com/docs/models/gpt-5-mini)
- [GPT-5 nano Model](https://platform.openai.com/docs/models/gpt-5-nano)
- [GPT-5.2 Model](https://platform.openai.com/docs/models/gpt-5.2)
- [OpenAI API Pricing](https://openai.com/api/pricing/)
