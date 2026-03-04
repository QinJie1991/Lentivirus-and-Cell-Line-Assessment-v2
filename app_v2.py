"""
病毒与细胞系评估系统 v2.5 - 最终整合版
Final Integrated Version with Gene Autocomplete
"""

import streamlit as st
import requests
import json
import time
import re
import html
import http.client
import uuid
import logging
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime
from functools import lru_cache
import difflib

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("virus_assessment")

# 页面配置
st.set_page_config(
    page_title="病毒与细胞系评估系统 v2.5",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== 安全配置 ====================

class SecurityConfig:
    """安全配置"""
    
    GENE_NAME_PATTERN = re.compile(r'^[a-zA-Z][a-zA-Z0-9]*(-?[a-zA-Z0-9]+)*$')
    MAX_GENE_LENGTH = 50
    
    ALLOWED_ORGANISMS = {
        'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'
    }
    
    @staticmethod
    def sanitize_input(text: str, max_length: int = 100) -> str:
        """清理用户输入"""
        if not text:
            return ""
        text = text.strip()[:max_length]
        text = ''.join(char for char in text if ord(char) >= 32 and char not in ['<', '>', '"', "'"])
        return text
    
    @staticmethod
    def validate_gene_name(gene_name: str) -> Tuple[bool, str]:
        """验证基因名格式"""
        if not gene_name:
            return False, "基因名不能为空"
        
        if len(gene_name) > SecurityConfig.MAX_GENE_LENGTH:
            return False, f"基因名过长（最大{SecurityConfig.MAX_GENE_LENGTH}字符）"
        
        if '..' in gene_name or '//' in gene_name or gene_name.startswith('-'):
            return False, "基因名包含非法字符组合"
        
        if not SecurityConfig.GENE_NAME_PATTERN.match(gene_name):
            return False, "基因名格式无效（必须以字母开头）"
        
        return True, ""

# ==================== 业务配置 ====================

class Config:
    MAX_LENTIVIRUS_CAPACITY = 8000
    WARNING_CAPACITY = 6000
    
    @staticmethod
    def get_ncbi_email():
        return st.secrets.get("NCBI_EMAIL", "user@example.com")
    
    @staticmethod
    def get_qwen_model():
        return st.secrets.get("AI_MODEL", "qwen-plus")

# ==================== 数据模型 ====================

@dataclass
class LiteratureEvidence:
    pmid: str
    title: str
    abstract: str
    evidence_text: str
    evidence_type: str
    confidence: float

@dataclass
class AIAnalysisResult:
    conclusion: str
    risk_level: str
    citations: List[LiteratureEvidence]
    raw_response: str
    validation_status: str

@dataclass
class HardRuleCheck:
    rule_name: str
    passed: bool
    reason: str
    source: str
    overrideable: bool = False

# ==================== API频率限制器（带Debounce） ====================

class APIRateLimiter:
    """API请求频率限制器"""
    def __init__(self, requests_per_second: float = 3.0):
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0
    
    def wait(self):
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()

ncbi_limiter = APIRateLimiter(3.0)
qwen_limiter = APIRateLimiter(5.0)

class DebounceTimer:
    """防抖计时器（用于基因输入）"""
    def __init__(self, delay: float = 0.3):
        self.delay = delay
        self.last_call = 0
    
    def should_trigger(self) -> bool:
        current = time.time()
        if current - self.last_call >= self.delay:
            self.last_call = current
            return True
        return False

# ==================== 基因自动完成服务（整合Clinical Tables API） ====================

class GeneAutocompleteService:
    """基因名自动完成服务 - 基于NCBI Clinical Tables API"""
    
    def __init__(self):
        self.clinical_tables_url = "https://clinicaltables.nlm.nih.gov/api/ncbi_genes/v3/search"
        self.cache = {}
        self.cache_timeout = 3600
    
    @st.cache_data(ttl=3600, show_spinner=False)
    def get_suggestions(_self, query: str, organism: str = "human", limit: int = 8) -> List[Dict]:
        """
        获取基因名建议（带防抖和缓存）
        使用NCBI Clinical Tables API（比E-utilities更快，支持模糊匹配）
        """
        if not query or len(query) < 2:
            return []
        
        try:
            # 物种映射
            organism_map = {
                "human": "Homo sapiens",
                "mouse": "Mus musculus",
                "rat": "Rattus norvegicus"
            }
            organism_name = organism_map.get(organism, organism)
            
            params = {
                "terms": query,
                "maxList": limit,
                "df": "symbol,name,chromosome,gene_id,type_of_gene",
                "q": f"organism:\"{organism_name}\""
            }
            
            response = requests.get(_self.clinical_tables_url, params=params, timeout=5)
            response.raise_for_status()
            
            data = response.json()
            
            if data and len(data) >= 3:
                results = []
                headers = data[0]
                rows = data[2]
                
                for row in rows:
                    gene_info = dict(zip(headers, row))
                    results.append({
                        "symbol": html.escape(gene_info.get("symbol", "")),
                        "name": html.escape(gene_info.get("name", "")),
                        "gene_id": gene_info.get("gene_id", ""),
                        "chromosome": html.escape(gene_info.get("chromosome", "")),
                        "type": html.escape(gene_info.get("type_of_gene", ""))
                    })
                return results
            
            return []
            
        except Exception as e:
            logger.warning(f"Gene suggestion error: {e}")
            return []

# ==================== 简化的基因输入组件 ====================

class GeneInputComponent:
    """简化的基因输入组件（避免过度复杂的状态管理）"""
    
    def __init__(self, gene_service: GeneAutocompleteService):
        self.gene_service = gene_service
        self.debounce = DebounceTimer(0.3)  # 300ms防抖
    
    def render(self, organism: str, key_prefix: str = "gene") -> Optional[str]:
        """
        渲染基因输入组件
        返回选中的基因symbol，或None
        """
        # 使用简洁的session state key
        input_key = f"{key_prefix}_input"
        selected_key = f"{key_prefix}_selected"
        
        # 初始化
        if input_key not in st.session_state:
            st.session_state[input_key] = ""
        if selected_key not in st.session_state:
            st.session_state[selected_key] = ""
        
        # 输入框
        col1, col2 = st.columns([1, 0.1])
        
        with col1:
            user_input = st.text_input(
                "基因名（支持自动完成）",
                placeholder="输入2个字符以上获取建议，如: TP, EG, KR...",
                value=st.session_state[input_key],
                key=f"{key_prefix}_text"
            )
        
        with col2:
            st.write("")  # 对齐
            st.write("")
            if st.button("🔄", key=f"{key_prefix}_refresh", help="清除选择"):
                st.session_state[input_key] = ""
                st.session_state[selected_key] = ""
                st.rerun()
        
        selected_gene = None
        
        # 获取建议（带防抖）
        if len(user_input) >= 2 and self.debounce.should_trigger():
            with st.spinner(""):
                suggestions = self.gene_service.get_suggestions(user_input, organism)
                st.session_state[f"{key_prefix}_suggestions"] = suggestions
        
        # 显示建议
        suggestions = st.session_state.get(f"{key_prefix}_suggestions", [])
        
        if suggestions and not st.session_state[selected_key]:
            st.markdown("**💡 点击选择基因：**")
            
            # 创建按钮网格
            cols = st.columns(min(len(suggestions), 4))
            for i, gene in enumerate(suggestions[:8]):  # 最多显示8个
                with cols[i % 4]:
                    btn_label = f"{gene['symbol']}\n<small>{gene['name'][:30]}...</small>"
                    if st.button(btn_label, key=f"{key_prefix}_btn_{i}", use_container_width=True):
                        st.session_state[selected_key] = gene['symbol']
                        st.session_state[input_key] = gene['symbol']
                        st.session_state[f"{key_prefix}_info"] = gene
                        st.rerun()
        
        # 显示已选基因信息卡片
        if st.session_state[selected_key]:
            gene_info = st.session_state.get(f"{key_prefix}_info", {})
            if gene_info:
                st.markdown(f"""
                <div style='padding: 12px; border-radius: 8px; background-color: #e3f2fd; 
                            border: 1px solid #2196F3; margin: 10px 0;'>
                    <h4 style='color: #1976D2; margin: 0;'>✅ 已选择: {gene_info.get('symbol', '')}</h4>
                    <p style='margin: 5px 0; font-size: 0.9em;'>
                        <b>全称:</b> {gene_info.get('name', '')}<br>
                        <b>染色体:</b> {gene_info.get('chromosome', '')} | 
                        <b>类型:</b> {gene_info.get('type', '')}<br>
                        <a href='https://www.ncbi.nlm.nih.gov/gene/{gene_info.get('gene_id', '')}' 
                           target='_blank' style='color: #1976D2;'>🔗 NCBI查看详情</a>
                    </p>
                </div>
                """, unsafe_allow_html=True)
            
            selected_gene = st.session_state[selected_key]
        
        return selected_gene

# ==================== NCBI E-utilities客户端（用于深度数据） ====================

class NCBIClient:
    """NCBI E-utilities客户端（用于获取转录本、文献等深度数据）"""
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.email = Config.get_ncbi_email()
    
    def _make_request(self, endpoint: str, params: Dict, retmode: str = "json") -> Optional[Dict]:
        ncbi_limiter.wait()
        
        params.update({
            'tool': 'LentivirusAssessment_v2',
            'email': self.email
        })
        if self.api_key:
            params['api_key'] = self.api_key
        
        url = f"{self.BASE_URL}/{endpoint}"
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json() if retmode == "json" else response.text
            
        except requests.exceptions.RequestException as e:
            logger.error(f"NCBI request failed: {type(e).__name__}")
            return None
        except Exception as e:
            logger.error(f"Unexpected NCBI error: {type(e).__name__}")
            return None
    
    @st.cache_data(ttl=3600, show_spinner=False)
    def fetch_gene_data(_self, gene_name: str, organism: str) -> Tuple[Dict, List[Dict]]:
        """获取基因详细数据（用于硬性规则检查）"""
        search_params = {
            'db': 'gene',
            'term': f"{gene_name}[Gene] AND {organism}[Organism]",
            'retmode': 'json',
            'retmax': 1
        }
        
        result = _self._make_request('esearch.fcgi', search_params)
        if not result:
            return {}, []
        
        gene_ids = result.get('esearchresult', {}).get('idlist', [])
        if not gene_ids:
            return {}, []
        
        gene_id = gene_ids[0]
        
        summary_params = {
            'db': 'gene',
            'id': gene_id,
            'retmode': 'json'
        }
        
        result = _self._make_request('esummary.fcgi', summary_params)
        if not result:
            return {}, []
        
        summary = result.get('result', {}).get(gene_id, {})
        
        gene_info = {
            'id': gene_id,
            'name': gene_name,
            'description': summary.get('description', ''),
            'organism': organism,
            'summary': summary.get('summary', '')
        }
        
        transcripts = _self._fetch_transcripts(gene_id)
        return gene_info, transcripts
    
    def _fetch_transcripts(self, gene_id: str) -> List[Dict]:
        """获取转录本长度（用于载体容量检查）"""
        try:
            search_params = {
                'db': 'nuccore',
                'term': f"{gene_id}[GeneID] AND (NM_[Title] OR XM_[Title])",
                'retmode': 'json',
                'retmax': 10
            }
            
            result = self._make_request('esearch.fcgi', search_params)
            if not result:
                return []
            
            ids = result.get('esearchresult', {}).get('idlist', [])
            if not ids:
                return []
            
            summary_params = {
                'db': 'nuccore',
                'id': ','.join(ids),
                'retmode': 'json'
            }
            
            result = self._make_request('esummary.fcgi', summary_params)
            if not result:
                return []
            
            docs = result.get('result', {})
            transcripts = []
            
            for uid in ids:
                try:
                    doc = docs.get(uid, {})
                    acc = doc.get('accessionversion', '')
                    length = doc.get('slen', 0)
                    
                    if (acc.startswith('NM_') or acc.startswith('XM_')) and length > 0:
                        transcripts.append({
                            'id': acc,
                            'length': int(length),
                            'title': str(doc.get('title', ''))[:200]
                        })
                except Exception:
                    continue
            
            return transcripts
            
        except Exception as e:
            logger.error(f"Transcript fetch error: {type(e).__name__}")
            return []
    
    def search_literature(self, gene_name: str, experiment_type: str) -> List[Dict]:
        """检索PubMed文献（用于AI分析）"""
        all_papers = []
        seen_pmids = set()
        
        keyword_map = {
            'knockout': ['knockout', 'CRISPR', 'deletion', 'deficiency'],
            'knockdown': ['knockdown', 'shRNA', 'siRNA', 'RNAi'],
            'overexpression': ['overexpression', 'transgenic', 'expression']
        }
        
        keywords = keyword_map.get(experiment_type.lower(), ['expression'])
        keywords.extend(['cell', 'viability'])
        
        for keyword in keywords:
            try:
                query = f"{gene_name} {keyword}"
                
                search_params = {
                    'db': 'pubmed',
                    'term': query,
                    'retmode': 'json',
                    'retmax': 3,
                    'sort': 'relevance'
                }
                
                result = self._make_request('esearch.fcgi', search_params)
                if not result:
                    continue
                
                pmids = result.get('esearchresult', {}).get('idlist', [])
                if not pmids:
                    continue
                
                new_pmids = [p for p in pmids if p not in seen_pmids]
                if not new_pmids:
                    continue
                
                fetch_params = {
                    'db': 'pubmed',
                    'id': ','.join(new_pmids),
                    'retmode': 'json'
                }
                
                result = self._make_request('esummary.fcgi', fetch_params)
                if not result:
                    continue
                
                docs = result.get('result', {})
                
                for pmid in new_pmids:
                    try:
                        doc = docs.get(pmid, {})
                        title = doc.get('title', '')
                        abstract = doc.get('abstract', '') or doc.get('sorttitle', '')
                        
                        if not title:
                            continue
                        
                        if not abstract:
                            abstract = f"[Title only] {title[:100]}"
                        
                        all_papers.append({
                            'pmid': str(pmid),
                            'title': html.escape(str(title)[:300]),
                            'abstract': html.escape(str(abstract)[:600]),
                            'keyword': html.escape(keyword),
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                        })
                        seen_pmids.add(pmid)
                        
                    except Exception:
                        continue
                
            except Exception:
                continue
        
        return all_papers

# ==================== 硬性规则引擎 ====================

class HardRulesEngine:
    """硬性规则引擎（物理限制和生物学事实）"""
    
    BASE_ESSENTIAL = {
        'ACTB', 'GAPDH', 'HSP90', 'HSPA5', 'RPL11', 'RPS3', 
        'PCNA', 'TOP2A', 'AURKB', 'PLK1', 'MYC'
    }
    
    BASE_TOXIC = {'BAX', 'BAK1', 'CASP3', 'TP53', 'MYC', 'FAS', 'TNF'}
    
    def check_all(self, gene_info: Dict, transcripts: List[Dict], 
                  experiment_type: str) -> Tuple[bool, List[HardRuleCheck]]:
        checks = []
        
        # 规则1: 载体容量
        check = self._check_vector_capacity(transcripts)
        checks.append(check)
        
        # 规则2: 必需基因
        if 'knockout' in experiment_type.lower():
            check = self._check_essential_gene(gene_info.get('name', ''))
            checks.append(check)
        
        # 规则3: 毒性基因
        if 'overexpression' in experiment_type.lower():
            check = self._check_known_toxic(gene_info.get('name', ''))
            checks.append(check)
        
        return all(c.passed for c in checks), checks
    
    def _check_vector_capacity(self, transcripts: List[Dict]) -> HardRuleCheck:
        """检查载体容量限制"""
        valid_lengths = [t.get('length', 0) for t in transcripts if t.get('length', 0) > 0]
        
        if not valid_lengths:
            return HardRuleCheck(
                rule_name="载体容量检查",
                passed=True,
                reason="转录本长度信息不可用，跳过检查",
                source="数据缺失",
                overrideable=True
            )
        
        max_length = max(valid_lengths)
        
        if max_length > Config.MAX_LENTIVIRUS_CAPACITY:
            return HardRuleCheck(
                rule_name="载体容量超限",
                passed=False,
                reason=f"最大转录本长度 {max_length:,}bp 超过物理极限 {Config.MAX_LENTIVIRUS_CAPACITY:,}bp",
                source="病毒学物理限制",
                overrideable=False
            )
        elif max_length > Config.WARNING_CAPACITY:
            return HardRuleCheck(
                rule_name="载体容量警告",
                passed=True,
                reason=f"转录本长度 {max_length:,}bp 接近极限",
                source="病毒学最佳实践",
                overrideable=True
            )
        else:
            return HardRuleCheck(
                rule_name="载体容量检查",
                passed=True,
                reason=f"转录本长度 {max_length:,}bp 在合理范围内",
                source="病毒学标准",
                overrideable=True
            )
    
    def _check_essential_gene(self, gene_name: str) -> HardRuleCheck:
        if gene_name.upper().strip() in self.BASE_ESSENTIAL:
            return HardRuleCheck(
                rule_name="必需基因检查",
                passed=False,
                reason=f"{html.escape(gene_name)} 是细胞必需基因，敲除可能导致细胞死亡",
                source="DepMap Database",
                overrideable=False
            )
        return HardRuleCheck(
            rule_name="必需基因检查",
            passed=True,
            reason=f"{html.escape(gene_name)} 未在必需基因核心列表中",
            source="DepMap Database",
            overrideable=True
        )
    
    def _check_known_toxic(self, gene_name: str) -> HardRuleCheck:
        if gene_name.upper().strip() in self.BASE_TOXIC:
            return HardRuleCheck(
                rule_name="毒性基因检查",
                passed=False,
                reason=f"{html.escape(gene_name)} 是已知促凋亡/毒性基因",
                source="文献数据库",
                overrideable=False
            )
        return HardRuleCheck(
            rule_name="毒性基因检查",
            passed=True,
            reason="非已知高毒性基因",
            source="文献数据库",
            overrideable=True
        )

# ==================== 文献验证算法（改进版） ====================

class CitationValidator:
    """文献引用验证器（使用序列相似度）"""
    
    @staticmethod
    def validate(evidence_text: str, source_abstract: str, threshold: float = 0.6) -> Tuple[bool, float]:
        """验证引用文本是否真实存在于原文中"""
        if not evidence_text or not source_abstract:
            return False, 0.0
        
        evidence_clean = evidence_text.lower().strip()
        abstract_clean = source_abstract.lower().strip()
        
        # 直接包含
        if evidence_clean in abstract_clean:
            return True, 1.0
        
        # 滑动窗口相似度
        best_ratio = 0.0
        window_size = min(len(evidence_clean), len(abstract_clean))
        
        if window_size < 10:
            return False, 0.0
        
        for i in range(0, len(abstract_clean) - window_size + 1, max(1, window_size // 4)):
            window = abstract_clean[i:i+window_size]
            ratio = difflib.SequenceMatcher(None, evidence_clean, window).quick_ratio()
            best_ratio = max(best_ratio, ratio)
            
            if best_ratio >= threshold:
                return True, best_ratio
        
        # 词覆盖检查
        evidence_words = set(evidence_clean.split())
        abstract_words = set(abstract_clean.split())
        
        if evidence_words:
            coverage = len(evidence_words & abstract_words) / len(evidence_words)
            if coverage >= 0.8:
                return True, max(best_ratio, coverage)
        
        return best_ratio >= threshold, best_ratio

# ==================== 阿里云千问分析器 ====================

class QwenRAGAnalyzer:
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.model = Config.get_qwen_model()
        self.validator = CitationValidator()
    
    def analyze(self, gene_name: str, experiment_type: str, literature: List[Dict]) -> AIAnalysisResult:
        """基于真实文献的AI分析"""
        if not self.api_key or not literature:
            return AIAnalysisResult(
                conclusion="未配置AI分析或未检索到文献",
                risk_level="unknown",
                citations=[],
                raw_response="",
                validation_status="skipped"
            )
        
        try:
            context = self._build_context(literature)
            prompt = self._build_prompt(gene_name, experiment_type, context)
            
            with st.spinner("🤖 AI分析文献中（约10-30秒）..."):
                response = self._call_api(prompt)
            
            if not response:
                return AIAnalysisResult(
                    conclusion="AI服务暂时不可用",
                    risk_level="error",
                    citations=[],
                    raw_response="",
                    validation_status="api_error"
                )
            
            parsed = self._parse_response(response, literature)
            validated = self._validate_citations(parsed, literature)
            
            return validated
            
        except Exception as e:
            error_id = str(uuid.uuid4())[:8]
            logger.error(f"AI analysis error [{error_id}]: {type(e).__name__}")
            return AIAnalysisResult(
                conclusion=f"AI分析异常，请联系管理员（错误ID: {error_id}）",
                risk_level="error",
                citations=[],
                raw_response="",
                validation_status="exception"
            )
    
    def _build_context(self, literature: List[Dict]) -> str:
        contexts = []
        for i, paper in enumerate(literature[:8], 1):
            contexts.append(f"""
【文献{i}】PMID: {paper['pmid']}
标题: {paper['title'][:100]}
摘要: {paper['abstract'][:400]}...
""")
        return "\n---\n".join(contexts)
    
    def _build_prompt(self, gene_name: str, exp_type: str, context: str) -> str:
        safe_gene = html.escape(gene_name)
        
        return f"""作为分子生物学专家，基于以下文献评估{safe_gene}在{exp_type}实验中的风险。

【规则】
1. 只能基于提供的文献判断
2. 每个结论必须标注【PMID: XXXXX】
3. 必须引用文献原文片段作为证据
4. 无证据时回答"文献未提及"
5. 风险等级：high/medium/low/none

【文献】
{context}

【输出JSON格式】
{{
    "overall_risk": "high/medium/low/none",
    "summary": "基于【PMID: xxx】...",
    "risks": [
        {{
            "type": "细胞毒性/增殖影响/其他",
            "level": "high/medium/low",
            "evidence": "文献原文片段",
            "pmid": "12345678"
        }}
    ],
    "suggestions": [
        {{"action": "建议策略", "basis": "基于【PMID: xxx】"}}
    ]
}}"""
    
    def _call_api(self, prompt: str) -> Optional[str]:
        try:
            import dashscope
            dashscope.api_key = self.api_key
            qwen_limiter.wait()
            
            response = dashscope.Generation.call(
                model=self.model,
                messages=[{'role': 'user', 'content': prompt}],
                result_format='message',
                temperature=0.1,
                max_tokens=2000
            )
            
            if response.status_code == 200:
                return response.output.choices[0].message.content
            return None
            
        except ImportError:
            return self._call_http(prompt)
        except Exception as e:
            logger.error(f"Qwen API error: {type(e).__name__}")
            return None
    
    def _call_http(self, prompt: str) -> Optional[str]:
        conn = None
        try:
            conn = http.client.HTTPSConnection("dashscope.aliyuncs.com", timeout=60)
            qwen_limiter.wait()
            
            payload = json.dumps({
                "model": self.model,
                "input": {"messages": [{"role": "user", "content": prompt}]},
                "parameters": {"temperature": 0.1, "max_tokens": 2000}
            })
            
            headers = {
                'Authorization': f'Bearer {self.api_key}',
                'Content-Type': 'application/json'
            }
            
            conn.request("POST", "/api/v1/services/aigc/text-generation/generation", 
                        payload, headers)
            res = conn.getresponse()
            data = json.loads(res.read().decode("utf-8"))
            
            if 'output' in data:
                return data['output'].get('text') or data['output'].get('choices', [{}])[0].get('message', {}).get('content')
            return None
            
        except Exception as e:
            logger.error(f"HTTP call error: {type(e).__name__}")
            return None
        finally:
            if conn:
                try:
                    conn.close()
                except:
                    pass
    
    def _parse_response(self, response: str, source_literature: List[Dict]) -> AIAnalysisResult:
        if not response:
            return AIAnalysisResult("AI返回空", "error", [], "", "empty")
        
        try:
            json_match = re.search(r'\{.*\}', response, re.DOTALL)
            if not json_match:
                return AIAnalysisResult("AI返回格式错误", "error", [], "", "no_json")
            
            data = json.loads(json_match.group())
            
            citations = []
            for risk in data.get('risks', []):
                pmid = str(risk.get('pmid', ''))
                evidence = str(risk.get('evidence', ''))[:200]
                
                source = next((p for p in source_literature if p['pmid'] == pmid), {})
                
                citations.append(LiteratureEvidence(
                    pmid=pmid,
                    title=source.get('title', 'Unknown')[:100],
                    abstract=source.get('abstract', '')[:200],
                    evidence_text=evidence,
                    evidence_type=str(risk.get('level', 'unknown')),
                    confidence=0.9 if risk.get('level') == 'high' else 0.7
                ))
            
            return AIAnalysisResult(
                conclusion=str(data.get('summary', ''))[:500],
                risk_level=str(data.get('overall_risk', 'unknown')),
                citations=citations,
                raw_response="",
                validation_status="parsed"
            )
            
        except Exception as e:
            logger.error(f"Parse error: {type(e).__name__}")
            return AIAnalysisResult("解析错误", "error", [], "", "parse_error")
    
    def _validate_citations(self, result: AIAnalysisResult, source_literature: List[Dict]) -> AIAnalysisResult:
        if not result.citations:
            return result
        
        available_pmids = {p['pmid'] for p in source_literature}
        valid_citations = []
        
        for cite in result.citations:
            if cite.pmid not in available_pmids:
                logger.warning(f"Hallucinated PMID: {cite.pmid}")
                continue
            
            source = next((p for p in source_literature if p['pmid'] == cite.pmid), {})
            is_valid, similarity = self.validator.validate(
                cite.evidence_text, 
                source.get('abstract', ''),
                threshold=0.6
            )
            
            if is_valid:
                cite.evidence_type = f"verified({similarity:.2f})"
                cite.confidence = similarity
                valid_citations.append(cite)
            elif similarity > 0.4:
                cite.evidence_type = f"paraphrased({similarity:.2f})"
                cite.confidence = similarity * 0.7
                valid_citations.append(cite)
        
        result.citations = valid_citations
        
        if not valid_citations and result.risk_level not in ['unknown', 'error']:
            result.validation_status = "unverified"
            result.risk_level = "uncertain"
            result.conclusion += " [警告：AI引用无法验证]"
        
        return result

# ==================== 主评估引擎 ====================

class HybridAssessmentEngine:
    def __init__(self, qwen_api_key: Optional[str] = None, ncbi_api_key: Optional[str] = None):
        self.hard_rules = HardRulesEngine()
        self.rag_analyzer = QwenRAGAnalyzer(qwen_api_key)
        self.ncbi = NCBIClient(ncbi_api_key)
    
    def assess(self, gene_name: str, organism: str, experiment_type: str) -> Dict:
        """执行完整评估"""
        result = {
            'timestamp': datetime.now().isoformat(),
            'gene': gene_name,
            'experiment': experiment_type,
            'decision_hierarchy': {},
            'final_recommendation': '',
            'primary_basis': ''
        }
        
        # 获取基因数据
        with st.spinner("🔍 检索基因详细数据..."):
            gene_info, transcripts = self.ncbi.fetch_gene_data(gene_name, organism)
        
        if not gene_info:
            return {'error': f'无法获取基因 {html.escape(gene_name)} 的信息'}
        
        result['gene_info'] = {
            'id': gene_info.get('id', ''),
            'name': gene_info.get('name', ''),
            'description': gene_info.get('description', '')[:200]
        }
        
        # 硬性规则检查
        with st.spinner("⚙️ 执行硬性规则检查..."):
            hard_passed, hard_checks = self.hard_rules.check_all(
                gene_info, transcripts, experiment_type
            )
        
        result['decision_hierarchy']['hard_rules'] = {
            'passed': hard_passed,
            'checks': [self._to_dict(c) for c in hard_checks]
        }
        
        # 检查阻断
        blocking = [c for c in hard_checks if not c.passed and not c.overrideable]
        if blocking:
            result['final_recommendation'] = 'BLOCKED'
            result['primary_basis'] = '硬性物理/生物学限制'
            result['reason'] = blocking[0].reason
            return result
        
        # AI文献分析
        with st.spinner("📚 检索文献并AI分析..."):
            literature = self.ncbi.search_literature(gene_name, experiment_type)
            
            ai_result = self.rag_analyzer.analyze(gene_name, experiment_type, literature) if literature else \
                       AIAnalysisResult("未检索到文献", "unknown", [], "", "no_literature")
        
        result['decision_hierarchy']['literature_based'] = {
            'literature_count': len(literature),
            'literatures': literature[:5],
            'ai_analysis': self._to_dict(ai_result)
        }
        
        # 生成建议
        if ai_result.risk_level == 'high':
            rec = "⚠️ 高风险：建议更换方案或采用诱导系统"
        elif ai_result.risk_level == 'medium':
            rec = "⚡ 中等风险：可进行，但需优化条件"
        elif ai_result.risk_level == 'low':
            rec = "✅ 低风险：可按标准流程进行"
        else:
            rec = "❓ 风险未知：缺乏文献支持，建议预实验"
        
        result['final_recommendation'] = rec
        result['primary_basis'] = f"基于{len(ai_result.citations)}篇文献验证" if ai_result.citations else "基于有限数据"
        
        return result
    
    def _to_dict(self, obj):
        if hasattr(obj, '__dataclass_fields__'):
            return {k: self._to_dict(v) for k, v in obj.__dict__.items()}
        elif isinstance(obj, list):
            return [self._to_dict(i) for i in obj]
        elif isinstance(obj, str):
            return obj[:500]
        return obj

# ==================== UI渲染 ====================

def render_header():
    st.markdown("""
    <h1 style='text-align: center; color: #1f77b4;'>
        🧬 病毒与细胞系评估系统 v2.5
    </h1>
    <p style='text-align: center; color: #666;'>
        智能基因自动完成 | 混合决策架构 | 安全加固
    </p>
    """, unsafe_allow_html=True)

def render_sidebar():
    with st.sidebar:
        st.header("⚙️ 配置")
        
        secret_qwen = st.secrets.get("DASHSCOPE_API_KEY") or st.secrets.get("QWEN_API_KEY")
        secret_ncbi = st.secrets.get("NCBI_API_KEY")
        
        st.info("API配置状态")
        
        user_qwen = st.text_input("阿里云千问 Key", type="password", 
                                  placeholder="留空使用Secrets")
        user_ncbi = st.text_input("NCBI Key", type="password", 
                                  placeholder="可选")
        
        final_qwen = user_qwen.strip() if user_qwen else secret_qwen
        final_ncbi = user_ncbi.strip() if user_ncbi else secret_ncbi
        
        st.divider()
        st.caption("""
        🔒 安全措施：
        - 输入验证与XSS防护
        - 300ms防抖机制
        - 强制文献引用验证
        - 异常信息脱敏
        """)
        
        return final_ncbi, final_qwen

def render_quick_genes(organism: str, gene_component: GeneInputComponent):
    """常用基因快捷按钮"""
    st.markdown("**⚡ 常用基因快速选择：**")
    
    common_genes = {
        "human": ["TP53", "EGFR", "KRAS", "PTEN", "BRAF", "MYC", "BCL2", "CASP3", "BAX"],
        "mouse": ["Trp53", "Egfr", "Kras", "Pten", "Braf", "Myc", "Bcl2", "Casp3", "Bax"],
        "rat": ["Tp53", "Egfr", "Kras", "Pten", "Braf", "Myc", "Bcl2", "Casp3", "Bax"]
    }
    
    genes = common_genes.get(organism, [])
    cols = st.columns(5)
    
    for i, gene in enumerate(genes):
        with cols[i % 5]:
            if st.button(gene, key=f"quick_{gene}", use_container_width=True):
                st.session_state["main_gene_input_selected"] = gene
                st.session_state["main_gene_input_input"] = gene
                # 触发自动完成的建议获取
                st.rerun()

def main():
    """主函数"""
    render_header()
    ncbi_key, qwen_key = render_sidebar()
    
    # 初始化服务
    gene_service = GeneAutocompleteService()
    gene_component = GeneInputComponent(gene_service)
    
    # 物种选择
    organism = st.selectbox(
        "物种",
        ["human", "mouse", "rat"],
        format_func=lambda x: {"human": "人类", "mouse": "小鼠", "rat": "大鼠"}.get(x, x)
    )
    
    # 常用基因快捷按钮（放在输入框上方）
    render_quick_genes(organism, gene_component)
    
    # 基因输入（带自动完成）
    st.markdown("---")
    gene = gene_component.render(organism, key_prefix="main_gene_input")
    
    # 实验类型
    exp_type = st.selectbox(
        "实验类型",
        ["knockout", "knockdown", "overexpression"],
        format_func=lambda x: {
            "knockout": "敲除 (CRISPR)",
            "knockdown": "敲低 (shRNA)",
            "overexpression": "过表达"
        }.get(x, x)
    )
    
    analyze = st.button("🔍 开始智能评估", type="primary", use_container_width=True)
    
    if analyze:
        if not gene:
            st.error("请选择或输入一个基因")
            return
        
        if not qwen_key:
            st.error("请配置阿里云千问API Key")
            return
        
        # 验证
        is_valid, error_msg = SecurityConfig.validate_gene_name(gene)
        if not is_valid:
            st.error(f"输入验证失败: {error_msg}")
            return
        
        # 映射物种
        organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus", "rat": "Rattus norvegicus"}
        organism_clean = organism_map.get(organism, organism)
        
        gene_clean = SecurityConfig.sanitize_input(gene, 50)
        exp_type_clean = SecurityConfig.sanitize_input(exp_type, 50)
        
        # 执行评估
        try:
            engine = HybridAssessmentEngine(qwen_api_key=qwen_key, ncbi_api_key=ncbi_key)
            
            with st.spinner("正在进行混合决策评估..."):
                result = engine.assess(gene_clean, organism_clean, exp_type_clean)
            
            if 'error' in result:
                st.error(html.escape(result['error']))
                return
            
            # 显示结果
            st.divider()
            st.markdown(f"## 🎯 评估报告 - {html.escape(result['gene'])}")
            
            rec = result['final_recommendation']
            rec_color = {"❌": "#ffebee", "⚠️": "#fff3e0", "⚡": "#fff8e1", "✅": "#e8f5e9"}.get(rec[:2], "#f5f5f5")
            
            st.markdown(f"""
            <div style='padding: 20px; background-color: {rec_color}; 
                        border-radius: 10px; text-align: center; margin: 20px 0;'>
                <h3>{html.escape(rec)}</h3>
                <small>{html.escape(result.get('primary_basis', ''))}</small>
            </div>
            """, unsafe_allow_html=True)
            
            # 详细结果
            hierarchy = result['decision_hierarchy']
            
            # 硬性规则
            st.markdown("### 🚦 硬性规则检查")
            for check in hierarchy['hard_rules']['checks']:
                icon = "✅" if check['passed'] else "❌" if not check['overrideable'] else "⚠️"
                color = "green" if check['passed'] else "red" if not check['overrideable'] else "orange"
                st.markdown(f"""
                <div style='padding: 10px; border-left: 4px solid {color}; 
                            background-color: #f8f9fa; margin: 5px 0;'>
                    <b>{icon} {html.escape(check['rule_name'])}</b><br/>
                    {html.escape(check['reason'])}<br/>
                    <small>来源: {html.escape(check['source'])}</small>
                </div>
                """, unsafe_allow_html=True)
            
            # AI分析
            lit_data = hierarchy.get('literature_based', {})
            if lit_data:
                st.markdown("### 📚 AI文献分析")
                st.write(f"检索到 {lit_data.get('literature_count', 0)} 篇文献")
                
                ai_analysis = lit_data.get('ai_analysis', {})
                if ai_analysis:
                    risk = ai_analysis.get('risk_level', 'unknown')
                    st.markdown(f"**风险等级**: {risk.upper()}")
                    st.markdown(f"**结论**: {html.escape(ai_analysis.get('conclusion', ''))}")
                    
                    citations = ai_analysis.get('citations', [])
                    if citations:
                        st.markdown("**验证的引用：**")
                        for cite in citations:
                            st.markdown(f"""
                            - **PMID**: {html.escape(cite.get('pmid', ''))} 
                              [{html.escape(cite.get('evidence_type', ''))}]
                            - 证据: "{html.escape(cite.get('evidence_text', '')[:100])}..."
                            """)
            
            # 审计信息
            with st.expander("🔍 审计信息"):
                st.json(result)
                
        except Exception as e:
            error_id = str(uuid.uuid4())[:8]
            logger.exception(f"Unhandled error: {e}")
            st.error(f"系统错误（ID: {error_id}），请联系管理员")

if __name__ == "__main__":
    main()
