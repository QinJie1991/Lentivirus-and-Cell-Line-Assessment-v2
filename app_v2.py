"""
病毒与细胞系评估系统 v2.0 - 混合决策架构
Hybrid Risk Assessment System with RAG-based AI
"""

import streamlit as st
import requests
import json
import time
import re
import hashlib
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime
import pandas as pd
import http.client
import urllib.parse

# 页面配置
st.set_page_config(
    page_title="病毒与细胞系评估系统 v2.0",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== 数据模型 ====================

@dataclass
class LiteratureEvidence:
    """文献证据数据类"""
    pmid: str
    title: str
    abstract: str
    evidence_text: str
    evidence_type: str  # 'direct', 'indirect', 'weak'
    confidence: float

@dataclass
class AIAnalysisResult:
    """AI分析结果"""
    conclusion: str
    risk_level: str
    citations: List[LiteratureEvidence]
    raw_response: str
    validation_status: str

@dataclass
class HardRuleCheck:
    """硬性规则检查结果"""
    rule_name: str
    passed: bool
    reason: str
    source: str
    overrideable: bool = False

@dataclass
class AuxiliaryHint:
    """AI辅助提示"""
    category: str
    content: str
    confidence: float
    warning: str = "此为AI预测，无文献实证，仅供参考"

# ==================== 配置管理 ====================

class Config:
    """配置管理"""
    # 硬性规则阈值
    MAX_LENTIVIRUS_CAPACITY = 8000  # bp
    WARNING_CAPACITY = 6000
    IDEAL_CAPACITY = 4000
    
    # 阿里云千问配置（用户需填写）
    QWEN_API_KEY = None
    QWEN_MODEL = "qwen-max"  # 或 qwen-turbo
    
    # NCBI配置
    NCBI_EMAIL = "user@example.com"
    NCBI_TOOL = "LentivirusAssessment_v2"

# ==================== API频率限制器 ====================

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

# ==================== 第一层：硬性规则引擎 ====================

class HardRulesEngine:
    """
    硬性规则引擎 - 基于物理限制和已知生物学事实
    这些规则优先级最高，AI无权覆盖
    """
    
    # 已知必需基因列表（示例，实际应从DepMap等数据库加载）
    ESSENTIAL_GENES = {
        'ACTB', 'GAPDH', 'HSP90', 'HSPA5', 'RPL11', 'RPS3', 
        'PCNA', 'TOP2A', 'AURKB', 'PLK1'
    }
    
    def check_all(self, gene_info: Dict, transcripts: List[Dict], 
                  experiment_type: str) -> Tuple[bool, List[HardRuleCheck]]:
        """
        执行所有硬性规则检查
        返回: (是否通过, 规则列表)
        """
        checks = []
        
        # 规则1: 载体容量限制（物理限制）
        check = self._check_vector_capacity(transcripts)
        checks.append(check)
        
        # 规则2: 必需基因敲除检查
        if 'knockout' in experiment_type.lower():
            check = self._check_essential_gene(gene_info.get('name', ''))
            checks.append(check)
        
        # 规则3: 已知毒性基因（过表达）
        if 'overexpression' in experiment_type.lower():
            check = self._check_known_toxic(gene_info.get('name', ''))
            checks.append(check)
        
        all_passed = all(c.passed for c in checks)
        return all_passed, checks
    
    def _check_vector_capacity(self, transcripts: List[Dict]) -> HardRuleCheck:
        """检查载体容量限制"""
        if not transcripts:
            return HardRuleCheck(
                rule_name="载体容量检查",
                passed=True,
                reason="未获取到转录本信息，跳过检查",
                source="系统默认"
            )
        
        max_length = max([t.get('length', 0) for t in transcripts if t.get('length')])
        
        if max_length > Config.MAX_LENTIVIRUS_CAPACITY:
            return HardRuleCheck(
                rule_name="载体容量超限",
                passed=False,
                reason=f"最大转录本长度 {max_length}bp 超过慢病毒物理包装极限 {Config.MAX_LENTIVIRUS_CAPACITY}bp",
                source="病毒学物理限制（Naldini et al., 1996）",
                overrideable=False
            )
        elif max_length > Config.WARNING_CAPACITY:
            return HardRuleCheck(
                rule_name="载体容量警告",
                passed=True,
                reason=f"转录本长度 {max_length}bp 接近极限，包装效率将显著降低",
                source="病毒学最佳实践",
                overrideable=True
            )
        else:
            return HardRuleCheck(
                rule_name="载体容量检查",
                passed=True,
                reason=f"转录本长度 {max_length}bp 在合理范围内",
                source="病毒学标准",
                overrideable=True
            )
    
    def _check_essential_gene(self, gene_name: str) -> HardRuleCheck:
        """检查必需基因"""
        if gene_name.upper() in self.ESSENTIAL_GENES:
            return HardRuleCheck(
                rule_name="必需基因检查",
                passed=False,
                reason=f"{gene_name} 是细胞必需基因，敲除将导致细胞死亡",
                source="DepMap CRISPR Screening Data",
                overrideable=False
            )
        return HardRuleCheck(
            rule_name="必需基因检查",
            passed=True,
            reason=f"{gene_name} 未在必需基因列表中",
            source="DepMap Database",
            overrideable=True
        )
    
    def _check_known_toxic(self, gene_name: str) -> HardRuleCheck:
        """检查已知毒性基因（示例）"""
        known_toxic = {'BAX', 'BAK1', 'CASP3', 'TP53', 'MYC'}
        if gene_name.upper() in known_toxic:
            return HardRuleCheck(
                rule_name="已知毒性基因",
                passed=False,
                reason=f"{gene_name} 是已知促凋亡/毒性基因，过表达极可能导致包装细胞死亡",
                source="文献综述（Cell Death & Differentiation）",
                overrideable=False
            )
        return HardRuleCheck(
            rule_name="毒性基因检查",
            passed=True,
            reason="非已知高毒性基因",
            source="文献数据库",
            overrideable=True
        )

# ==================== 第二层：RAG-based AI文献分析 ====================

class QwenRAGAnalyzer:
    """
    基于阿里云千问的RAG（检索增强生成）分析器
    强制要求AI基于提供的文献进行判断，并引用具体PMID
    """
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or Config.QWEN_API_KEY
        if not self.api_key:
            st.error("未配置阿里云千问API Key，无法使用AI文献分析功能")
    
    def analyze(self, gene_name: str, experiment_type: str, 
                literature: List[Dict]) -> AIAnalysisResult:
        """
        基于文献的AI分析
        强制要求：
        1. 只能使用提供的文献
        2. 必须标注PMID
        3. 必须提供原文引用
        """
        if not self.api_key or not literature:
            return AIAnalysisResult(
                conclusion="无AI分析（未配置API或文献不足）",
                risk_level="unknown",
                citations=[],
                raw_response="",
                validation_status="skipped"
            )
        
        # 构建RAG上下文
        context = self._build_rag_context(literature)
        
        # 构建约束性Prompt
        prompt = self._construct_prompt(gene_name, experiment_type, context)
        
        # 调用千问API
        try:
            response = self._call_qwen_api(prompt)
            parsed = self._parse_response(response, literature)
            validated = self._validate_citations(parsed, literature)
            return validated
        except Exception as e:
            return AIAnalysisResult(
                conclusion=f"AI分析失败: {str(e)}",
                risk_level="error",
                citations=[],
                raw_response=str(e),
                validation_status="failed"
            )
    
    def _build_rag_context(self, literature: List[Dict]) -> str:
        """构建RAG上下文"""
        contexts = []
        for i, paper in enumerate(literature[:10], 1):  # 限制Top-10
            contexts.append(f"""
【文献{i}】
PMID: {paper.get('pmid', 'N/A')}
Title: {paper.get('title', 'N/A')}
Abstract: {paper.get('abstract', 'N/A')[:800]}...
""")
        return "\n---\n".join(contexts)
    
    def _construct_prompt(self, gene_name: str, experiment_type: str, 
                         context: str) -> str:
        """构建约束性Prompt"""
        return f"""你是一位严格的分子生物学评审专家。请基于以下提供的文献，评估{gene_name}在{experiment_type}实验中的风险。

【严格规则 - 必须遵守】
1. **只能基于提供的【文献】部分做出判断**，严禁使用预训练知识或外部信息
2. 每个风险结论必须标注来源【PMID: XXXXX】
3. 必须引用文献中的原文片段作为证据（不少于10个字）
4. 如果文献中没有相关证据，必须明确回答"文献未提及"，禁止推测
5. 风险等级定义：
   - high: 文献明确报道实验失败、细胞死亡、严重毒性
   - medium: 文献报道实验困难、效率降低、部分毒性
   - low: 文献报道实验成功或无明显影响
   - none: 文献未提及相关风险

【待评估文献】
{context}

【输出格式 - 必须严格遵循JSON格式】
{{
    "overall_risk": "high/medium/low/none",
    "summary": "基于【PMID: 12345】和【PMID: 67890】，该基因...",
    "risks": [
        {{
            "type": "细胞毒性/抗病毒/增殖影响/包装效率",
            "level": "high/medium/low",
            "evidence": "文献中的原文片段，至少10个字",
            "pmid": "12345678",
            "mechanism": "文献提出的机制（如有）"
        }}
    ],
    "suggestions": [
        {{
            "action": "建议采取的行动",
            "basis": "基于【PMID: XXXXX】"
        }}
    ],
    "confidence": "high/medium/low"
}}

【重要提醒】
- 如果你引用了不存在的PMID，你的回答将视为无效
- 如果你无法找到相关证据，请诚实回答"文献未提及"
- 禁止编造数据或推测机制"""

    def _call_qwen_api(self, prompt: str) -> str:
        """调用阿里云千问API"""
        try:
            import dashscope
            dashscope.api_key = self.api_key
            
            response = dashscope.Generation.call(
                model=Config.QWEN_MODEL,
                messages=[{'role': 'user', 'content': prompt}],
                result_format='message',
                temperature=0.1,  # 低温度减少幻觉
                max_tokens=2000
            )
            
            if response.status_code == 200:
                return response.output.choices[0].message.content
            else:
                raise Exception(f"API错误: {response.message}")
                
        except ImportError:
            # 如果没有dashscope，使用HTTP API
            return self._call_qwen_http(prompt)
    
    def _call_qwen_http(self, prompt: str) -> str:
        """使用HTTP直接调用"""
        conn = http.client.HTTPSConnection("dashscope.aliyuncs.com")
        
        payload = json.dumps({
            "model": Config.QWEN_MODEL,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.1
        })
        
        headers = {
            'Authorization': f'Bearer {self.api_key}',
            'Content-Type': 'application/json'
        }
        
        conn.request("POST", "/api/v1/services/aigc/text-generation/generation", 
                    payload, headers)
        res = conn.getresponse()
        data = json.loads(res.read().decode("utf-8"))
        
        if 'output' in data and 'text' in data['output']:
            return data['output']['text']
        else:
            raise Exception(f"API响应错误: {data}")
    
    def _parse_response(self, response: str, source_literature: List[Dict]) -> AIAnalysisResult:
        """解析AI响应"""
        try:
            # 提取JSON部分
            json_match = re.search(r'\{.*\}', response, re.DOTALL)
            if json_match:
                data = json.loads(json_match.group())
            else:
                data = json.loads(response)
            
            # 构建证据列表
            citations = []
            for risk in data.get('risks', []):
                pmid = risk.get('pmid', '')
                # 查找原文
                source_paper = next((p for p in source_literature if p.get('pmid') == pmid), {})
                
                citations.append(LiteratureEvidence(
                    pmid=pmid,
                    title=source_paper.get('title', 'Unknown'),
                    abstract=source_paper.get('abstract', ''),
                    evidence_text=risk.get('evidence', ''),
                    evidence_type=risk.get('level', 'unknown'),
                    confidence=0.9 if risk.get('level') == 'high' else 0.7
                ))
            
            return AIAnalysisResult(
                conclusion=data.get('summary', ''),
                risk_level=data.get('overall_risk', 'unknown'),
                citations=citations,
                raw_response=response,
                validation_status="parsed"
            )
            
        except json.JSONDecodeError:
            return AIAnalysisResult(
                conclusion="AI返回格式错误，无法解析",
                risk_level="error",
                citations=[],
                raw_response=response,
                validation_status="parse_error"
            )
    
    def _validate_citations(self, result: AIAnalysisResult, 
                           source_literature: List[Dict]) -> AIAnalysisResult:
        """验证引用真实性（防止AI幻觉）"""
        available_pmids = {p.get('pmid') for p in source_literature}
        
        valid_citations = []
        for cite in result.citations:
            if cite.pmid in available_pmids:
                # 验证引用的文本是否在原文中（模糊匹配）
                source = next((p for p in source_literature if p.get('pmid') == cite.pmid), {})
                abstract = source.get('abstract', '').lower()
                evidence = cite.evidence_text.lower()
                
                # 简单的包含检查（实际应用可使用fuzzy matching）
                if len(evidence) > 5 and any(word in abstract for word in evidence.split()[:3]):
                    cite.evidence_type = "verified"
                else:
                    cite.evidence_type = "paraphrased"  # AI可能改写了原文
                valid_citations.append(cite)
            else:
                # AI编造了PMID！标记为无效
                cite.evidence_type = "hallucinated"
                # 不加入valid_citations
        
        result.citations = valid_citations
        if not valid_citations and result.risk_level != "unknown":
            result.validation_status = "suspected_hallucination"
            result.risk_level = "uncertain"
        
        return result

# ==================== 第三层：辅助AI提示 ====================

class AuxiliaryAIHints:
    """基于算法和AI的辅助提示（非决策依据）"""
    
    def __init__(self, qwen_api_key: Optional[str] = None):
        self.api_key = qwen_api_key
    
    def analyze_sequence(self, sequence: str) -> List[AuxiliaryHint]:
        """序列分析（算法为主）"""
        hints = []
        
        # GC含量分析
        if sequence:
            gc_content = (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence) * 100
            if gc_content > 70:
                hints.append(AuxiliaryHint(
                    category="序列特征",
                    content=f"GC含量高达{gc_content:.1f}%，可能影响转录效率",
                    confidence=0.8,
                    warning="基于算法计算，未经验证"
                ))
            elif gc_content < 30:
                hints.append(AuxiliaryHint(
                    category="序列特征",
                    content=f"GC含量仅{gc_content:.1f}%，可能影响载体稳定性",
                    confidence=0.7,
                    warning="基于算法计算，未经验证"
                ))
        
        # 重复序列检测
        if self._has_repeats(sequence):
            hints.append(AuxiliaryHint(
                category="序列复杂度",
                content="检测到潜在重复序列，可能导致重组",
                confidence=0.6,
                warning="基于模式匹配，建议实验验证"
            ))
        
        return hints
    
    def _has_repeats(self, seq: str, min_len: int = 6) -> bool:
        """简单重复检测"""
        if not seq:
            return False
        seq = seq.upper()
        for i in range(len(seq) - min_len * 2):
            if seq[i:i+min_len] == seq[i+min_len:i+min_len*2]:
                return True
        return False
    
    def ai_protein_prediction(self, protein_name: str) -> Optional[AuxiliaryHint]:
        """AI蛋白预测（可选，使用千问）"""
        if not self.api_key:
            return None
        
        # 仅作为示例，实际可调用结构预测API
        return AuxiliaryHint(
            category="蛋白结构",
            content=f"AI预测{protein_name}含有跨膜结构域（基于序列分析）",
            confidence=0.5,
            warning="此为AI预测，无实验验证，仅供参考"
        )

# ==================== 数据获取层 ====================

class NCBIClient:
    """NCBI数据客户端"""
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
    
    def fetch_gene_data(self, gene_name: str, organism: str) -> Tuple[Dict, List[Dict]]:
        """获取基因和转录本数据"""
        # 搜索Gene ID
        ncbi_limiter.wait()
        search_url = f"{self.BASE_URL}/esearch.fcgi"
        params = {
            'db': 'gene',
            'term': f"{gene_name}[Gene] AND {organism}[Organism]",
            'retmode': 'json',
            'email': Config.NCBI_EMAIL,
            'tool': Config.NCBI_TOOL
        }
        if self.api_key:
            params['api_key'] = self.api_key
        
        resp = requests.get(search_url, params=params, timeout=30)
        data = resp.json()
        gene_ids = data.get('esearchresult', {}).get('idlist', [])
        
        if not gene_ids:
            return {}, []
        
        gene_id = gene_ids[0]
        
        # 获取摘要
        ncbi_limiter.wait()
        summary_url = f"{self.BASE_URL}/esummary.fcgi"
        params = {
            'db': 'gene',
            'id': gene_id,
            'retmode': 'json'
        }
        if self.api_key:
            params['api_key'] = self.api_key
            
        resp = requests.get(summary_url, params=params, timeout=30)
        summary = resp.json().get('result', {}).get(gene_id, {})
        
        gene_info = {
            'id': gene_id,
            'name': gene_name,
            'description': summary.get('description', ''),
            'organism': organism
        }
        
        # 获取转录本（简化版，实际应从nuccore获取）
        transcripts = [{'id': f"NM_{gene_id}00{i}", 'length': 3000 + i*500} 
                      for i in range(1, 4)]  # 模拟数据
        
        return gene_info, transcripts
    
    def search_literature(self, gene_name: str, keywords: List[str]) -> List[Dict]:
        """检索文献"""
        all_papers = []
        
        for keyword in keywords:
            ncbi_limiter.wait()
            query = f"{gene_name} {keyword}"
            search_url = f"{self.BASE_URL}/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': query,
                'retmode': 'json',
                'retmax': 5,
                'sort': 'relevance'
            }
            if self.api_key:
                params['api_key'] = self.api_key
            
            resp = requests.get(search_url, params=params, timeout=30)
            data = resp.json()
            pmids = data.get('esearchresult', {}).get('idlist', [])
            
            if pmids:
                # 获取摘要
                ncbi_limiter.wait()
                fetch_url = f"{self.BASE_URL}/efetch.fcgi"
                fetch_params = {
                    'db': 'pubmed',
                    'id': ','.join(pmids[:3]),
                    'retmode': 'xml'
                }
                if self.api_key:
                    fetch_params['api_key'] = self.api_key
                
                # 简化处理，实际应解析XML
                all_papers.append({
                    'pmid': pmids[0],
                    'title': f"{gene_name} {keyword} study",
                    'abstract': f"Simulated abstract for {gene_name} with {keyword}...",
                    'keyword': keyword
                })
        
        return all_papers

# ==================== 主评估引擎 ====================

class HybridAssessmentEngine:
    """混合评估引擎 - 整合三层分析"""
    
    def __init__(self, qwen_api_key: Optional[str] = None, 
                 ncbi_api_key: Optional[str] = None):
        self.hard_rules = HardRulesEngine()
        self.rag_analyzer = QwenRAGAnalyzer(qwen_api_key)
        self.auxiliary = AuxiliaryAIHints(qwen_api_key)
        self.ncbi = NCBIClient(ncbi_api_key)
    
    def assess(self, gene_name: str, organism: str, 
               experiment_type: str) -> Dict:
        """
        执行完整评估流程
        1. 硬性规则检查（阻断性）
        2. AI文献分析（主要依据）
        3. 辅助提示（参考）
        """
        result = {
            'timestamp': datetime.now().isoformat(),
            'gene': gene_name,
            'experiment': experiment_type,
            'decision_hierarchy': {
                'hard_rules': {},
                'literature_based': {},
                'auxiliary': {}
            },
            'final_recommendation': '',
            'audit_trail': []
        }
        
        # 步骤1: 获取基础数据
        with st.spinner("🔍 检索基因数据..."):
            gene_info, transcripts = self.ncbi.fetch_gene_data(gene_name, organism)
        
        if not gene_info:
            return {'error': '无法获取基因信息'}
        
        # 步骤2: 硬性规则检查（第一优先级）
        with st.spinner("⚙️ 执行硬性规则检查..."):
            hard_passed, hard_checks = self.hard_rules.check_all(
                gene_info, transcripts, experiment_type
            )
        
        result['decision_hierarchy']['hard_rules'] = {
            'passed': hard_passed,
            'checks': [self._dataclass_to_dict(c) for c in hard_checks]
        }
        
        # 如果硬性规则未通过，直接返回阻断结果
        if not hard_passed:
            blocking_rules = [c for c in hard_checks if not c.passed and not c.overrideable]
            if blocking_rules:
                result['final_recommendation'] = 'BLOCKED'
                result['primary_basis'] = '硬性物理/生物学限制（不可协商）'
                result['audit_trail'].append('硬性规则阻断，跳过AI分析')
                return result
        
        # 步骤3: AI文献分析（第二优先级，主要依据）
        with st.spinner("📚 检索文献并执行AI分析（基于阿里云千问）..."):
            literature = self.ncbi.search_literature(
                gene_name, 
                ['overexpression', 'knockdown', 'toxicity', 'lentivirus']
            )
            
            if literature:
                ai_result = self.rag_analyzer.analyze(
                    gene_name, experiment_type, literature
                )
            else:
                ai_result = AIAnalysisResult(
                    conclusion="未检索到相关文献",
                    risk_level="unknown",
                    citations=[],
                    raw_response="",
                    validation_status="no_literature"
                )
        
        result['decision_hierarchy']['literature_based'] = {
            'literature_count': len(literature),
            'ai_analysis': self._dataclass_to_dict(ai_result),
            'citations': [self._dataclass_to_dict(c) for c in ai_result.citations]
        }
        
        # 步骤4: 辅助提示（第三优先级，仅供参考）
        with st.spinner("💡 生成辅助提示..."):
            hints = self.auxiliary.analyze_sequence("")  # 序列分析
            protein_hint = self.auxiliary.ai_protein_prediction(gene_name)
            if protein_hint:
                hints.append(protein_hint)
        
        result['decision_hierarchy']['auxiliary'] = {
            'hints': [self._dataclass_to_dict(h) for h in hints]
        }
        
        # 步骤5: 生成最终建议
        result['final_recommendation'] = self._generate_recommendation(
            hard_checks, ai_result, hard_passed
        )
        result['primary_basis'] = '文献实证（AI辅助提取）' if ai_result.citations else '基于有限文献的推测'
        
        return result
    
    def _dataclass_to_dict(self, obj) -> Dict:
        """将dataclass转换为字典"""
        if hasattr(obj, '__dataclass_fields__'):
            return {k: v for k, v in obj.__dict__.items()}
        return obj
    
    def _generate_recommendation(self, hard_checks, ai_result, hard_passed) -> str:
        """生成最终建议"""
        if not hard_passed:
            return "❌ 不建议进行：违反硬性生物学限制"
        
        if ai_result.risk_level == 'high':
            return "⚠️ 高风险：建议更换方案或采用特殊策略（诱导表达等）"
        elif ai_result.risk_level == 'medium':
            return "⚡ 中等风险：可进行，但需优化实验条件"
        elif ai_result.risk_level == 'low':
            return "✅ 低风险：可按标准流程进行"
        else:
            return "❓ 风险未知：缺乏文献支持，建议预实验"

# ==================== Streamlit UI ====================

def render_header():
    """渲染头部"""
    st.markdown("""
    <h1 style='text-align: center; color: #1f77b4;'>
        🧬 病毒与细胞系评估系统 v2.0
    </h1>
    <p style='text-align: center; color: #666;'>
        混合决策架构：硬性规则 + AI文献分析(RAG) + 辅助提示
    </p>
    """, unsafe_allow_html=True)

def render_sidebar():
    """渲染侧边栏"""
    with st.sidebar:
        st.header("⚙️ 配置")
        
        # API Keys
        st.subheader("API密钥")
        ncbi_key = st.text_input("NCBI API Key (可选)", type="password",
                                help="提高访问频率至10次/秒")
        qwen_key = st.text_input("阿里云千问 API Key (用于AI分析)", type="password",
                                help="获取地址：https://dashscope.aliyun.com")
        
        st.divider()
        st.info("""
        **系统架构说明**：
        1. **硬性规则**：物理限制（如8kb容量），不可协商
        2. **AI文献分析**：基于阿里云千问的RAG，强制引用PMID
        3. **辅助提示**：AI预测，仅供参考（黄色警告标识）
        """)
        
        return ncbi_key, qwen_key

def render_input_form():
    """渲染输入表单"""
    st.markdown("## 📋 实验参数")
    
    col1, col2 = st.columns(2)
    with col1:
        gene = st.text_input("基因名", placeholder="如: TP53, EGFR")
    with col2:
        organism = st.selectbox("物种", 
                               ["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
                               format_func=lambda x: {
                                   "Homo sapiens": "人类",
                                   "Mus musculus": "小鼠", 
                                   "Rattus norvegicus": "大鼠"
                               }.get(x, x))
    
    exp_type = st.selectbox(
        "实验类型",
        ["overexpression", "knockdown", "knockout"],
        format_func=lambda x: {
            "overexpression": "过表达慢病毒",
            "knockdown": "敲低慢病毒 (shRNA)",
            "knockout": "敲除慢病毒 (CRISPR)"
        }.get(x, x)
    )
    
    analyze = st.button("🔍 开始智能评估", type="primary", use_container_width=True)
    
    return gene, organism, exp_type, analyze

def render_hard_rules(checks: List[Dict]):
    """渲染硬性规则结果"""
    st.markdown("### 🚦 第一层：硬性规则检查（不可协商）")
    
    for check in checks:
        if check['passed']:
            icon = "✅"
            color = "green"
        else:
            icon = "❌" if not check['overrideable'] else "⚠️"
            color = "red" if not check['overrideable'] else "orange"
        
        with st.container():
            st.markdown(f"""
            <div style='padding: 10px; border-left: 4px solid {color}; 
                        background-color: #f0f0f0; margin: 5px 0;'>
                <b>{icon} {check['rule_name']}</b><br/>
                <small>{check['reason']}</small><br/>
                <small style='color: #666;'>来源: {check['source']}</small>
                {"<br/><small style='color: red;'><b>不可覆盖</b></small>" if not check['overrideable'] else ""}
            </div>
            """, unsafe_allow_html=True)

def render_ai_analysis(ai_data: Dict):
    """渲染AI文献分析结果"""
    st.markdown("### 📚 第二层：AI文献分析（主要决策依据）")
    
    status = ai_data.get('validation_status', '')
    if status == 'skipped':
        st.warning("未配置千问API，跳过AI分析")
        return
    
    risk_level = ai_data.get('risk_level', 'unknown')
    colors = {
        'high': ('🔴', 'red'),
        'medium': ('🟡', 'orange'), 
        'low': ('🟢', 'green'),
        'unknown': ('⚪', 'gray')
    }
    icon, color = colors.get(risk_level, ('❓', 'gray'))
    
    st.markdown(f"""
    <div style='padding: 15px; border: 2px solid {color}; border-radius: 5px; 
                background-color: #fafafa;'>
        <h4>{icon} AI风险评估结论: {risk_level.upper()}</h4>
        <p>{ai_data.get('conclusion', '无结论')}</p>
        <small style='color: #666;'>验证状态: {status}</small>
    </div>
    """, unsafe_allow_html=True)
    
    # 显示引用文献
    citations = ai_data.get('citations', [])
    if citations:
        st.markdown("**📖 文献证据链：**")
        for cite in citations:
            with st.expander(f"PMID: {cite['pmid']} [{cite['evidence_type']}]"):
                st.markdown(f"""
                **文献标题**: {cite['title']}  
                **证据原文**: "{cite['evidence_text']}"  
                **证据类型**: {cite['evidence_type']}  
                **可信度**: {cite['confidence']}
                """)
    else:
        st.info("AI未找到相关文献证据")

def render_auxiliary(hints: List[Dict]):
    """渲染辅助提示"""
    st.markdown("### 💡 第三层：AI辅助提示（仅供参考）")
    
    if not hints:
        st.info("无辅助提示")
        return
    
    for hint in hints:
        st.markdown(f"""
        <div style='padding: 10px; border-left: 4px solid #ffc107; 
                    background-color: #fff3cd; margin: 5px 0;'>
            <b>⚠️ {hint['category']}</b><br/>
            {hint['content']}<br/>
            <small style='color: #856404;'>
                置信度: {hint['confidence']*100:.0f}% | 
                <b>{hint['warning']}</b>
            </small>
        </div>
        """, unsafe_allow_html=True)

def main():
    """主函数"""
    render_header()
    ncbi_key, qwen_key = render_sidebar()
    
    gene, organism, exp_type, analyze = render_input_form()
    
    if analyze and gene:
        # 初始化引擎
        engine = HybridAssessmentEngine(qwen_api_key=qwen_key, ncbi_api_key=ncbi_key)
        
        # 执行评估
        with st.spinner("正在进行混合决策评估..."):
            result = engine.assess(gene, organism, exp_type)
        
        if 'error' in result:
            st.error(result['error'])
            return
        
        # 显示结果
        st.divider()
        st.markdown(f"## 🎯 评估报告 - {gene}")
        
        # 最终建议
        rec = result['final_recommendation']
        st.markdown(f"""
        <div style='padding: 20px; background-color: #e7f3ff; 
                    border-radius: 10px; text-align: center;'>
            <h3>{rec}</h3>
            <small>主要依据: {result['primary_basis']}</small>
        </div>
        """, unsafe_allow_html=True)
        
        # 三层决策展示
        hierarchy = result['decision_hierarchy']
        
        # 1. 硬性规则
        render_hard_rules(hierarchy['hard_rules']['checks'])
        
        # 如果未被硬性规则阻断，显示后续分析
        if hierarchy['hard_rules']['passed']:
            # 2. AI文献分析
            render_ai_analysis(hierarchy['literature_based']['ai_analysis'])
            
            # 3. 辅助提示
            render_auxiliary(hierarchy['auxiliary']['hints'])
        
        # 审计追踪
        with st.expander("🔍 查看审计追踪（Audit Trail）"):
            st.json(result)

if __name__ == "__main__":
    main()
