# Lentivirus-and-Cell-Line-Assessment-v2
This is an application used for assessing cell line and lentivirus.
# 病毒与细胞系评估系统 v2.5

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-app-url.streamlit.app)

基于混合决策架构的智能评估工具，整合 **NCBI Clinical Tables API**（基因自动完成）与 **阿里云千问**（AI文献分析），为慢病毒包装和细胞系构建提供科学的风险评估。

## 🚀 核心特性

### 1. 智能基因输入（新增）
- **实时自动完成**：基于 NCBI Clinical Tables API，输入2个字符即提示官方基因符号
- **多物种支持**：人类（Human）、小鼠（Mouse）、大鼠（Rat）
- **基因信息卡片**：显示染色体位置、基因类型、NCBI链接
- **常用基因快捷栏**：TP53、EGFR、KRAS等一键选择
- **300ms防抖机制**：避免频繁API调用，提升响应速度

### 2. 三层决策架构
第一层：硬性规则（物理/生物学限制）
载体容量检查（8kb物理极限）
必需基因检查（DepMap数据库）
毒性基因检查（文献数据库）

第二层：AI文献分析（RAG-based）
NCBI E-utilities 检索真实文献
阿里云千问语义分析
强制PMID引用与验证

第三层：辅助提示（算法预测）
序列特征分析（GC含量等）

### 3. 安全加固
- **输入验证**：严格基因名格式检查（防注入/XSS）
- **XSS防护**：所有HTML输出强制转义
- **异常脱敏**：使用 UUID 替代 MD5，错误信息不暴露内部细节
- **引用验证**： difflib.SequenceMatcher 验证AI引用真实性（相似度>0.6）


# 运行
streamlit run app.py
