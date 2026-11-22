# 公共复现档案安全与内容审计

审计日期：2026-11-22  
审计范围：本目录中的拟公开文件；原研究目录及投稿文件保持只读。

## 结论

本地 staging 已通过内容、安全、路径、语法、图形冒烟测试和数字一致性检查。档案仅包含复现代码、配置模板、公开来源元数据、小型派生汇总表和 Figure 2–4 的精确源数据与绘图脚本；未包含原始或大规模 GWAS 文件、受限数据、个人信息、密钥、内部提示词或稿件二进制文件。

公开推送与归档条件已满足：作者确认代码采用 MIT，原创派生表格、元数据和文档采用 CC BY 4.0；Git 提交身份为 `duanxian0301 <duanxian0301@gmail.com>`。作者授权 `v1.0.0` GitHub/Zenodo 归档，creator 仅填写 `Jing Shen`，单位、ORCID、基金、通讯作者标识及论文 DOI 按作者要求不填写。

## 已纳入

- 冻结的 UKB 61×30 因子载荷和 13 个保留因子的 QC/派生摘要；
- Karjalainen 61×15 独立 GFA 载荷及结构对应摘要；
- S2/S3、位点审计和 LAVA 的小型正式汇总；
- Figure 2–4 精确源数据和绘图脚本；
- 核心 GFA、GLS、独立结构比较、单性状基线、位点/LAVA 与审计脚本；
- 61 性状公开 accession、输入清单、文件和字段字典；
- 环境版本、配置示例、分析说明、判定定义与已知限制。

## 已排除

- UK Biobank 受控资源、个体级数据及任何可识别信息；
- 大规模、harmonized 或 munged GWAS summary statistics；
- LD reference、PLINK reference、LAVA correlation file、缓存及序列化对象；
- 稿件、补充材料、期刊 PDF、内部提示词、GROK/BOT 记录及临时日志；
- NDD、Aging 和 MVMR 可行性工作，因为这些不属于稿件阳性结果链；
- 两个带有本地下载/删除逻辑的过滤脚本，避免公开不可移植和潜在破坏性操作。

公开 GWAS 的 accession 和链接见 `metadata/input_manifest.tsv`；读者应从原始数据仓库获取并遵守其使用条款。

## 已完成的检查

- 发布检查：16/16 PASS；包括必需文件、单文件大小、禁止二进制、UTF-8、本地路径、密钥/内部词、LFS 指针、矩阵维度、保留因子数和图源哈希。
- 稿件/补充数字核对：35/35 PASS（`validation/submission_numeric_qa.tsv`）。
- Python 语法：10 个文件，0 失败；R parse-only：8 个文件，0 失败。
- Figure 2–4：均由纳入的 TSV 成功生成 PDF、SVG 和约 600 dpi PNG。渲染产物为本地冒烟测试输出并由 `.gitignore` 排除；公开仓库保留源数据和脚本。
- 远端检查：目标 GitHub 仓库可访问且当前为空；未发现需合并的既有提交。

## 关键科学边界

- 固定 UKB 载荷在 Karjalainen 数据上的资源特异 GLS 投影，与 Karjalainen 独立拟合的 GFA 因子不是同一 estimand。
- 结构对应和固定投影一致性不能写成完全独立复现。
- 单性状比较为描述性基线；共享 GWAS 输入使比较相关，不能据此宣称独立统计增量。
- 本档案不提供从零开始重跑受限或大规模上游分析所需的 GWAS/LD 文件，只提供真实结果链、配置接口和安全派生产物。

## 执行命令（不含新科学分析）

```text
python validation/check_release.py
python validation/smoke_test_figures.py
Python ast.parse over staged *.py
Rscript --vanilla parse(file=...) over staged *.R
git ls-remote --symref https://github.com/duanxian0301/NMR_Latent_Transport_Disease.git HEAD
```

## Gate

GitHub 推送 gate：PASS。许可证、提交身份、发布检查和 manifest 均已完成。  
Zenodo/release gate：AUTHORIZED。`CITATION.cff` 和 `.zenodo.json` 采用作者批准的最小元数据；仅在实际归档成功后记录真实 DOI。
