# Oxeiptosis v04 model reproduction package

本地复算包针对本轮论文的数值核验；不是从公共数据库全新下载、从零部署的完整研究流水线。

## Run

需要 R 和 survival。审查环境为 R 4.6.1、survival 3.8-6，详情见 expected/R_session_v04.txt。校验脚本仅需 Python 标准库。解压后保留本目录结构；以下命令在本目录运行。Windows 请先进入本目录后使用下面的相对脚本名；当前主机直接向 Rscript 传递含中文的绝对路径会出现区域设置/路径解码错误。

```text
python verify_inputs.py
Rscript reproduce.R smoke
Rscript reproduce.R full
```

smoke 模式拟合三队列、五构造、M0/M1/M2 共 45 项 OS 模型。full 模式另外拟合 300 项复发风险集诊断及 6,000 项固定成员的随机背景 M1 模型。默认结果写入 run_smoke 或 run_full；第二个参数可指定新的结果目录。输入和 expected 只读。输出目录会覆盖同名本地运行结果，重新运行需留档时请指定不同目录。

当前移植脚本的 smoke 模式已与本轮 R 结果逐行核对。full 模式所用计算主体与 audit_source/verify_models.R 相同；6,345 项完整计算已在原工作目录实际执行，移植包全量模式未另行再跑。随机矩阵使用 gzip 仅作无损压缩，6,000 个成员组合未重抽。

## Contents and interpretation

- inputs：三队列 OS 模型帧、复发诊断帧、随机成员注册表和随机评分矩阵。sample_id 是公开数据的编码标识，不是新增临床受试者信息。
- expected：本轮 R 完整结果与版本记录。30 项复发模型带有收敛警告，保留作诊断，不能作为有效推断。6,000 项随机背景 M1 模型无收敛警告。
- audit_source：本地原矩阵重建、既往解析函数来源、论文与图表生成代码。保留了工作区相对/绝对依赖，不能直接当作此包内的一键全流程命令。
- provenance：冻结输入、额外来源、原矩阵重建依赖哈希、数值和版面记录。数据源仍留在原路径，不因打包迁移或重写。
- manifest.json：本包发布时文件的 SHA256（不含 manifest 自身与后续 run_* 结果）。

OS 主评分 OX4 为等权 KEAP1/PGAM5/AIFM1/OTUD1。各基因标准化 ddof=0，合成分数标准化 ddof=1。GSE39582 参考总体为 566，OS 拟合总体为 561；不要在复算时把它改成拟合子集 SD 后仍沿用本结果。

这里复现的是已声明构造的关联分析，不是 oxeiptosis 功能活性的校准。复发风险集保留来源与零时间假设；TCGA 新肿瘤事件不视同 RFS。核验算法实现不等于独立人工评定。论文定向文献核验不等于全量系统综述。

本目录是待公开发布的复现包候选版本；最终创建者、许可证、引用信息和公共仓库/Zenodo 标识尚未填入。原始公共矩阵来源与访问说明见上一级 supplementary_methods_v04.md 和参考文献。正式发布前应由作者核对 Data Availability、公共数据使用条款和最终版本号。
