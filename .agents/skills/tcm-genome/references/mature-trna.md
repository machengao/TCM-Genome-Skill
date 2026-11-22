# 从已有位点重建成熟 tRNA

输入为 `trna.py` 的 tRNAscan-SE 预测结果目录，含 `trna.fasta`、`trna.tsv`、`summary.json`。保留输入位点 ID、同序列的不同拷贝和原有区室。自动读取输入目录中的 `reference_annotation/`。也可用 `--annotation` 指定同版本的既有注释目录。注释重叠、坐标差异及类型冲突仅作为对照信息。用 `--genome` 提供同版本未压缩的基因组 FASTA，可逐碱基复核两套输入。只有裸 FASTA 时，应先找回位点表和参考坐标。

从技能目录运行，使用 Python 3.11+、Infernal 的 `cmalign` 和 tRNAscan-SE 的通用 `TRNAinf.cm`。模型包含受体茎七个外侧配对位点、鉴别碱基和反密码子标记。Linux 发行版若安装在其他位置，用 `--model` 指定路径。在 Windows 上可加 `--wsl-distro TCM-Genome-Ubuntu`，模型路径使用该 WSL 发行版中的路径。

```sh
python scripts/mature_trna.py --input previous-prediction --annotation previous-annotation --genome reference.fna --output mature-results --label "物种 · 区室 · 参考号"
```

1. 根据上游提取方式还原 DNA 坐标。GenBank 分段遵循 feature 的转录顺序，GFF 分段按链方向排列。现有 FASTA 已经按标注外显子拼接。比对同参考同链的注释时，检查重叠关系、共有碱基和端点。成熟序列主体沿用扫描预测的外显子、边界和链方向。
2. 对剪接后序列进行全局结构比对。枚举输入末端覆盖 CCA 的 0、1、2、3 个碱基的解释，并验证结构位点 1–7、66–72 和其后的鉴别碱基。仅根据预测序列评估成熟边界。末端或反密码子无法定位、存在多解时，该条记录进入复核表。
3. 保留结构上归于 CCA 的输入碱基，按 A/C/G/U 输出 5′ → 3′ 序列，并补齐余下的 CCA。添加只表示相对于选定边界的加工推断。His 单独核对模型位置 1 前的 G−1；细胞器保留有坐标支持的 G−1，未找到时待复核。典型核 tRNA-His 推断补 G−1。真实样本中的末端、修饰及编辑状态仍取决于转录本证据。
4. 有已核实的同物种位点替换时，传入 `--edits edits.json`。数组元素包含 `locus_id`、版本化 `seq_id`、整数 `genome_position`、按转录方向给出的单碱基 `from` 和 `to`（A/C/G/U）及 `source`（原始记录或论文及位点）。脚本先核对参考碱基，再用编辑型参与结构比对。无证据时使用参考碱基，编辑状态为 `unassessed`。碱基化学修饰不在四字母序列中展开。

## 查阅结果

打开 `report.html` 查阅各位点的结果和复核原因。`mature_trna.fasta` 只收录边界有唯一解释的计算候选。`mature_trna.tsv` 包含所有输入位点、原始序列、边界依据、CCA 来源及反密码子起始位置。提供参考基因组时，`downstream_match_nt` 表示所选边界外与拟补尾序列前缀相同的连续碱基数。所有输出坐标从 1 开始，反密码子位置在成熟序列上计数。

`base_map.tsv` 把每一成熟碱基关联到参考坐标、链和加工来源，推断添加的碱基其基因组坐标为空。`alignment_id` 对应 `alignment_inputs.fasta` 和 `alignment.sto` 中获选的结构比对序列，即去掉 G−1 和 CCA 后的主体。`summary.json` 记录输入校验值、模型校验值及 Infernal 版本。位点数与不同序列数分别统计。

`model_loops` 给出成熟体的 D 环、反密码子环及 T 环范围。`anticodon_cut_after` 记录反密码子环内五个结构位点在成熟序列上的切后坐标。读取[片段建库](trna-fragments.md)继续处理。

非典型受体茎、严重截短的细胞器 tRNA 和复杂编辑可能无法用此通用模型重建。可用原始结构或转录本数据复核具体位点。方法依据见 [tRNAscan-SE 2.0](https://doi.org/10.1093/nar/gkab688)、[成熟体参考的 CCA 与 G−1 处理](https://doi.org/10.1038/srep41184) 和[植物细胞器 tRNA-His 加工研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC2995067/)。
