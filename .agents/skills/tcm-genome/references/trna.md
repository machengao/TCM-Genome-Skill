# 从注释提取 tRNA

1. 按[药典基原提取与归并](pharmacopoeia.md)查明基原。按[多库组装检索](genome-sources.md)查找用户指定的核、叶绿体或线粒体参考，核对选用版本、材料和组装成分。用户只要组装目录时保留元数据。用户要求 tRNA 时下载所选参考及配套注释。
2. 每套基因组独立运行脚本。NCBI 细胞器序列登录号通过 EFetch 取 GenBank 记录。版本化的 NCBI 组装登录号通过 Datasets 取 FASTA、GFF3 和序列报告。GWH 登录号通过组装详情获取 FASTA 和 GFF3。也支持本地 GenBank 或同一组装的 FASTA/GFF3。
3. 有任何可用 tRNA 原注释的参考走坐标提取。不精确、缺链方向或假基因标签的条目留在表中。版本、序列 ID、注释坐标或下载不一致时处理错误并核对来源。只有注释缺失或已解析且没有可用 tRNA 位点才启动预测。疑似注释不全但已含可用位点时，如实报告原注释范围。

## 命令

从 Skill 目录运行。需要 Python 3.11+ 和 Biopython，可用 `python -m pip install "biopython>=1.84,<2"` 安装。远程 GenBank 和 GWH 需要 curl，远程 NCBI 组装需要 NCBI Datasets CLI。回退需要支持 `-j` GFF3 输出的 tRNAscan-SE 2.x。植物核基因组使用真核模式，植物细胞器使用细胞器模式。处理其他类群的细胞器时先确认适用的预测模型。

本机 Windows 已配置 D 盘环境时，先在同一 PowerShell 会话运行 `& 'D:\Bioinformatics\Environments\TCM-Genome-Skill\Use-TcmGenome.ps1'`。启用环境后，预测会按 `TCM_TRNASCAN_WSL_DISTRO` 调用指定的 WSL 发行版。输入和输出需放在该发行版能访问的 Windows 路径中。也可以用 `--scanner` 显式指定本机可执行文件。

```text
python scripts/trna.py --accession NC_006290.1 --compartment chloroplast --output ginseng-cp
python scripts/trna.py --accession MZ389476.1 --compartment mitochondrion --output ginseng-mt
python scripts/trna.py --fasta genome.fna --gff genomic.gff --compartment nuclear --output trna-nuclear
```

已用[批量下载](genome-sources.md#批量下载参考)保存参考时，用 `--accession` 和 `--download-dir` 指定登录号和缓存根目录，可直接复用序列、注释和分区报告。`--threads` 设置单次运行的预测线程，默认 1。需要同时预测多个参考时，可分别启动命令并为每次指定独立结果目录，根据可用 CPU、内存和临时磁盘空间分配线程。

本地 FASTA 和 GWH 组装若包含多个分区，先确认序列清单，再用可重复的 `--sequence-id` 选择目标序列。NCBI 组装优先使用序列报告识别细胞器成分。保留脚本的输入 SHA-256 与源版本。

## 结果和解释

- `trna.tsv` 逐位点提供序列 ID、链、1-based 闭区间、名称、可用性和证据来源。`trna.fasta` 提供按转录方向、按已标注外显子拼接的 DNA 序列。保留每个独立坐标，包括重复拷贝和跨区室的相同序列。
- `summary.json` 记录执行路径和输入校验值。预测另保留工具原始表和 GFF3。直接提取的是发布的注释，预测走计算证据。成熟 RNA 末端加工、CCA、RNA 编辑和实验确认在后续步骤单列。
- tRNAscan-SE 成功完成且明确报告零命中时，保留零条预测的结果与摘要。来源 GFF3 仍按正常规则校验。
- 保留注释位点数和实际导出序列数。原数据库没有标记的残缺位点或误注释仍可能被直接提取，人工复核结果应独立记录。
