# 多库组装检索

## 从白话确定必查来源

- 每个目标物种都检索 NCBI Datasets 和 GWH，即使其中一个已经命中。
- 从用户请求提取附加必查库及适用范围。适用于整批的库用 `--require` 指定；按类群或物种指定的库写入相应基原的 `required_sources`。例如“所有基原加查 ENA，真菌也必须查 MycoCosm”对应全局 ENA，真菌记录另设 MycoCosm。
- NCBI 和 GWH 都没有确认组装时按类群补查。真菌查 MycoCosm 和 ENA，植物查 Phytozome 和 ENA，藻类查 PhycoCosm 和 ENA，微生物查 IMG 和 ENA，动物及其他类群查 ENA。已明确异名时，用接受名和异名补查。矿物等不含基原生物的条目记录为不适用。
- 对于尚无脚本适配的指定库，尝试其公开 API 或官方网页，并记录实际查询词和结果。接口或网页不可用时，将该库标为待查或访问失败。继续交付其余来源的结果。

## 执行检索

从 Skill 目录执行。脚本依赖 Python 3.10+ 与 curl，输入为已完成药典基原提取和名称核对的 JSON 列表。把药材名、条文链接等信息留在各自的输入行，报告会保留。

```json
[
  {
    "name": "Polyporus umbellatus",
    "aliases": ["Dendropolyporus umbellatus"],
    "tax_id": 158314,
    "group": "fungi",
    "required_sources": ["MycoCosm"]
  }
]
```

```text
python scripts/genome_sources.py taxa.json --require ENA --output genome-report.json
```

- 用户明确说“线粒体”时传 `--sequence-type mitochondrion`，“叶绿体”传 `--sequence-type chloroplast`。两个都要就同时传入；还要常规组装则再加 `--sequence-type nuclear`。未指定时按常规组装检索。
- 细胞器检索使用 [NCBI Datasets Organelle](https://www.ncbi.nlm.nih.gov/datasets/organelle/) 的类型字段和 GWH 直接提交记录详情页的 Genome Representation 标注。GWH 镜像目录的部分组装由 NCBI 细胞器入口另查。保留登录号、类型依据和出处。跨库合并细胞器记录时应核对登录号或明确的交叉引用。其他指定库的细胞器查询交由公开页面补查。
- 常规组装记录有时包含细胞器序列。只有检查具体序列清单后才能确认核序列和细胞器边界。标为 `assembly_unverified` 的记录尚未做组件级核对。当前流程仅查询元数据。

- `group` 可取 `fungi`、`plants`、`animals`、`algae`、`microbes` 或 `other`。只有已核实的同物异名写入 `aliases`。争议名称写入 `candidate_names`，按该名命中的记录等待核实。只有有证据的分类编号写入 `tax_id`。
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/get-genome-metadata/) 和 [ENA Portal](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access.html) 读组装元数据。GWH 用[公开组装目录](https://ngdc.cncb.ac.cn/gwh/browse/assembly)的网页查询发现登录号，再用[官方组装详情 API](https://ngdc.cncb.ac.cn/gwh/api_documents)核对直接提交记录。网页查询格式变更时按查询失败处理。
- [JGI Data Portal API](https://files.jgi.doe.gov/apidoc/) 提供 MycoCosm、Phytozome、PhycoCosm 和 IMG 的数据集查询。PhycoCosm 文件接口无结果时，再查[公开发布目录](https://phycocosm.jgi.doe.gov/phycocosm/home/releases)。数据集线索须回查平台或论文，确认有对应基原的核基因组和组装版本后再计为组装。

## 核对并交付

- 逐库保留检索词、结果和失败状态。`hit` 表示检得名称匹配的组装，`candidate` 表示有待核资料，`no_hit` 表示相应查询成功完成且未命中。`error`、`access_required`、`rate_limited`、`manual_required` 和 `not_run` 都表示覆盖有缺口。报告可在部分来源失败时保留已取得的数据。
- 核对提交物种、菌株、组装级别、登录号版本及出处。显式配对的 GCA/GCF 和跨库共享的同一版本登录号可合并，保留各库证据。只有未带版本的编号、样本号或名称相同时，先核对记录版本再归并。按用户指定类型呈现组装与细胞器序列。原始读段和转录组不计为基因组组装。
- 对药典历史用名、错用名和分类争议复查样本证据。覆盖核对按每个必查库与待核项报告。只有所有所需查询都完成，才标记对应范围已经查全。

## 批量下载参考

从组装检索结果选定版本化登录号。对已确认跨库配对的同一参考选一个下载入口，保留其他库的出处。用户明确要求下载时，在同一个 PowerShell 会话启用 D 盘环境，再从 Skill 目录执行：

```powershell
& 'D:\Bioinformatics\Environments\TCM-Genome-Skill\Use-TcmGenome.ps1'
python scripts/genome_download.py --accession-file selected-accessions.txt --output D:\Bioinformatics\Databases\PublicDatasets
```

清单每行写一个登录号，也可直接把登录号作为位置参数传入。每个已完成的参考保存到输出目录下的同名子目录，配有资产清单和 SHA-256。再次运行时复用原目录的已完成参考，失败项可重试。控制台会显示本次 JSON 进度报告的位置，也可用 `--report` 指定新路径。

下载从一个任务起步。脚本每隔约五秒用实际收到的字节测总吞吐，逐档试增并发。提升不足或下载报错时降档，让已开始的下载跑完。默认上限为四个，可用 `--max-downloads` 调整。短文件的测速样本不足时保守运行。本地 Python 进程负责测量和排队。
