# 从预测成熟 tRNA 建立 tRF 与 t-half 库

从技能目录运行。输入使用成熟重建模块输出的整套结果，其中 `summary.json` 的 `sequence_source` 为 `trnascan_prediction`。历史结果若曾以发布注释修正末端或外显子，应从原预测目录重新成熟化。

```sh
python scripts/trna_fragments.py --input mature-results --output fragment-library
```

- 对每条通过重建的成熟 tRNA，沿 5′ → 3′ 方向穷举长度 19 到 25 nt 的连续窗口。贴住母体 5′ 端标为 `tRF-5`，贴住 3′ 端标为 `tRF-3`，两端均在内部标为 `i-tRF`。另外记录两端对应的 D 环、反密码子环、T 环、CCA、G−1 或其他区域。对于长度至少 25 nt 的母体，窗口总数为 `7 × (母体长度 − 21)`。
- t-half 按成熟体自身的结构比对注释切割。以反密码子三碱基对应的结构位置 34–36 为中心，同时取相邻的结构位点 33、37，分别在五个位点之后切割。结构编号和成熟序列下标分别记录。实际切点从 `anticodon_cut_after` 读取，采用 1-based 成熟坐标。切在第 k 个碱基后，对应 5′ 半段 [1,k] 和 3′ 半段 [k+1,L]，共用同一 `pair_id`。成熟体上无法完整映射五个位点时仅生成短片段，位点列入 `unresolved_cut_loci`。
- 两个库的 TSV 保留每个母体、窗口位置、长度、类别和序列。打开 `report.html` 可逐位点查看 t-half 配对。两个 FASTA 分别按 RNA 序列去重，以序列的 SHA-256 生成稳定 ID，查询 TSV 的 `sequence_id` 可找回全部来源。坐标回溯基因组时，用 TSV 的 `parent_id` 和成熟包中的 `base_map.tsv`。证据字段 `enumerated_candidate` 表示计算枚举，生物样本检出情况需另行取证。
