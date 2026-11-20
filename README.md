# TCM-Genome-Skill

A Codex skill for building a traceable link between Chinese Pharmacopoeia entries, the organisms they specify, and public genome references. It also extracts or predicts tRNA loci from selected references and can help trace published evidence for genes involved in natural product biosynthesis.

## What it does

- **Read the pharmacopoeia.** Use the official online 2025 Chinese Pharmacopoeia, Volume I, or text supplied by the user. Read source-organism descriptions for medicinal materials and ingredient sections for formulas. Keep links to the original entries and record coverage gaps.
- **Build a species search list.** Check scientific names and synonyms against taxonomic sources. Preserve every medicinal material–source relationship and search sources that resolve to the same taxon together. Keep unresolved identifications visible.
- **Find genome references.** Search NCBI and GWH for each taxon, plus any sources explicitly requested. Request nuclear assemblies, mitochondrial genomes, and chloroplast genomes separately. For nuclear assemblies, the workflow can consult sources such as ENA, MycoCosm, or Phytozome when the primary sources yield no confirmed match.
- **Download selected references.** Save reusable bundles by accession, with checksums and a local progress report. The batch downloader starts with one transfer and adjusts concurrency based on observed throughput, up to four downloads by default.
- **Extract or predict tRNA loci.** Use compatible database annotations when they contain usable tRNA loci. When there are none, run tRNAscan-SE and retain its original outputs. A successful scan with zero hits produces a valid zero-result report.

## Use in Codex

Open the repository in Codex and start with one of these prompts using [`$tcm-genome`](.agents/skills/tcm-genome/SKILL.md).

> $tcm-genome Find genome references for every fungal medicinal material in the 2025 Chinese Pharmacopoeia, Volume I. Check NCBI and GWH, report coverage, and export a metadata catalog to Excel.

> $tcm-genome Find nuclear, mitochondrial, and chloroplast genome references for ginseng. List the accession versions and evidence for each compartment.

> $tcm-genome Extract tRNA loci from the ginseng chloroplast reference NC_006290.1. Use the published annotation where available and run tRNAscan-SE if no usable tRNA annotations are present.

> $tcm-genome Find the published evidence for ADS and CYP71AV1 in artemisinin biosynthesis in Artemisia annua. Include traceable gene identifiers and the reactions supported by each source.

## Run the scripts locally

Use Python 3.11 or later and put `curl` on your PATH. Install Biopython for tRNA processing. Install NCBI Datasets CLI to download NCBI assembly packages and tRNAscan-SE 2.x with Infernal if prediction is needed. On Windows, install the scanner in WSL and set `TCM_TRNASCAN_WSL_DISTRO` to the name of the distribution that provides it.

From the repository root, enter the skill directory and read one pharmacopoeia entry.

```sh
cd .agents/skills/tcm-genome
python -m pip install "biopython>=1.84,<2"
python scripts/pharmacopoeia.py materia-medica --title 人参 --output herbs.json
```

The genome search script takes a JSON list of taxa whose identities have been checked. For example, save the following as `taxa.json`.

```json
[
  {"name": "Panax ginseng", "tax_id": 4054, "group": "plants"}
]
```

Use the metadata report to choose a reference accession. This example shows the commands for a chloroplast lookup, downloading the selected reference, and processing its tRNA loci.

```sh
python scripts/genome_sources.py taxa.json --sequence-type chloroplast --output genome-report.json
python scripts/genome_download.py NC_006290.1 --output reference-cache
python scripts/trna.py --accession NC_006290.1 --download-dir reference-cache --compartment chloroplast --output ginseng-chloroplast-trna
```

Use `--accession-file` to download a list of selected accessions and `--max-downloads` to set the concurrency ceiling. Use `--threads` to set the thread count for an individual prediction.

## Results and coverage

The reader and search scripts write JSON reports. Codex can turn those reports into an Excel catalog with source links and query status. For tRNA processing, `trna.tsv` records loci, coordinates, strands, and evidence, `trna.fasta` contains the corresponding genomic DNA sequences in transcriptional orientation, and `summary.json` records the processing route and input checksums.

The fungal medicinal-material pilot reviewed all 616 entries in the medicinal materials and decoction pieces section of Volume I. The ginseng chloroplast test yielded 37 sequences from the published annotation and 33 sequences in a separate prediction run on FASTA alone. Database coverage, annotation availability, and taxonomic resolution vary by organism. Review the per-source query status before treating a catalog as complete.

See the skill references for [pharmacopoeia extraction](.agents/skills/tcm-genome/references/pharmacopoeia.md), [genome search and downloads](.agents/skills/tcm-genome/references/genome-sources.md), and [tRNA processing](.agents/skills/tcm-genome/references/trna.md).
