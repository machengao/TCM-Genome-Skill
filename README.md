# TCM-Genome-Skill

An Agent Skill and a set of command-line tools for following Chinese Pharmacopoeia entries from source organisms to public genome references. The workflow also extracts or predicts tRNA loci from selected references and helps trace published evidence for genes involved in natural product biosynthesis.

## What it does

- **Read the pharmacopoeia.** Use the official online 2025 Chinese Pharmacopoeia, Volume I, or text supplied by the user. Read source-organism descriptions for medicinal materials and ingredient sections for formulas. Keep links to the original entries and record coverage gaps.
- **Build a species search list.** Check scientific names and synonyms against taxonomic sources. Preserve every medicinal material–source relationship and search sources that resolve to the same taxon together. Keep unresolved identifications visible.
- **Find genome references.** Search NCBI and GWH for each taxon, plus any sources explicitly requested. Request nuclear assemblies, mitochondrial genomes, and chloroplast genomes separately. For nuclear assemblies, the workflow can consult sources such as ENA, MycoCosm, or Phytozome when the primary sources yield no confirmed match.
- **Follow up in the literature.** The search report lists genome types that still need a reference under `literature_followup`. Review the original papers and trace their assemblies to accession records or public files. When the available assemblies are comparably suitable, journal impact can help prioritize which studies to review first.
- **Download selected references.** Save reusable bundles by accession, with checksums and a local progress report. The batch downloader starts with one transfer and adjusts concurrency based on observed throughput, up to four downloads by default.
- **Extract or predict tRNA loci.** Use compatible database annotations when they contain usable tRNA loci. When there are none, run tRNAscan-SE and retain its original outputs. A successful scan with zero hits produces a valid zero-result report.

## Use with an AI agent

Open the repository root as your Codex workspace. Codex discovers the skill under `.agents/skills/tcm-genome/`. Type [`$tcm-genome`](.agents/skills/tcm-genome/SKILL.md) followed by a request. [Codex skill discovery](https://learn.chatgpt.com/docs/build-skills#where-codex-loads-local-skills)

> $tcm-genome Find genome references for every fungal medicinal material in the 2025 Chinese Pharmacopoeia, Volume I. Check NCBI and GWH, report coverage, and export a metadata catalog to Excel.

> $tcm-genome Find nuclear, mitochondrial, and chloroplast genome references for ginseng. List the accession versions and evidence for each compartment.

> $tcm-genome Extract tRNA loci from the ginseng chloroplast reference NC_006290.1. Use the published annotation where available and run tRNAscan-SE if no usable tRNA annotations are present.

> $tcm-genome Find the published evidence for ADS and CYP71AV1 in artemisinin biosynthesis in Artemisia annua. Include traceable gene identifiers and the reactions supported by each source.

For another agent that supports Agent Skills, load the entire `tcm-genome` directory in a location that the agent recognizes. The agent needs terminal access to run the scripts and web access to follow up in the literature.

## Run the scripts locally

The scripts run directly from a terminal. Clone the repository into a folder of your choice, use Python 3.11 or later, and put `curl` on your PATH for online queries. If your system uses `python3`, use that command to create the virtual environment.

```sh
git clone https://github.com/machengao/TCM-Genome-Skill.git
cd TCM-Genome-Skill
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on Linux or macOS, or `.\.venv\Scripts\Activate.ps1` in Windows PowerShell.

- For **NCBI assembly downloads**, install [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) and put `datasets` on your PATH.
- For **tRNA extraction and prediction**, install Biopython as shown below. Prediction also requires [tRNAscan-SE 2.x with Infernal](https://github.com/UCSC-LoweLab/tRNAscan-SE/blob/master/INSTALL) and support for the `-j` option.
- For **prediction from Windows Python**, install the scanner in [WSL](https://learn.microsoft.com/en-us/windows/wsl/install). Find the distribution name with `wsl --list --verbose` and set it in the PowerShell session that runs the scripts. For a distribution named Ubuntu, run `$env:TCM_TRNASCAN_WSL_DISTRO = 'Ubuntu'`. When running Python inside Linux or macOS, make `tRNAscan-SE` available on your PATH or use `--scanner`.

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

Choose your own input and output paths and quote paths that contain spaces. Relative paths in the examples resolve from the directory where you run the command. Pass the same download directory to `--download-dir` to reuse the downloaded references, keeping the accession subdirectories and their manifests together.

## Results and coverage

The reader and search scripts write JSON reports. An agent or spreadsheet tool can use those reports to build an Excel catalog with source links and query status. For tRNA processing, `trna.tsv` records loci, coordinates, strands, and evidence, `trna.fasta` contains the corresponding genomic DNA sequences in transcriptional orientation, and `summary.json` records the processing route and input checksums.

The fungal medicinal-material pilot reviewed all 616 entries in the medicinal materials and decoction pieces section of Volume I. The ginseng chloroplast test yielded 37 sequences from the published annotation and 33 sequences in a separate prediction run on FASTA alone. Database coverage, annotation availability, and taxonomic resolution vary by organism. Review the per-source query status before treating a catalog as complete.

See the skill references for [pharmacopoeia extraction](.agents/skills/tcm-genome/references/pharmacopoeia.md), [genome search and downloads](.agents/skills/tcm-genome/references/genome-sources.md), and [tRNA processing](.agents/skills/tcm-genome/references/trna.md).
