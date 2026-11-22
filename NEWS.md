---
title: "changeLog"
---

# SPACCa [<img src="https://sebastian-gregoricchio.github.io/SPACCa/SPACCa_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/SPACCa)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/SPACCa?style=social)](https://github.com/sebastian-gregoricchio/SPACCa/fork)


<br>

#### [v0.1.3](https://github.com/sebastian-gregoricchio/SPACCa/releases/tag/0.1.3) - November 22<sup>th</sup> 2026
* **Paired-end data: filtered BAMs keep complete read pairs only** (DNA-mapping and ChIP-analyses workflows). The MAPQ filter acts on single reads, so a read whose mate was removed stayed in the BAM as an orphan that still carried the old proper-pair flag, mate position and TLEN: R1 and R2 counts differed in the filtered BAMs, and MACS (`-f BAMPE`) and `bamCoverage --extendReads` treated orphans as complete fragments. After the MAPQ filter, mates are now re-paired (`samtools collate` + `samtools fixmate -m`) and only primary, properly paired alignments are kept (`samtools view -f 2 -F 2828`) before duplicate marking. Secondary and supplementary alignments are removed as well. Single-end data are not affected. Output file names are unchanged; the step is streamed and leaves no intermediate files.
* **CNA correction fixed**: copy numbers are rounded to integers (round(ploidy x 2<sup>log2</sup>)), gains and losses beyond `CNA_threshold` are both corrected (CN = 0 is not corrected), the signal is divided by CN/ploidy and `CNA_threshold` is no longer truncated to an integer. Applies to both the RPGC and the GC-corrected CNA-corrected bigWigs. New config parameter `copy_number_variation > ploidy` (default 2).

<br>

#### [v0.1.2](https://github.com/sebastian-gregoricchio/SPACCa/releases/tag/0.1.2) - July 5<sup>th</sup> 2026
* The analyses snakefile now allos for samples without input/control sample.

<br>

#### [v0.1.1](https://github.com/sebastian-gregoricchio/SPACCa/releases/tag/0.1.1) - July 5<sup>th</sup> 2026
* Minor bug fixing for handling sample names in peak calling rule.
* Manual update.

<br>

#### [v0.1.0](https://github.com/sebastian-gregoricchio/SPACCa/releases/tag/0.1.0) - April 30<sup>th</sup> 2025
* First release (some manuals are still missing but the workflow is functional).
