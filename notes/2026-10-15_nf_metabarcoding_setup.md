# 2026-10-15 — nf-metabarcoding on the chimera sweep

## Aim

Measure how chimeras affect diversity estimates at the end of the complete
nf-metabarcoding pipeline. `-chimera_perc` (0 to 9 %) is the only parameter that
differs between samples. No other artefact is simulated: no sequencing errors, no
homopolymer errors. Every non-chimeric read is identical to its amplicon, so any OTU
beyond the 20 simulated taxa comes from a chimera.

## Sweep design: 10 chimera levels × 5 replicates = 50 samples

- **One community for all samples.** Grinder draws it once (`config/grinder_community.txt`:
  20 taxa, power law 0.1, seed 42) among the SILVA references carrying both primers. It is
  saved in `results/chimera_sweep/community/` (ranks, abundance file, the 20 reference
  sequences) and passed to every sample with `-abundance_file`.
- **One seed per sample:** `-random_seed` = 100 × level + replicate (1 … 5, 101 … 105,
  …, 901 … 905). Replicates of a level are independent draws of the same community.
- **One pipeline run per chimera level** (`chimera_sweep/chimera_00` …
  `chimera_sweep/chimera_09`, 5 samples each). Part B of the pipeline pools the samples
  of a run (see below), so pooling levels would mix their chimeras. Each run instead
  behaves like a study at a single chimera rate.
- Samples are named `chimera_<level>_r<replicate>`, for example `chimera_05_r3`.
  `results/chimera_sweep/summary.tsv` lists the reads, chimeric reads and distinct
  sequences of every sample.
- Each sample reads only the 20 community references instead of the whole of SILVA
  (seconds instead of ~140 s per sample). With an abundance file, grinder skips every
  reference missing from it anyway; checked: `chimera_09_r3` is byte-identical either way.
- **Reproducibility needs a fixed Perl hash seed.** Grinder builds chimeras by walking
  Perl hashes, and Perl shuffles their order in every process. A first run without a
  fixed hash seed gave identical non-chimeric reads but different chimeras for the
  same `-random_seed`. The script now sets `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`, and
  each sample is byte-identical when regenerated (checked on `chimera_09_r3`).

### Generated sweep (2026-10-15, 3 min 20 s, 77 MB)

All 50 samples passed the script's checks: 2,000 reads each, the 20 community
references, and exactly 20 distinct non-chimeric sequences (one per amplicon, no
errors). Grinder decides per read whether it becomes chimeric, so realized rates
fall a little below nominal at the high levels. Use `chimeric_pct` from `summary.tsv`
rather than the nominal level. Chimeric reads have 2 parents (4,044 reads), 3 (448)
or 4 (12).

| Level (%) | Chimeric reads, mean % (range of 5) | Distinct chimeras per sample (mean) | Distinct chimeras in the run | … with total abundance ≥ 2 | … in ≥ 2 replicates |
|---|---|---|---|---|---|
| 0 | 0.0 | 0 | 0 | 0 | 0 |
| 1 | 1.2 (0.9–1.4) | 22 | 87 | 17 | 15 |
| 2 | 2.1 (1.6–2.6) | 35 | 130 | 38 | 30 |
| 3 | 2.8 (2.6–2.9) | 48 | 184 | 37 | 34 |
| 4 | 4.0 (3.6–4.5) | 63 | 226 | 64 | 53 |
| 5 | 5.2 (4.7–5.9) | 81 | 281 | 84 | 73 |
| 6 | 6.2 (5.0–7.0) | 95 | 318 | 102 | 94 |
| 7 | 6.9 (6.0–7.5) | 102 | 319 | 115 | 109 |
| 8 | 8.2 (7.3–9.2) | 115 | 378 | 115 | 108 |
| 9 | 8.6 (7.5–9.2) | 126 | 420 | 140 | 125 |

"Total abundance ≥ 2" counts the chimeras `uchime_denovo` will see in a run
(`chimera_minsize` 2). The others are singletons in the whole run, which the abundance
and spread filters remove whatever chimera detection does.

## What was set up

Pipeline: [nf-metabarcoding](https://github.com/frederic-mahe/nf-metabarcoding)
(Frédéric Mahé), release **v0.6.0**, pinned in `scripts/05_run_nf_metabarcoding.sh`.
Nextflow pulls it from GitHub, so it is not copied into this project. It needs a Unix
machine (another laptop, or the analysis server `mgheno@194.254.138.192`), because it
does not run on the Windows workstation.

| Step | Script | Where | Output |
|------|--------|-------|--------|
| simulate amplicons | `scripts/02_run_grinder_chimera_sweep.sh` | WSL (grinder 0.5.4) | `results/chimera_sweep/` |
| reference for taxonomy (once) | `scripts/04_build_stampa_reference.sh` | nextflow machine | `results/reference/..._16S_V4_stampa.fasta` |
| stage reads, per level | `scripts/03_prepare_nf_inputs.sh chimera_sweep/chimera_05` | nextflow machine | `results/nf_input/chimera_sweep/chimera_05/` + `samplesheet.csv` |
| run pipeline, per level | `scripts/05_run_nf_metabarcoding.sh chimera_sweep/chimera_05` | nextflow machine | `results/nf_metabarcoding/chimera_sweep/chimera_05/` |
| (server only) copy / fetch | `scripts/sync_server.sh push` / `pull` | workstation | `~/metaB_sim/` on the server |

Parameters: `config/grinder_community.txt` (community), `config/grinder_params.txt`
(reads), `config/primers_16S_V4.fasta` (primers), `config/nf_metabarcoding_params.yml`
(analysis, pipeline defaults written out) and `config/nf_metabarcoding_server.config`
(resource cap on a shared machine).

## Marker: 16S rRNA V4, 515F-Y / 806RB

- 515F-Y `GTGYCAGCMGCCGCGGTAA` (Parada et al. 2016) and 806RB `GGACTACNVGGGTWTCTAAT`
  (Apprill et al. 2015): the Earth Microbiome Project pair for bacteria and archaea.
- In SILVA, 525,178 of 905,628 sequences (58 %) match both primers exactly. The
  region between the primers is 253 nt (median), 250–257 nt for 98 % of them
  (min 12, max 1,643), so amplicons with primers are about 292 nt. Only 27 amplicons
  exceed 500 nt: grinder would truncate those reads (`-read_dist 500`) and the
  pipeline would drop them, much as real overlong amplicons are lost.
- Checked on a one-sample grinder test (`-chimera_perc 9`, 143 s): all 20 reads are
  whole amplicons (292–294 nt) starting with 515F-Y and ending with the reverse
  complement of 806RB. Once their primers are removed, the 5 distinct non-chimeric
  reads are identical to their reference's region in the stampa reference.
- Chloroplasts and mitochondria also match (SILVA files them under Bacteria), as in
  real 16S data. SILVA is not filtered.

## Decisions

- **One primer file for simulation and analysis.** `config/grinder_params.txt`
  points grinder to `config/primers_16S_V4.fasta` (`-forward_reverse`). Scripts 04
  and 05 read the primers from that same file. The pipeline therefore always trims the
  primers grinder used, and its reference covers the same region.
- **Grinder in amplicon mode.** Grinder keeps only the references that match both
  primers exactly. IUPAC codes are allowed, and U matches T: BioPerl's
  `AmpliconSearch` builds RNA-aware regexes, so SILVA needs no DNA conversion for grinder.
  Each amplicon runs from primer to primer, primers included.
- **Whole-amplicon reads (`-read_dist 500`, `-unidirectional 1`).** The pipeline's
  cutadapt step needs both primers in each read and discards reads missing either one.
  Grinder starts amplicon reads at the forward primer and shortens a read to its amplicon
  when the read is longer (`Grinder.pm`, `next_single_read`). A read length above the
  longest amplicon therefore gives whole amplicons, which play the part of merged read
  pairs. Reads stay single-end, so the pipeline skips its merging step.
- **Qualities `-qual_levels 40 12`.** Correct bases get Q40: the top quality level of
  NovaSeq X (it reports 12, 24 and 40), and close to the quality of merged pairs (vsearch
  caps it at 41). Erroneous bases would get Q12, the lowest level, but none are
  simulated. The first value must stay ≥ 37. Below that, every sequence exceeds the
  pipeline's `max_ee` (0.0002 expected errors per base) and the table comes out empty.
- **Samplesheet (`--input`) rather than `--fastq_folder`.** Sample IDs are explicit
  (`chimera_00` ... `chimera_09`) instead of being guessed by the pipeline's R1/R2
  file-name patterns, and `-resume` notices changed inputs (the pipeline README
  recommends it for reproducible runs).
- **U → T when staging.** Grinder writes RNA reads, because SILVA is RNA. Script 03
  converts the sequence lines to DNA and keeps headers (grinder's `reference=`,
  `amplicon=`, `description=`) and qualities unchanged.
- **Reference cut to the amplified region.** Stampa scores hits with vsearch
  `--iddef 1`, which counts terminal gaps. A 253-nt read aligned to a ~1,450-nt
  full-length SILVA sequence would score below `stampa_id` 0.5 and get no assignment.
  Script 04 therefore keeps only the region between the primers (primers excluded,
  as for the trimmed reads), found with the same exact IUPAC matching as grinder. Every
  simulated reference is thus in the taxonomic reference.

## To keep in mind when reading results

- **Most chimeras are seen once.** Grinder test with the current design (20 taxa,
  2,000 reads, `-chimera_perc 9`, 144 s per sample): 201 chimeric reads (182 bimeras,
  18 trimeras, 1 quadrimera), forming 145 distinct chimeras, 118 of which occur once.
  `uchime_denovo` only sees sequences with a total abundance ≥ 2 (`chimera_minsize`).
  The occurrence table then keeps a sequence only if its abundance is ≥ 3 or it occurs
  in ≥ 2 samples. A chimera seen once in the whole run is therefore dropped by the
  filters, not by chimera detection.
- **Samples are pooled in Part B.** Global dereplication, clustering, chimera detection
  and the abundance/spread filter run on all samples of a run together, so what happens
  to a chimera in one sample depends on the other samples of the run. In an earlier
  test, 33 of the 86 distinct chimeras of a 5 % sample also occurred in a 9 % sample.
  A chimera seen once in each has a spread of 2: it is tested by uchime and, if missed,
  kept as an OTU in both samples. This is why each chimera level is its own run: the 5
  replicates of a level share chimeras the way samples of one real study would.
- **Chimeras recur across replicates.** With 20 parents and breakpoints limited to
  shared 10-mers (`-chimera_kmer 10`), the same chimera can arise in several
  replicates even with independent seeds, as real PCR chimeras do.
- **What the final table contains.** The pipeline publishes only the filtered table
  (`occurrence_table/chimera_05_5_samples_table.tsv` for level 5, plus its
  taxonomy-annotated copy):
  - OTUs flagged by uchime (`Y`) are already removed, and the `chimera` column is
    reset to `N` after mumu.
  - mumu can also merge a chimera into a parent it co-occurs with (≥ 84 % similar).
  - Extra OTUs in the final table are chimeras that escaped both steps.
  - The pre-filter tables stay in the work directory (`results/nf_work/<run>/`).
    The uchime results (`.uchime`) are published.
- **Qualities are uniform.** With no errors simulated, every base is Q40, so the
  expected-error filter (`max_ee`) removes nothing (ee/length = 0.0001 ≤ 0.0002).

## History

- First test run (same day) used grinder's shotgun mode: random 100-nt windows without
  primers. The pipeline would have needed `no_trimming`, would have clustered windows
  rather than taxa, and could not have detected chimeras between windows. Replaced by
  amplicon mode before any pipeline run. Its outputs were moved to
  `results/chimera_sweep_shotgun_2026-10-15/`.
- The first amplicon design (10 samples, one per level, all with seed 42, one pooled
  run) was replaced by 10 levels × 5 replicates, one run per level, before any
  pipeline run.
