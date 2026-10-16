# metaB sim

```
metaB_sim/
├── README.md              # what the project is, how to run it (1 paragraph to start)
├── LICENSE
├── .gitignore             # keeps big files out of git
├── .gitattributes         # LF line endings for scripts/configs (they run on the Linux server)
│
├── data/
│   ├── reference/         # inputs you didn't create
│   │   └── SILVA_144_SSURef_NR99_tax_silva_trunc.fasta
│   └── (no other files; keep this read-only in spirit)
│
├── scripts/               # all code lives here
│   ├── 01_make_community.py              # builds the abundance ranks (ranks file)
│   ├── 02_run_grinder_chimera_sweep.sh   # runs grinder (0-9 % chimeras x 5 replicates) + checks
│   ├── 03_prepare_nf_inputs.sh           # Grinder reads -> nf-metabarcoding inputs (U->T, samplesheet)
│   ├── 04_build_stampa_reference.sh      # SILVA -> amplified region, stampa format (taxonomy)
│   ├── 05_run_nf_metabarcoding.sh        # runs nf-metabarcoding (on the server)
│   └── sync_server.sh                    # push project to the server / pull results back
│
├── config/                # parameters, so you never hard-code numbers in scripts
│   ├── grinder_community.txt             # simulated community (richness, abundance model), drawn once
│   ├── grinder_params.txt                # per-sample reads: depth, amplicon mode, qualities, errors
│   ├── primers_16S_V4.fasta              # 515F-Y / 806RB, used by grinder and the pipeline
│   ├── nf_metabarcoding_params.yml       # nf-metabarcoding analysis parameters
│   └── nf_metabarcoding_server.config    # nf-metabarcoding resources on the server
│
├── results/               # everything the pipeline generates (can be deleted and rebuilt)
│   ├── community/
│   │   └── grinder-ranks.txt  # 44 MB, don't commit
│   ├── simulated/
│   │   └── grinder-reads.fa
│   ├── chimera_sweep/         # Grinder output: community/, summary.tsv, and one folder
│   │                          #   per chimera level (chimera_00 ... chimera_09, 5 samples each)
│   ├── nf_input/<run>/        # staged reads + samplesheet
│   ├── reference/             # SILVA cut to the amplified region, stampa format
│   ├── nf_metabarcoding/<run>/  # pipeline results (occurrence tables, logs, reports)
│   └── nf_work/<run>/         # nextflow work dir (needed by -resume, delete when done)
│
└── notes/                 # your lab notebook: decisions, parameters tried, questions
    ├── 2026-10-14_first_run.md
    └── 2026-10-15_nf_metabarcoding_setup.md
```

## Running nf-metabarcoding

Grinder simulates 16S rRNA V4 amplicons (primers 515F-Y / 806RB,
`config/primers_16S_V4.fasta`): 10 chimera levels (0 to 9 %) × 5 replicates, all
from the same 20-taxon community. The reads are then processed with
[nf-metabarcoding](https://github.com/frederic-mahe/nf-metabarcoding) v0.6.0, one
pipeline run per chimera level, on a Unix machine with `nextflow` (Java ≥ 17) and
either conda or a container engine (nextflow does not run on Windows). Setup
decisions and known limits are in `notes/2026-10-15_nf_metabarcoding_setup.md`.

Simulate the reads (grinder is installed in WSL Ubuntu), from the project root:

```bash
wsl -d Ubuntu -- bash scripts/02_run_grinder_chimera_sweep.sh
```

On the machine running nextflow, from the project root, inside `tmux`:

```bash
bash scripts/04_build_stampa_reference.sh
for run in results/chimera_sweep/chimera_??; do
    bash scripts/03_prepare_nf_inputs.sh "chimera_sweep/${run##*/}"
    bash scripts/05_run_nf_metabarcoding.sh "chimera_sweep/${run##*/}"
done
```

To use the analysis server instead, `bash scripts/sync_server.sh push chimera_sweep`
copies the project there and `bash scripts/sync_server.sh pull chimera_sweep` fetches
the results back.
