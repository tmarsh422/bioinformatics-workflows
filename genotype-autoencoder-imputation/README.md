[README.md](https://github.com/user-attachments/files/29477756/README.md)
# Data

This directory is where the pipeline expects to find IMPUTE2 inputs.
**Nothing in `data/raw/` is committed to this repo** (see `.gitignore`)
because genotype data is either large, restricted, or both.

## Expected structure

```
data/raw/
├── Softwares/
│   └── impute_v2.3.0_Windows/
│       └── impute2.exe
├── RefPanels/
│   ├── genetic_map_chr7_combined_b37.txt
│   ├── 1000GP_Phase3_chr7.hap/
│   └── 1000GP_Phase3_chr7.legend/
└── GeneFiles/
    ├── twin_chr7.gens
    ├── twin_chr7.sample
    └── twin_chr7.strand.txt
```

## Where to get each piece

- **IMPUTE2 binary**: download from the
  [IMPUTE2 homepage](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html).
- **1000 Genomes Phase 3 reference panels**: publicly available from the
  same IMPUTE2 site under "Reference-format haplotypes." These are
  shareable/public data, so it's reasonable to include a download
  script here if you want a fully reproducible setup.
- **GeneFiles (study sample genotypes)**: this is study-specific
  genotype data and is **not included**. Substitute your own
  `.gens`/`.sample`/`.strand` files, or generate a small synthetic
  example if you want others to be able to run the pipeline end-to-end
  without real genotype data.

## Paths

All file locations are configured in `config.yml` at the project root
— update `paths.base_directory` and the relative paths under it if
your layout differs from the structure above.
