# CardioClassifier v2

A Python reimplementation of CardioClassifier, a disease-specific ACMG/AMP variant
classification tool for inherited cardiac conditions. The original (Whiffin et al.
2018, *Genetics in Medicine*) was a Perl/PHP tool built around a static, locally
cached snapshot of gene lists, population frequencies, and clinical databases. This
version queries the same kinds of resources live (Ensembl VEP, gnomAD, ClinVar,
ClinGen's CSpec registry) and scores evidence with the Tavtigian Bayesian point
system instead of the original's categorical ACMG rule tables.

## What it does

Given an HGVS coding-sequence variant and a diagnosis, it:

- Annotates the variant via Ensembl VEP (consequence, canonical transcript, REVEL
  and SpliceAI scores)
- Pulls population frequency from gnomAD (joint genome + exome, FAF95)
- Checks ClinVar for same-residue precedent (PS1/PM5)
- Applies gene- and disease-specific ACMG/AMP thresholds from ClinGen's CSpec
  registry where one exists, otherwise falls back to generic 2015 defaults
- Asks for a handful of case-level facts that can't be looked up automatically
  (de novo status, segregation, phasing, phenotype specificity)
- Returns a full evidence breakdown, a combined Tavtigian score, and a posterior
  probability of pathogenicity, with the VUS tier further split into VUS-low,
  VUS-mid, and VUS-high

## Files

| File | What it does |
|---|---|
| `annotationn.ipynb` | Variant annotation (VEP + gnomAD) |
| `classifierr.ipynb` | The 28 ACMG/AMP evidence rules and the Bayesian scorer |
| `mainn.ipynb` | Interactive driver, run this one |

`mainn.ipynb` loads the other two automatically if they aren't already in memory,
you don't need to run them separately first.

## Setup

```
pip install -r requirements.txt
```

Developed and run on Python 3.9. Then open `mainn.ipynb` in Jupyter and run all
cells. It will ask for a diagnosis, an HGVS variant, and a few case-level questions,
then print the full classification.

## Data files

`classifierr.ipynb` reads these local reference files, all included in this repo:

| File | Where it comes from | What's in it | Role in the pipeline |
|---|---|---|---|
| `Cardiac_G2P_cleaned_HCM_syndromic.csv` | CardiacG2P (Josephs et al. 2023), public resource | Gene-disease validity, inheritance pattern, allelic requirement, disease mechanism, variant classes | Confirms whether a (gene, diagnosis) pair is valid before running any rules; gates PVS1 (loss-of-function mechanism check), BS2/BP2 (penetrance), and PM3 (biallelic/recessive check) |
| `hotspot_regions.csv` | Walsh et al. 2019 / CSpec, published paper | Mutational hotspot coordinates and gene-wide low-benign-missense flags | Feeds PM1 (position-specific hotspot) and PP2 (gene-wide missense constraint) |
| `acmg_curations.csv` | Legacy database, literature-curated (PMID-cited) | Per-variant curated evidence: rule code, strength, evidence text, source | Supplies curator-only and functional evidence codes (PS3/BS3, segregation, etc.) for variants that have already been curated, so the user isn't asked to re-derive evidence that's already published |
| `clingen_curations.csv` | Legacy database, ClinGen's own prior classifications | Gene, cDNA, ClinGen's classification, rule code, strength | Same purpose as above, specifically for evidence drawn from ClinGen's own curation records |
| `report_comment_functional.csv` | Legacy database, literature-curated (every row PMID-cited) | Gene, variant, functional-study outcome, free-text comment | Supplies PS3/BS3 (functional study support or refute) when a variant matches a previously curated functional study |
| `original_case_control_studies.csv` | Legacy database, named cohort studies (e.g. LMM) | Gene, cDNA, disease, study source, case count, sample count, classification | Supplementary case-level evidence source, matched by (gene, cDNA), used alongside the JUL files rather than replacing them |
| `JUL_HCMgenes.tsv` | HCM case cohort, aggregate allele counts | Genomic position, ref/alt allele, AC (allele count), AN (allele number) | The "case" side of PS4's odds-ratio calculation (case cohort vs. gnomAD controls) for HCM diagnoses |
| `JUL_DCMgenes.tsv` | DCM case cohort, aggregate allele counts | Same structure as above | Same role as above, for DCM diagnoses |
| `titin_exon.csv` | Legacy database, TTN exon structure | Exon coordinates, PSI (percent spliced in) value, domain/phase info | TTN-specific PVS1 gate: a TTN truncating variant only counts toward PVS1 if it falls in a constitutively-expressed exon (PSI ≥ 90%), since TTN is too large and alternatively spliced for generic PVS1 to apply safely |
| `titin_transcript.csv` | Legacy database, TTN transcript definitions | Transcript name, Ensembl/RefSeq/UniProt IDs, protein length | Identifies which transcript's exon numbering and PSI values are relevant, matched against VEP's chosen canonical TTN transcript |
| `titin_exon_transcript.csv` | Legacy database, exon-transcript join table | Links each exon to the transcript(s) it belongs to | Resolves which exon number in a given transcript a variant's position falls in, needed to look up that exon's PSI value |

The three TTN files work together as one unit: `titin_transcript.csv` picks the
transcript, `titin_exon_transcript.csv` maps that transcript's exons, and
`titin_exon.csv` gives each exon's PSI value, that's the whole chain the PVS1
gate walks for a TTN truncating variant.

Everything else (VEP, gnomAD, ClinVar, ClinGen CSpec) is queried live, so you need
an internet connection to run it.

## Limitations

- No caching of API responses. A slow connection means a slow run, and results can
  shift over time as the underlying databases change, that's the point, not a bug.
- Curator-only evidence (segregation, de novo status, phasing, phenotype
  specificity, healthy-adult observation) has to be supplied by the user. It can't
  be inferred from public data.
- Gene-specific ACMG thresholds only exist for genes with a published ClinGen
  CSpec, currently 8 genes on this panel. Every other gene falls back to generic
  ACMG/AMP 2015 defaults.
- Built for inherited cardiac conditions specifically. Adapting it to another
  disease area means swapping the gene-disease validity source and re-checking
  every disease-specific rule and threshold.

## Background

This reimplements CardioClassifier (Whiffin et al. 2018, *Genetics in Medicine*
20(10):1246-54). Full methodology, validation against the original MYH7
benchmark, and a broader ClinVar-based evaluation are described in the
accompanying thesis.

## License

The code in this repository (`annotationn.ipynb`, `classifierr.ipynb`,
`mainn.ipynb`) is released under the MIT License, see `LICENSE`.

The reference data files in the table above are **not** covered by that licence.
Each is redistributed here for reproducibility under its original source's terms,
with attribution as listed in the Data files section. CardiacG2P and the Walsh et
al. 2019 hotspot data are public resources; the legacy-database and case-cohort
files (`acmg_curations.csv`, `clingen_curations.csv`,
`report_comment_functional.csv`, `original_case_control_studies.csv`,
`JUL_HCMgenes.tsv`, `JUL_DCMgenes.tsv`, the TTN files) are aggregate,
non-identifiable summaries only. If you reuse any of them, cite the original
source rather than this repository.
