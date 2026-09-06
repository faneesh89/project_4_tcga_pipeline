# TCGA Cancer Driver Gene Mutation Analysis Pipeline

A Python pipeline that analyzes real somatic mutation data from The Cancer
Genome Atlas (TCGA) to identify how frequently known cancer driver genes
are mutated across a patient cohort.

Built on real, publicly available TCGA-LUAD (Lung Adenocarcinoma) data
downloaded from the NCI Genomic Data Commons (GDC) — not synthetic or
toy data.

---

## What it does

Given a folder of patient MAF (Mutation Annotation Format) files, the
pipeline:

1. Parses all mutation records across every patient file
2. Filters to a curated set of 20 well-documented cancer driver genes
3. Classifies mutations by type (missense, nonsense, frameshift, silent, etc.)
4. Calculates what percentage of the *full patient cohort* has each gene
   mutated — correctly counting each patient once, even if they carry
   multiple mutations in the same gene
5. Generates an HTML report summarizing the results

## Results on the included dataset

Run against 30 real TCGA-LUAD patients (8,705 total somatic mutations):

| Gene | Patients Mutated | Frequency | Mutation Types |
|------|-------------------|-----------|-----------------|
| TP53  | 13/30 | 43.3% | Missense, Frame Shift Del |
| KRAS  | 11/30 | 36.7% | Missense |
| EGFR  | 7/30  | 23.3% | Missense, In Frame Del |
| STK11 | 3/30  | 10.0% | Nonsense, Missense |
| KEAP1 | 3/30  | 10.0% | Missense, Frame Shift Del |

These results are consistent with published mutation frequency profiles
for lung adenocarcinoma — TP53 and KRAS as the dominant drivers, with
EGFR, STK11, and KEAP1 appearing at frequencies matching the literature.

---

## Pipeline architecture

```
data/*.maf.gz
      │
      ▼
maf_parser.py            → reads compressed MAF files directly (no manual
                             extraction), maps ~140 raw columns down to the
                             9 fields the pipeline actually needs, via a
                             name→position lookup built from each file's
                             own header (not a fixed column index)
      │
      ▼
cancer_gene_filter.py    → keeps only mutations hitting one of 20 curated
                             driver genes; tags each with its role
                             (oncogene / tumor suppressor)
      │
      ▼
mutation_analyzer.py     → counts mutations per gene, classifies mutation
                             types, and tracks which *unique patients*
                             have each gene mutated
      │
      ▼
frequency_calculator.py  → converts patient counts into cohort-wide
                             percentages
      │
      ▼
report_generator.py      → writes the final HTML report
```

Each stage is a separate file with a single responsibility. `main.py`
wires them together and can be re-run against any folder of MAF files
with the same structure.

---

## Data source

- **Program:** The Cancer Genome Atlas (TCGA), accessed via the NCI
  Genomic Data Commons (GDC) Data Portal
- **Project:** TCGA-LUAD (Lung Adenocarcinoma)
- **File type:** Masked Somatic Mutation MAF (`aliquot_ensemble_masked`)
  — open-access files that have already been through GDC's ensemble
  variant-calling aggregation (results cross-checked across MuTect2,
  VarScan2, MuSE, and Pindel) and masking (germline variants and
  low-quality calls removed)
- 30 patient files were used for this analysis (a manageable, real subset
  rather than the full multi-thousand-case project)

Somatic MAFs are publicly distributable under GDC's data access
policies; no patient-identifying germline information is included.

---

## Key design decisions

**Curated 20-gene list instead of the full COSMIC Cancer Gene Census.**
COSMIC's Cancer Gene Census contains 719 genes, most with very low
mutation frequency and limited signal for a project at this scope.
Gene selection here prioritizes genes with strong, documented mutation
frequency evidence from TCGA PanCancer Atlas analyses (e.g. TP53 ~35%,
KRAS ~11%, PIK3CA ~13% across cancer types generally), plus several
LUAD-specific drivers (STK11, KEAP1). Full COSMIC access requires a
license/registration not available for this project; a documented
alternative would be deriving the gene list programmatically from a
licensed COSMIC export with a defined frequency threshold — noted here
as a possible future enhancement rather than built in, to keep the
pipeline's dependencies self-contained.

**Only 9 of ~140 MAF columns are used.** GDC's Masked Somatic Mutation
files carry annotations from multiple variant callers and population
databases layered on top of the core mutation record. Only the fields
needed to identify the gene, mutation type, location, and patient are
extracted (`Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`,
`Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`,
`Variant_Type`, `Tumor_Sample_Barcode`).

**Generator-based file parsing.** `maf_parser.py` uses `yield` rather
than building a full list in memory, so the pipeline can scale to much
larger MAF files without a memory bottleneck.

---

## A bug I found and fixed

The first working version of this pipeline calculated frequency using
only the patients who had at least one driver-gene mutation as the
denominator. Five of the 30 patients in this dataset had mutations
elsewhere in the genome but none in the 20 tracked driver genes — a
biologically real and expected outcome.

Because those 5 patients were still legitimately part of the analyzed
cohort, excluding them from the denominator inflated every reported
frequency by roughly 9 percentage points (e.g. TP53 initially showed
52.0% instead of the correct 43.3%).

The fix: `main.py` now counts the true total patient count from the
*unfiltered* mutation data, before any gene filtering happens, and
passes that number explicitly into the frequency calculation — rather
than relying on a count that could only ever include patients who
happened to survive the filter.

```python
all_patient_barcodes = set(
    m['Tumor_Sample_Barcode'] for m in all_mutations
    if m.get('Tumor_Sample_Barcode')
)
total_patients_analyzed = len(all_patient_barcodes)
```

## Handling patients with multiple mutations

A patient can carry more than one mutation in the same gene (biologically
expected for tumor suppressors, which often require both gene copies to
be damaged). `mutation_analyzer.py` tracks patients per gene using a
`set()`, not a counter, so a patient with two TP53 mutations is still
counted once toward TP53's patient frequency — while still being
reflected in the separate total mutation count.

---

## Known limitations

- Structural variants (large deletions, gene fusions, copy number
  changes) are not captured — MAF files record point mutations and
  small indels only.
- Severity classification (`high` / `moderate` / `low`) is a
  simplification based on mutation type alone. Real functional impact
  depends heavily on *where* in the protein a mutation lands and the
  gene's role (oncogene vs. tumor suppressor) — not modeled here.
- The 20-gene list is a curated subset, not a data-derived or
  exhaustive driver gene catalog.
- Frequencies reflect this specific 30-patient cohort only and should
  not be generalized as population-wide statistics.

## Possible future enhancements

- Statistical significance testing (e.g. comparing observed mutation
  rate against a background rate) rather than reporting raw frequency
  alone
- Programmatic gene list derived from a licensed COSMIC Cancer Gene
  Census export
- Docker containerization for environment reproducibility

---

## Running the pipeline

```bash
git clone <this-repo>
cd project_4_tcga_pipeline
# place .maf.gz files in data/
python3 main.py
```

Output is written to `output/cancer_report.html`.

## Requirements

- Python 3.x (standard library only — no external dependencies)

---

## Project context

This is the fourth in a series of four genomics pipeline projects,
each building on real public data:

1. FASTQ quality control clone
2. VCF variant analysis (ClinVar, 4.4M variants)
3. Pharmacogenomics variant interpretation (CPIC/PharmVar)
4. **TCGA cancer driver gene mutation analysis** (this project)