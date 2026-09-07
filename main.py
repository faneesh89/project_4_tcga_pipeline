# main.py
# ─────────────────────────────────────────────────────────────────────────────
# Runs the whole Project 4 pipeline end to end:
#
#   1. find all .maf.gz files in data/
#   2. parse them into mutation dicts       (maf_parser)
#   3. count TRUE total patients (before filtering — see note below)
#   4. keep only cancer driver genes        (cancer_gene_filter)
#   5. count mutations, types, patients     (mutation_analyzer)
#   6. calculate per-gene frequencies       (frequency_calculator)
#   7. write the HTML report                (report_generator)
#
# NOTE on the patient denominator:
# Some patients have mutations but NONE in our driver gene list. They still
# count as "analyzed" — so the frequency denominator must be the total
# patients in the raw data, not just those with driver mutations. Otherwise
# every percentage comes out inflated.
#
# Run with:  python3 main.py
# ─────────────────────────────────────────────────────────────────────────────

import glob
import sys

from maf_parser import read_multiple_mafs
from cancer_gene_filter import filter_driver_mutations
from mutation_analyzer import analyze_mutations
from frequency_calculator import calculate_frequencies, get_summary_line
from report_generator import generate_html_report


DATA_DIR = "data"
OUTPUT_FILE = "output/cancer_report.html"


def main():
    # ---- Step 1: find the input files ----
    maf_files = sorted(glob.glob(f"{DATA_DIR}/*.maf.gz"))

    if not maf_files:
        print(f"ERROR: no .maf.gz files found in '{DATA_DIR}/'")
        print("Check that your downloaded TCGA files are in that folder.")
        sys.exit(1)

    print(f"Found {len(maf_files)} MAF files in '{DATA_DIR}/'")
    print("-" * 60)

    # ---- Step 2: parse all files into mutation dicts ----
    all_mutations = list(read_multiple_mafs(maf_files))

    # ---- Step 3: count the TRUE patient total, BEFORE driver filtering ----
    all_patient_barcodes = set(
        m['Tumor_Sample_Barcode'] for m in all_mutations
        if m.get('Tumor_Sample_Barcode')
    )
    total_patients_analyzed = len(all_patient_barcodes)

    print("-" * 60)
    print(f"Total mutations parsed: {len(all_mutations)}")
    print(f"Total patients analyzed: {total_patients_analyzed}")

    # ---- Step 4: keep only mutations in known cancer driver genes ----
    driver_mutations = filter_driver_mutations(all_mutations)

    print(f"Driver-gene mutations found: {len(driver_mutations)}")

    if not driver_mutations:
        print("No mutations in the driver gene list were found.")
        print("The pipeline ran correctly — this dataset just had no hits.")
        sys.exit(0)

    # ---- Step 5: analyze — counts, types, per-patient tracking ----
    stats = analyze_mutations(driver_mutations)

    # how many patients had at least ONE driver mutation?
    patients_with_drivers = len(stats['all_patients'])
    patients_without = total_patients_analyzed - patients_with_drivers
    if patients_without > 0:
        print(f"  ({patients_without} patient(s) had no mutations in the "
              f"driver gene list — still counted in the denominator)")

    # ---- Step 6: turn counts into per-gene frequencies ----
    results = calculate_frequencies(
        stats, total_patients=total_patients_analyzed
    )

    print("-" * 60)
    print(get_summary_line(results, total_patients_analyzed))
    print("-" * 60)

    # quick console preview of the top genes
    print(f"{'Gene':<10} {'Patients':<12} {'Frequency':<12} {'Mutations'}")
    print("-" * 60)
    for r in results[:15]:
        patients = f"{r['patient_count']}/{r['total_patients']}"
        freq = f"{r['frequency_pct']}%"
        print(f"{r['gene']:<10} {patients:<12} {freq:<12} {r['mutation_count']}")

    # ---- Step 7: write the HTML report ----
    print("-" * 60)
    generate_html_report(
        results, stats,
        output_file=OUTPUT_FILE,
        total_patients=total_patients_analyzed
    )


if __name__ == "__main__":
    main()