# frequency_calculator.py
# ─────────────────────────────────────────────────────────────────────────────
# Turns the raw counts from mutation_analyzer.py into the actual clinically
# meaningful numbers: what PERCENTAGE of patients have each gene mutated.
#
# IMPORTANT — the denominator problem:
# mutation_analyzer.py only ever sees FILTERED driver mutations, so its
# 'all_patients' set can only contain patients who had at least one driver
# mutation. Patients who were analyzed but had NO driver-gene mutations
# would silently vanish from the denominator, inflating every percentage.
#
# So main.py counts the TRUE patient total from the unfiltered data and
# passes it in via total_patients. Only falls back to stats['all_patients']
# if no explicit total is given.
# ─────────────────────────────────────────────────────────────────────────────


def calculate_frequencies(stats, total_patients=None):
    """
    stats: the dict from mutation_analyzer.analyze_mutations()
    total_patients: the TRUE number of patients analyzed, including those
                    with zero driver-gene mutations. Pass this from main.py.
                    If omitted, falls back to counting only patients who had
                    driver mutations (which inflates percentages — see above).

    Returns a list of per-gene result dicts, sorted by frequency
    (most frequently mutated gene first).
    """
    if total_patients is None:
        # fallback — will UNDERCOUNT the denominator if any patients
        # had no driver mutations at all
        total_patients = len(stats.get('all_patients', set()))

    if total_patients == 0:
        return []

    results = []

    for gene, patient_set in stats.get('gene_patients', {}).items():
        patient_count = len(patient_set)
        frequency_pct = (patient_count / total_patients) * 100

        results.append({
            'gene': gene,
            'patient_count': patient_count,
            'total_patients': total_patients,
            'frequency_pct': round(frequency_pct, 1),
            'mutation_count': stats.get('gene_counts', {}).get(gene, 0),
            'mutation_types': stats.get('gene_type_breakdown', {}).get(gene, {}),
        })

    # most-frequently-mutated first
    results.sort(key=lambda r: r['frequency_pct'], reverse=True)

    return results


def format_mutation_types(mutation_types):
    """
    Turns {'Missense_Mutation': 12, 'Nonsense_Mutation': 4} into a
    readable string like "12 Missense, 4 Nonsense" for the report.
    """
    if not mutation_types:
        return "—"

    parts = []
    for vclass, count in sorted(mutation_types.items(),
                                 key=lambda x: x[1], reverse=True):
        short_name = vclass.replace('_Mutation', '').replace('_', ' ')
        parts.append(f"{count} {short_name}")

    return ", ".join(parts)


def get_summary_line(results, total_patients):
    """
    One-line overall summary for the top of the report.
    """
    if not results:
        return "No driver gene mutations found."

    top = results[0]
    return (f"Analyzed {total_patients} patients. "
            f"Most frequently mutated driver gene: {top['gene']} "
            f"({top['patient_count']}/{total_patients} patients, "
            f"{top['frequency_pct']}%).")