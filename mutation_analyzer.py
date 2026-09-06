# mutation_analyzer.py

# Takes filtered driver-gene mutations and analyzes them:
#   - counts mutations per gene
#   - counts mutation TYPES (Missense, Nonsense, Frame_Shift, Silent...)
#   - tracks which patients had mutations in which gene

SEVERITY_GROUPS = {
    'Frame_Shift_Del':   'high',
    'Frame_Shift_Ins':   'high',
    'Nonsense_Mutation': 'high',
    'Translation_Start_Site': 'high',
    'Nonstop_Mutation':  'high',
    'Missense_Mutation': 'moderate',
    'In_Frame_Del':      'moderate',
    'In_Frame_Ins':      'moderate',
    'Splice_Site':       'moderate',
    'Silent':            'low',
}


def get_severity(variant_classification):
    """
    Maps a Variant_Classification value to 'high', 'moderate', 'low',
    or 'other' if it's a classification we haven't categorized.
    """
    return SEVERITY_GROUPS.get(variant_classification, 'other')


def analyze_mutations(mutations):
    """
    Takes an iterable of driver-gene mutation dicts (output of
    cancer_gene_filter.filter_driver_mutations).

    Returns a stats dict:
    {
        'total_mutations': int,
        'gene_counts': {gene: count},                    # ← get_driver_stats() uses this
        'type_counts': {variant_classification: count},
        'severity_counts': {'high': n, 'moderate': n, ...},
        'gene_type_breakdown': {gene: {classification: count}},
        'gene_patients': {gene: set_of_patient_barcodes},  # ← used by frequency_calculator
        'all_patients': set_of_all_patient_barcodes,
    }
    """
    stats = {
        'total_mutations': 0,
        'gene_counts': {},
        'type_counts': {},
        'severity_counts': {},
        'gene_type_breakdown': {},
        'gene_patients': {},
        'all_patients': set(),
    }

    for mut in mutations:
        gene = mut.get('Hugo_Symbol')
        vclass = mut.get('Variant_Classification')
        patient = mut.get('Tumor_Sample_Barcode')

        stats['total_mutations'] += 1

        # --- count mutations per gene ---
        stats['gene_counts'][gene] = stats['gene_counts'].get(gene, 0) + 1

        # --- count mutation types overall ---
        stats['type_counts'][vclass] = stats['type_counts'].get(vclass, 0) + 1

        # --- count by severity group ---
        severity = get_severity(vclass)
        stats['severity_counts'][severity] = \
            stats['severity_counts'].get(severity, 0) + 1

        # --- per-gene breakdown of mutation types ---
        if gene not in stats['gene_type_breakdown']:
            stats['gene_type_breakdown'][gene] = {}
        gene_types = stats['gene_type_breakdown'][gene]
        gene_types[vclass] = gene_types.get(vclass, 0) + 1

        # --- track WHICH PATIENTS have a mutation in this gene ---
        if patient:
            if gene not in stats['gene_patients']:
                stats['gene_patients'][gene] = set()
            stats['gene_patients'][gene].add(patient)
            stats['all_patients'].add(patient)

    return stats


def summarize_gene(stats, gene):
    """
    Small helper — pulls together everything known about ONE gene
    from the stats dict, for readable reporting.
    """
    return {
        'gene': gene,
        'mutation_count': stats['gene_counts'].get(gene, 0),
        'patient_count': len(stats['gene_patients'].get(gene, set())),
        'mutation_types': stats['gene_type_breakdown'].get(gene, {}),
    }