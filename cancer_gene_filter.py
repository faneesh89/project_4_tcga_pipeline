# cancer_gene_filter.py

CANCER_GENES = {
    # ---- Tumor Suppressors ----
    'TP53':    'tumor_suppressor',
    'RB1':     'tumor_suppressor',
    'BRCA1':   'tumor_suppressor',
    'BRCA2':   'tumor_suppressor',
    'APC':     'tumor_suppressor',
    'PTEN':    'tumor_suppressor',
    'ARID1A':  'tumor_suppressor',
    'ATM':     'tumor_suppressor',
    'SMARCA4': 'tumor_suppressor',
    'STK11':   'tumor_suppressor',
    'KEAP1':   'tumor_suppressor',
    'NF1':     'tumor_suppressor',
    'CDKN2A':  'tumor_suppressor',

    # ---- Oncogenes ----
    'KRAS':    'oncogene',
    'PIK3CA':  'oncogene',
    'BRAF':    'oncogene',
    'EGFR':    'oncogene',
    'NRAS':    'oncogene',
    'KMT2C':   'oncogene',
    'KMT2D':   'oncogene',
}


def is_driver_gene(gene):
    """Check if a gene is a known cancer driver."""
    return gene in CANCER_GENES


def get_gene_type(gene):
    """Return 'tumor_suppressor', 'oncogene', or None for a given gene."""
    return CANCER_GENES.get(gene, None)


def filter_driver_mutations(mutations):
    """
    Filter a list of mutations to keep only those in driver genes.
    """
    drivers = []
    for mut in mutations:
        gene = mut.get('Hugo_Symbol')
        if gene in CANCER_GENES:
            mut['gene_type'] = CANCER_GENES[gene]
            drivers.append(mut)
    return drivers


def get_driver_stats(stats):
    """
    Extract driver gene statistics from mutation counts.
    """
    drivers = {}
    for gene, count in stats['gene_counts'].items():
        if gene in CANCER_GENES:
            drivers[gene] = {
                'count': count,
                'type': CANCER_GENES[gene]
            }
    return drivers