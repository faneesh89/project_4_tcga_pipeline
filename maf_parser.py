import gzip
import os

NEEDED_COLUMNS = [
    'Hugo_Symbol',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Reference_Allele',
    'Tumor_Seq_Allele2',
    'Variant_Classification',
    'Variant_Type',
    'Tumor_Sample_Barcode'
]

def read_maf(maf_file):
    if not os.path.exists(maf_file):
        raise FileNotFoundError(f'MAF file not found : {maf_file}')

    if maf_file.endswith('.gz'):
        opener = gzip.open(maf_file,'rt')
    else:
        opener = open(maf_file,'r')

    mutation_count = 0

    with opener as f:
        header_line = None
        for line in f:
            if line.startswith('#'):
                continue
            header_line = line.strip()
            break
        if header_line is None:
            print(f'WARNING : No header/data found in {maf_file}')
            return

        header_columns = header_line.split('\t')


        col_index = {}
        for name in NEEDED_COLUMNS:
            if name in header_columns:
                col_index[name] = header_columns.index(name)
            else:
                print(f"WARNING : expected column '{name}' not found" f"in {maf_file}")

        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            mutation = {}
            for name in NEEDED_COLUMNS:
                if name in col_index:
                    idx = col_index[name]
                    mutation[name] = fields[idx] if idx < len(fields) else None
                else:
                    mutation[name] = None
            mutation_count += 1
            yield mutation

    print(f"{maf_file} : {mutation_count} mutations parsed")

def read_multiple_mafs(maf_files):

    for maf_file in maf_files:
        for mutation in read_maf(maf_file):
            yield mutation
