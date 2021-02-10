#!/usr/bin/env python

#!/usr/bin/env python

from pathlib import Path

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler


def main(input_primer_info_tab: Path,
         output_primer_bed: Path):
    """Convert primer_info.tab to BED file for iVar primer trimming"""
    
    df = pd.read_table(input_primer_info_tab)
    alt_name = 'NC_045512.2'
    ref_name = 'MN908947.3'

    bed_entries = []
    for i, r in df.iterrows():
        if r.chrom == alt_name:
            ls, le = r.left_start, r.left_end
            rs, re = r.right_start, r.right_end
            bed_entries += [(ref_name, ls-1, le, f'cleanplex_{i+1}_LEFT', 'cleanplex_1', '+'),
                           (ref_name, rs-1, re, f'cleanplex_{i+1}_RIGHT', 'cleanplex_2', '-'),]

    df_out =  pd.DataFrame(bed_entries)
    df_out.to_csv(output_primer_bed, sep='\t', index=False, header=False)


if __name__ == '__main__':
    typer.run(main)
