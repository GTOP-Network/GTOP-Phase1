# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/30 15:02
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import logging
import os

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyranges as pr
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


PROJ_DIR = os.environ.get("PROJECT_ROOT")
output_2tool_gtf = f'{PROJ_DIR}/project/GMTiP-RNA/20251031/long_read/HPC/output_lrw/isoform_discovery/merged/gtf_4/merged.2tools.gtf'
output_2tool_with_gene_gtf = f'{PROJ_DIR}/project/GMTiP-RNA/20251031/long_read/HPC/output_lrw/isoform_discovery/merged/gtf_4/merged.2tools.with_gene1.gtf'

def assign_gene_id_manual(input_gtf, output_gtf):
    logger.info("Reading transcript GTF...")
    gtf = pr.read_gtf(input_gtf, duplicate_attr=True)
    df = gtf.df



# ============================
# Usage example
# ============================
if __name__ == "__main__":
    input_gtf = output_2tool_gtf
    output_gtf = output_2tool_with_gene_gtf
    assign_gene_id_manual(input_gtf, output_gtf)

