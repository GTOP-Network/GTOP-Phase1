# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/24 11:32
@Email   : xuechao@szbl.ac.cn
@Desc    : Detect ORF fa.
"""
import gzip
import logging
import os
import re
import numpy as np
import pandas as pd
from util import GTOP_Tissue_Meta, GeneAnnotation, PROJ_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)



def __write_faa(faa_path,id_set,out_f,tid_gene=None,transcript_index=0,gene_index=1):
    if faa_path.endswith('.gz'):
        open_func = gzip.open
        open_kwargs = {'mode': 'rt', 'encoding': 'utf-8'}
    else:
        open_func = open
        open_kwargs = {'mode': 'r', 'encoding': 'utf-8'}
    k=0
    current_id=None
    current_seq = []
    with open_func(faa_path, **open_kwargs) as in_f:
        for line in in_f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id in id_set:
                    gene_name=current_gene_id
                    if tid_gene:
                        gene_name = tid_gene[current_id]
                    out_f.write(f'>{current_id} {gene_name}\n')
                    out_f.write(''.join(current_seq) + '\n')
                    k+=1
                header = line.lstrip('>')
                parse_arr=[part.strip() for part in re.split(r'\s+|\|', header) if part]
                current_id = parse_arr[transcript_index]
                current_gene_id = parse_arr[gene_index]
                current_seq = []
            else:
                current_seq.append(line.strip('*').strip())
        if current_id in id_set:
            gene_name = current_gene_id
            if tid_gene:
                gene_name = tid_gene[current_id]
            out_f.write(f'>{current_id} {gene_name}\n')
            out_f.write(''.join(current_seq) + '\n')
            k+=1
    return k


def make_enhanced_cds_fa():
    gencode_v47_cds_protein_fa_path=f'{PROJ_DIR}/raw_data/GMTiP/ref/LRS/gencode.v47.pc_translations.fa.gz'
    main_dir=f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'
    gtop_cds_protein_fa_path=f'{main_dir}/sqanti3_clean/filter.clean.faa'
    enhanced_gtf_prefix=f'{main_dir}/enhanced_gtf/GTOP_novel-GENCODE_v47'
    output_faa=f'{main_dir}/enhanced_faa/GTOP_novel-GENCODE_v47.faa'

    os.makedirs(os.path.dirname(output_faa),exist_ok=True)
    ga=GeneAnnotation().load_df(enhanced_gtf_prefix,remove_ver=False)
    tid_gene=ga.get_isoform_gene_map(gene_key='gene_id')
    tid_type=ga.get_isoform_type_map()
    print(f'load {len(tid_gene)} isoforms; novel {len([x for x in tid_gene.keys() if not x.startswith("ENS")])}.')
    id_set=set([t for t in tid_gene.keys() if tid_type[t] == 'protein_coding'])
    all_k=0
    with open(output_faa, 'w', encoding='utf-8') as out_f:
        # write gencode isoform
        logger.info(f'start merge gencode')
        k=__write_faa(gencode_v47_cds_protein_fa_path, id_set, out_f, tid_gene, transcript_index=1)
        logger.info(f'merge gencode {k}')
        all_k+=k
        # write novel isoform
        k=__write_faa(gtop_cds_protein_fa_path,id_set,out_f,tid_gene,transcript_index=0)
        logger.info(f'merge novel {k}')
        all_k+=k
        logger.info(f'done {all_k} proteins; save to {output_faa}')


def tissue_based_protein_cds_fa(min_tpm=1,min_sample=None):
    # if min_sample is None ,use median to determine.
    main_dir=f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS'
    isoform_expr_path = f'{main_dir}/quantification/flair_quant/combined/GTOP.transcript.tpm.flair.tsv'
    enhanced_faa_path=f'{main_dir}/isoform_discovery/merged/enhanced_faa/GTOP_novel-GENCODE_v47.faa'
    faa_dir=f'{main_dir}/isoform_discovery/merged/tissue_faa_tpm_0.1'
    os.makedirs(f'{faa_dir}',exist_ok=True)
    df=pd.read_csv(f'{isoform_expr_path}',index_col=0,sep='\t')
    # # exclude ambitious tissue samples
    # df=df[[x for x in df.columns if 'CA241-5032' not in x]]
    logger.info(f'load {df.shape[1]} samples')
    gtm=GTOP_Tissue_Meta()
    tissue_id={}
    for x in df.columns:
        tis_code=x.split('-')[2]
        tis=gtm.get_meta_by_tissue_code(tis_code,'Tissue')
        if tis not in tissue_id:
            tissue_id[tis]=[]
        tissue_id[tis].append(x)
    for tis,sid in tissue_id.items():
        sdf=df.loc[:,sid]
        if min_sample is None:
            median_expr=sdf.median(axis=1,skipna=True)
            target_isoforms=set(median_expr.loc[median_expr>min_tpm].index.tolist())
        else:
            target_isoforms = set(sdf.index[(sdf > min_tpm).sum(axis=1) >= min_sample])

        logger.info(f'{tis} expressed {len(target_isoforms)} protein-coding genes')
        correct_tis=re.sub(r'\s+',"_",tis)
        faa_output=f'{faa_dir}/{correct_tis}.faa'
        with open(faa_output, 'w', encoding='utf-8') as out_f:
            # write gencode isoform
            logger.info(f'start write')
            k = __write_faa(enhanced_faa_path, target_isoforms, out_f, tid_gene=None, transcript_index=0, gene_index=1)
            logger.info(f'done {len(target_isoforms)} isoforms, {k} proteins; save to {faa_output}')



if __name__ == '__main__':
    ## median tpm>5.
    tissue_based_protein_cds_fa(min_tpm=5)
    ## median tpm>0.1.
    tissue_based_protein_cds_fa(min_tpm=0.1)
