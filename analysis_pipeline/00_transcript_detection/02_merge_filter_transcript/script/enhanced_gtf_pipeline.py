# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/5 09:24
@Email   : xuechao@szbl.ac.cn
@Desc    :  
Extract novel isoforms from SQANTI3 filtered GTF based on classification report.

Enhanced version:
 - Extracts novel isoforms (Isoform + not known categories)
 - Replaces gene_id and transcript_id in GTF using associated_gene / associated_transcript mappings
"""
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
from multiprocessing.spawn import is_forking
from pathlib import Path
from tkinter.constants import CURRENT

import numpy as np
import pandas as pd

import argparse
import pandas as pd
from numpy.ma.core import outer

from util import GeneAnnotation

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.environ.get("PROJECT_ROOT")

REF_DIR=f'{PROJ_DIR}/raw_data/GMTiP/ref/LRS'
REF_GENOME_FASTA = f"{REF_DIR}/genome.fa"
REF_GTF = f"{REF_DIR}/gencode.v47.annotation.gtf"
GFFREAD_BIN = f'/lustre/home/cxue/software/gffread-0.12.7.Linux_x86_64/gffread'


class GTFParser:
    def __init__(self, gtf_file=None):
        self.gtf_file = gtf_file
        self.gtf_df:pd.DataFrame
        self.comment = None
        self.cols = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    @staticmethod
    def parse_attributes(attr_str):
        attr_dict = {}
        for attr in attr_str.strip().split(";"):
            attr = attr.strip()
            if not attr:
                continue
            if " " in attr:
                key, val = attr.split(" ", 1)
                attr_dict[key] = val.strip('"')
        return attr_dict

    @staticmethod
    def build_attributes(attr_dict):
        return "; ".join(f'{k} "{v}"' for k, v in attr_dict.items()) + ";"

    def load_comment(self):
        comment_arr=[]
        with open(self.gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    comment_arr.append(line.strip())
                else:
                    break
        return comment_arr

    def create_from_arr(self,arr):
        if isinstance(arr, (list)):
            cols = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            df = pd.DataFrame(arr, columns=cols)
        else:
            df=arr
        df.loc[:,'transcript_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('transcript_id',np.nan))
        df.loc[:,'gene_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('gene_id',np.nan))
        self.gtf_df = df
        self.gtf_summary()
        return self

    def write_file(self, out_file, comment_str):
        df = self.gtf_df
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        print(df.shape)
        with open(out_file, "w") as fout:
            fout.write(comment_str)
            for _, row in df.iterrows():
                fields = [
                    row["chrom"],
                    row["source"],
                    row["feature"],
                    str(int(row["start"])),
                    str(int(row["end"])),
                    row["score"],
                    row["strand"],
                    row["frame"],
                    row["attributes"]
                ]
                fout.write("\t".join(fields) + "\n")
        print(f" * Wrote corrected GTF to {out_file} ({len(df)} records)")

    def load_df(self):
        cols = self.cols
        df = pd.read_csv(self.gtf_file, sep='\t', comment='#', names=cols, dtype=str)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df.loc[:,'transcript_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('transcript_id',np.nan))
        df.loc[:,'gene_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('gene_id',np.nan))
        df.loc[:,'gene_name']=df['attributes'].map(lambda d: self.parse_attributes(d).get('gene_name',np.nan))
        self.gtf_df = df
        self.gtf_summary()
        return self

    def get_unique_gene_id(self, only_gene_future=False):
        anno = self.gtf_df
        known_genes = set()
        if only_gene_future:
            anno_gene_mask = anno['feature'] == 'gene'
        else:
            anno_gene_mask = anno.index
        for attr in anno.loc[anno_gene_mask, 'attributes']:
            d = self.parse_attributes(attr)
            if 'gene_id' in d:
                known_genes.add(d['gene_id'])
        return known_genes

    def get_transcript_gene_id_map(self,gene_item='gene_id'):
        anno = self.gtf_df
        transcript_gene_map = {}
        anno_gene_mask = anno['feature'] == 'transcript'
        for attr in anno.loc[anno_gene_mask, 'attributes']:
            d = self.parse_attributes(attr)
            if gene_item in d and 'transcript_id' in d:
                transcript_gene_map[d['transcript_id']] = d[gene_item]
        return transcript_gene_map

    def get_transcript_id_attr_map(self):
        anno = self.gtf_df
        transcript_gene_map = {}
        anno_gene_mask = anno['feature'] == 'transcript'
        for attr in anno.loc[anno_gene_mask, 'attributes']:
            d = self.parse_attributes(attr)
            transcript_gene_map[d['transcript_id']] = attr
        return transcript_gene_map

    def get_transcript_length_map(self):
        anno = self.gtf_df
        transcript_length_map = {}
        anno_gene_mask = anno['feature'] == 'transcript'
        for _,row in anno.loc[anno_gene_mask, :].iterrows():
            d = self.parse_attributes(row['attributes'])
            if 'transcript_id' in d:
                transcript_length_map[d['transcript_id']] = abs(row['start'] - row['end'])
        return transcript_length_map

    def get_attribute_map(self,feature='gene',id='gene_id'):
        anno = self.gtf_df
        attr_map = {}
        anno_gene_mask = anno['feature'] == feature
        for _,row in anno.loc[anno_gene_mask, :].iterrows():
            d = self.parse_attributes(row['attributes'])
            attr_map[d[id]] = d
        return attr_map

    def get_position_map(self,feature='gene',id='gene_id'):
        anno = self.gtf_df
        attr_map = {}
        anno_gene_mask = anno['feature'] == feature
        for _,row in anno.loc[anno_gene_mask, :].iterrows():
            d = self.parse_attributes(row['attributes'])
            attr_map[d[id]] = [row['start'], row['end']]
        return attr_map

    def gtf_summary(self):
        df  = self.gtf_df
        name = ''
        if self.gtf_file is not None:
            name = os.path.basename(self.gtf_file)
        genes = set()
        transcripts = set()
        exons = 0
        for attr, feat in zip(df['attributes'], df['feature']):
            d = self.parse_attributes(attr)
            if 'gene_id' in d:
                genes.add(d['gene_id'])
            if 'transcript_id' in d:
                transcripts.add(d['transcript_id'])
            if feat == 'exon':
                exons += 1
        print(f"\n[{name} GTF summary]")
        print(f"  Unique genes:       {len(genes):,}")
        print(f"  Unique transcripts: {len(transcripts):,}")
        print(f"  Exons:              {exons:,}")

    def correct_transcript_boundaries(self):
        df = self.gtf_df.copy()
        corrected_rows = []
        exons = df[df["feature"] == "exon"]
        for tid, sub in exons.groupby("transcript_id"):
            min_start = sub["start"].min()
            max_end = sub["end"].max()
            corrected_rows.append((tid, min_start, max_end))
        corrected_df = pd.DataFrame(corrected_rows, columns=["transcript_id", "new_start", "new_end"])
        df = df.merge(corrected_df, on="transcript_id", how="left")
        df.loc[df["feature"] == "transcript", "start"] = df.loc[df["feature"] == "transcript", "new_start"]
        df.loc[df["feature"] == "transcript", "end"] = df.loc[df["feature"] == "transcript", "new_end"]
        df.drop(columns=["new_start", "new_end"], inplace=True)
        ## correct start/end
        mask = df['start'] > df['end']
        if mask.any():
            logger.info(f"Warning: found {mask.sum()} lines with start > end; fixing...")
            df.loc[mask, ['start', 'end']] = df.loc[mask, ['end', 'start']].values
        logger.info(f" * Corrected transcript start/end using exon boundaries")
        self.gtf_df = df
        return

    def sort(self):
        '''
        sort gtf
        :return:
        '''
        def sort_chr_rule(chr:str):
            chr_map={'X':23,'Y':24,'M':25}
            if chr.startswith('chr'):
                chr_num=chr.split('chr')[-1]
                try:
                    chr_val=int(chr_num)
                except:
                    chr_val=chr_map[chr_num]
            else:
                chr_val=26+int(hashlib.md5(chr.encode('utf-8')).hexdigest()[:8], 16)
            return chr_val

        gene_pos=self.get_position_map('gene','gene_id')
        tran_pos=self.get_position_map('transcript','transcript_id')
        df = self.gtf_df.copy()
        df['chr_val']=df['chrom'].map(lambda x:sort_chr_rule(x))
        sort_data=[]
        for _,row in df.iterrows():
            attr=self.parse_attributes(row['attributes'])
            gene_id=attr['gene_id']
            gene_v=gene_pos[gene_id]
            if 'transcript_id' in attr:
                trans_v=tran_pos[attr['transcript_id']]
            else:
                trans_v=[-1,-1]
            if row['feature'] == 'exon':
                exons_v=[row['start'], row['end']]
            else:
                exons_v=[-1,-1]
            sort_row=gene_v+trans_v+exons_v
            sort_data.append(sort_row)
        sort_cols=["gene_s", "gene_e", "transcript_s","transcript_e","exon_s","exon_e"]
        sdf=pd.DataFrame(sort_data, columns=sort_cols)
        sdf.index=df.index
        mdf=pd.concat([df,sdf], axis=1)
        mdf.sort_values(by=['chr_val']+sort_cols, inplace=True)
        mdf.drop(['chr_val']+sort_cols, axis=1, inplace=True)
        self.gtf_df = mdf
        return self

    def update_gene_coor(self):
        df_updated = self.gtf_df.copy()
        df_updated['gene_id'] = df_updated['attributes'].apply(
            lambda x: self.parse_attributes(x).get('gene_id', None)
        )
        gene_rows = df_updated[df_updated['feature'] == 'gene'].copy()
        transcript_rows = df_updated[df_updated['feature'] == 'transcript'].copy()
        if not transcript_rows.empty:
            transcript_stats = transcript_rows.groupby('gene_id').agg(
                min_transcript_start=('start', 'min'),
                max_transcript_end=('end', 'max')
            ).reset_index()
            gene_rows_updated = pd.merge(
                gene_rows,
                transcript_stats,
                on='gene_id',
                how='left'
            )
            mask = gene_rows_updated['min_transcript_start'].notna()
            gene_rows_updated.loc[mask, 'start'] = gene_rows_updated.loc[mask, 'min_transcript_start']
            gene_rows_updated.loc[mask, 'end'] = gene_rows_updated.loc[mask, 'max_transcript_end']
            gene_rows_updated = gene_rows_updated.drop(['min_transcript_start', 'max_transcript_end'], axis=1)
        else:
            gene_rows_updated = gene_rows
        df_updated = df_updated[df_updated['feature'] != 'gene']
        df_updated = pd.concat([df_updated, gene_rows_updated], ignore_index=True)
        if 'gene_id' in df_updated.columns:
            df_updated = df_updated.drop('gene_id', axis=1)
        self.gtf_df = df_updated
        pass

    def filter_transcripts_with_genes(self, transcript_ids):
        gtf_df = self.gtf_df
        target_genes = gtf_df.loc[gtf_df['transcript_id'].isin(transcript_ids), 'gene_id'].unique()
        mask = gtf_df['transcript_id'].isin(transcript_ids) | (
                    (gtf_df['feature'] == 'gene') & gtf_df['gene_id'].isin(target_genes))
        return gtf_df.loc[mask, [x for x in gtf_df.columns if x not in ['transcript_id', 'gene_id']]]

    def filter_genes(self, gene_list):
        gtf_df = self.gtf_df
        mask = (gtf_df['feature'] == 'gene') & gtf_df['gene_id'].isin(gene_list)
        return gtf_df.loc[mask, [x for x in gtf_df.columns if x not in ['transcript_id', 'gene_id']]]



def __make_enhanced_gtf(class_file, gtf:GTFParser, ref_gtf:GTFParser, enhanced_gtf_path, enhanced_comment,
                      novel_gtf_path=None, novel_gtf_comment='', gtop_gtf_path=None, gtop_gtf_comment=''):
    annot_df = pd.read_csv(class_file, sep='\t')
    annot_df = annot_df.loc[annot_df["filter_result"] == "Isoform",:]
    ## correct associated_gene for known isoforms (FSM/ISM), because some isoform with _ join multi-genes.
    ref_isoform_gene_map=ref_gtf.get_transcript_gene_id_map()
    print(annot_df['associated_transcript'])
    annot_df['associated_gene']=annot_df.apply(lambda row:row['associated_gene'] if row['associated_transcript']=='novel'
    else ref_isoform_gene_map[row['associated_transcript']], axis=1)

    for col in ["isoform", "structural_category", "filter_result"]:
        if col not in annot_df.columns:
            raise ValueError(f"Missing required column '{col}' in class file")
    mnd_cata_full={np.nan:'noncoding_RNA',True:'nonsense_mediated_decay',False:'protein_coding'}

    annot_df['predicted_NMD.full']=annot_df['predicted_NMD'].map(lambda x:mnd_cata_full[x])
    # known = ["full-splice_match", "incomplete-splice_match"]
    annot_df_novel = annot_df.loc[annot_df['associated_transcript'] == 'novel', :]
    novel_isoform_ids = annot_df_novel['isoform'].unique()
    known_isoform_ids_in_all_isoforms=annot_df.loc[annot_df['associated_transcript'] != 'novel', 'associated_transcript'].unique()

    # get known gene id and other info, e.g., gene name, gene type.
    known_gene_attr=ref_gtf.get_attribute_map(feature='gene',id='gene_id')
    isoform_assoc_gene = dict(zip(annot_df_novel['isoform'], annot_df_novel['associated_gene']))
    # get transcript type
    novel_isoform_type = dict(zip(annot_df_novel["isoform"], annot_df_novel["predicted_NMD.full"]))
    isoform_length = gtf.get_transcript_length_map()
    isoform_gene = gtf.get_transcript_gene_id_map()
    logger.info(f'start make gene/transcript type')
    # formal transcript type
    formal_novel_isoform_type={}
    gene_type_map_arr={}
    for iso, gt in novel_isoform_type.items():
        if gt == 'noncoding_RNA':
            if isoform_length[iso] >= 200:
                gt='lncRNA'
            else:
                gt='other_ncRNA'
        formal_novel_isoform_type[iso]=gt
        gene=isoform_gene[iso]
        if gene not in gene_type_map_arr:
            gene_type_map_arr[gene]=[]
        gene_type_map_arr[gene].append(gt)
    # get gene type
    gene_type = {}
    for g,gts in gene_type_map_arr.items():
        gene_gt=None
        for gt in ['protein_coding','lncRNA','nonsense_mediated_decay','other_ncRNA']:
            if gts.count(gt) > 0:
                gene_gt=gt
                break
        gene_type[g]=gene_gt
    logger.info(f'start make novel isoform gtf')

    ldf=gtf.gtf_df
    sdf=ldf.loc[ldf['transcript_id'].isin(novel_isoform_ids), :].copy()
    parsed_attrs = sdf["attributes"].apply(GTFParser.parse_attributes)
    def change_id(attr:{}):
        tid=attr["transcript_id"]
        gid=attr["gene_id"]
        attr['transcript_type']=formal_novel_isoform_type[tid]
        annot_gid=isoform_assoc_gene[tid]
        if annot_gid in known_gene_attr.keys():
            for k,v in known_gene_attr[annot_gid].items():
                attr[k]=v
        # for novel gene
        else:
            attr['gene_id'] = gid
            attr['gene_name'] = gid
            attr['gene_type'] = gene_type[gid]
        return attr

    annot_geneid_in_novel_isoforms=set([isoform_assoc_gene[tid] for tid in sdf['transcript_id'] if isoform_assoc_gene[tid] in known_gene_attr])

    changed_attrs=parsed_attrs.apply(change_id)
    sdf['attributes']=changed_attrs.apply(GTFParser.build_attributes)

    novel_isoform_gtf=GTFParser()
    novel_isoform_gtf.create_from_arr(sdf)
    novel_isoform_gtf.gtf_summary()
    # add novel gene feature line
    ndf = novel_isoform_gtf.gtf_df.copy()
    novel_gene_rows=[]
    k=0
    logger.info(f'start make novel gene gtf record in novel isoform ')
    for gid, gdf in ndf.groupby(by='gene_id'):
        if gid not in known_gene_attr.keys():
            attr_map={'gene_id':gid,'gene_name':gid,'gene_type':gene_type[gid]}
            row = {
                'chrom': gdf.iloc[0]['chrom'],
                'source': 'GTOP',
                'feature': 'gene',
                'start': gdf['start'].min(),
                'end': gdf['end'].max(),
                'score': '.',
                'strand': gdf.iloc[0]['strand'],
                'frame': '.',
                'attributes': GTFParser.build_attributes(attr_map)
            }
            novel_gene_rows.append(row)
        k+=1
    logger.info(f'Novel isoform: all gene: {k}; novel gene: {len(novel_gene_rows)}')
    novel_genes_gtf_record_in_novel_isoform_df = pd.DataFrame(novel_gene_rows)
    # save novel gtf
    logger.info(f'Annotated genes: {len(annot_geneid_in_novel_isoforms)} in GTOP novel isoforms.')
    annot_genes_gtf_record_in_novel_isoform_df = ref_gtf.filter_genes(annot_geneid_in_novel_isoforms)
    logger.info(f'Get {annot_genes_gtf_record_in_novel_isoform_df.shape[0]} annotated genes record for ref gtf.')
    if novel_gtf_path is not None:
        merge_df=pd.concat([annot_genes_gtf_record_in_novel_isoform_df,novel_genes_gtf_record_in_novel_isoform_df,
                            novel_isoform_gtf.gtf_df],ignore_index=True) #ref_gtf.gtf_df,
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        # sort merged GTF
        logger.info(f'Novel: start sorting gtf')
        out_gtf.sort()
        logger.info(f'Novel: start writing gtf')
        out_gtf.write_file(novel_gtf_path,novel_gtf_comment)
        out_gtf.gtf_summary()
        logger.info(f"Novel: Writing GTF: {novel_gtf_path}")
    # save gtop gtf: novel + known (replace with GENCODE v47 transcript annotation).
    if gtop_gtf_path is not None:
        annot_isoform_df=ref_gtf.filter_transcripts_with_genes(known_isoform_ids_in_all_isoforms)
        logger.info(f'Known isoforms ID: {len(known_isoform_ids_in_all_isoforms)} in all GTOP isoforms.')
        merge_df=pd.concat([annot_isoform_df, annot_genes_gtf_record_in_novel_isoform_df,
                            novel_genes_gtf_record_in_novel_isoform_df,
                            novel_isoform_gtf.gtf_df],ignore_index=True) #ref_gtf.gtf_df,
        ## remove duplicate gene record
        mask = pd.Series(False, index=merge_df.index)
        gene_df = merge_df[merge_df['feature'] == 'gene']
        dup_gene = gene_df.duplicated(subset=['attributes'], keep='first')
        mask.loc[dup_gene.index] = dup_gene
        merge_df = merge_df[~((merge_df['feature'] == 'gene') & mask)]
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        out_gtf.gtf_summary()
        # sort merged GTF
        logger.info(f'GTOP: start sorting gtf')
        out_gtf.sort()
        logger.info(f'GTOP: start writing gtf')
        out_gtf.write_file(gtop_gtf_path,gtop_gtf_comment)
        logger.info(f"GTOP: Writing GTF: {gtop_gtf_path}")
    # save enhanced gtf: novel + GENCODE v47 transcript annotation.
    if enhanced_gtf_path is not None:
        merge_df=pd.concat([ref_gtf.gtf_df,novel_genes_gtf_record_in_novel_isoform_df,
                            novel_isoform_gtf.gtf_df],ignore_index=True) #ref_gtf.gtf_df,
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        # update gene coor with novel isoforms.
        logger.info(f'Enhanced: start updating gene coor')
        out_gtf.update_gene_coor()
        # sort merged GTF
        logger.info(f'Enhanced: start sorting gtf')
        out_gtf.sort()
        logger.info(f'Enhanced: start writing gtf')
        out_gtf.write_file(enhanced_gtf_path,enhanced_comment)
        out_gtf.gtf_summary()
        logger.info(f"Enhanced: Writing GTF: {enhanced_gtf_path}")

def __make_GTOP_gtf(class_file, gtf:GTFParser, ref_gtf:GTFParser, gtop_gtf_path=None, gtop_gtf_comment=''):
    annot_df = pd.read_csv(class_file, sep='\t')
    annot_df = annot_df.loc[annot_df["filter_result"] == "Isoform",:]
    ## correct associated_gene for known isoforms (FSM/ISM), because some isoform with _ join multi-genes.
    ref_isoform_gene_map=ref_gtf.get_transcript_gene_id_map()
    annot_df['associated_gene']=annot_df.apply(lambda row:row['associated_gene'] if row['associated_transcript']=='novel'
    else ref_isoform_gene_map[row['associated_transcript']], axis=1)

    for col in ["isoform", "structural_category", "filter_result"]:
        if col not in annot_df.columns:
            raise ValueError(f"Missing required column '{col}' in class file")
    mnd_cata_full={np.nan:'noncoding_RNA',True:'nonsense_mediated_decay',False:'protein_coding'}

    annot_df['predicted_NMD.full']=annot_df['predicted_NMD'].map(lambda x:mnd_cata_full[x])
    # known = ["full-splice_match", "incomplete-splice_match"]
    isoform_ids = annot_df['isoform'].unique()
    print(f'load {len(isoform_ids)} isoforms.')
    # get known gene id and other info, e.g., gene name, gene type.
    known_gene_attr=ref_gtf.get_attribute_map(feature='gene',id='gene_id')
    isoform_assoc_gene = dict(zip(annot_df['isoform'], annot_df['associated_gene']))
    # get transcript type
    isoform_type = dict(zip(annot_df["isoform"], annot_df["predicted_NMD.full"]))
    isoform_length = gtf.get_transcript_length_map()
    isoform_gene = gtf.get_transcript_gene_id_map()
    logger.info(f'start make gene/transcript type')
    # formal transcript type
    formal_isoform_type={}
    gene_type_map_arr={}
    for iso, gt in isoform_type.items():
        if gt == 'noncoding_RNA':
            if isoform_length[iso] >= 200:
                gt='lncRNA'
            else:
                gt='other_ncRNA'
        formal_isoform_type[iso]=gt
        gene=isoform_gene[iso]
        if gene not in gene_type_map_arr:
            gene_type_map_arr[gene]=[]
        gene_type_map_arr[gene].append(gt)
    # get gene type
    gene_type = {}
    for g,gts in gene_type_map_arr.items():
        gene_gt=None
        for gt in ['protein_coding','lncRNA','nonsense_mediated_decay','other_ncRNA']:
            if gts.count(gt) > 0:
                gene_gt=gt
                break
        gene_type[g]=gene_gt
    logger.info(f'start make novel isoform gtf')
    data=[]
    ldf=gtf.gtf_df
    sdf=ldf.loc[ldf['transcript_id'].isin(isoform_ids), :].copy()
    parsed_attrs = sdf["attributes"].apply(GTFParser.parse_attributes)
    def change_id(attr:{}):
        tid=attr["transcript_id"]
        gid=attr["gene_id"]
        annot_gid=isoform_assoc_gene[tid]
        if annot_gid in known_gene_attr.keys():
            for k,v in known_gene_attr[annot_gid].items():
                attr[k]=v
        # for novel gene
        else:
            attr['gene_id'] = gid
            attr['gene_name'] = gid
            attr['gene_type'] = gene_type[gid]
        return attr

    annot_geneid_in_isoforms=set([isoform_assoc_gene[tid] for tid in sdf['transcript_id'] if isoform_assoc_gene[tid] in known_gene_attr])

    changed_attrs=parsed_attrs.apply(change_id)
    sdf['attributes']=changed_attrs.apply(GTFParser.build_attributes)

    isoform_gtf=GTFParser()
    isoform_gtf.create_from_arr(sdf)
    isoform_gtf.gtf_summary()
    # add novel gene feature line
    ndf = isoform_gtf.gtf_df.copy()
    gene_rows=[]
    k=0
    logger.info(f'start make novel gene gtf record in novel isoform ')
    for gid, gdf in ndf.groupby(by='gene_id'):
        if gid not in known_gene_attr.keys():
            attr_map={'gene_id':gid,'gene_name':gid,'gene_type':gene_type[gid]}
            row = {
                'chrom': gdf.iloc[0]['chrom'],
                'source': 'GTOP',
                'feature': 'gene',
                'start': gdf['start'].min(),
                'end': gdf['end'].max(),
                'score': '.',
                'strand': gdf.iloc[0]['strand'],
                'frame': '.',
                'attributes': GTFParser.build_attributes(attr_map)
            }
            gene_rows.append(row)
        k+=1
    logger.info(f'GTOP isoform: all gene: {k}; novel gene: {len(gene_rows)}')
    novel_genes_gtf_record_in_novel_isoform_df = pd.DataFrame(gene_rows)
    # save novel gtf
    logger.info(f'Annotated genes: {len(annot_geneid_in_isoforms)} in GTOP isoforms.')
    annot_genes_gtf_record_in_isoform_df = ref_gtf.filter_genes(annot_geneid_in_isoforms)
    logger.info(f'Get {annot_genes_gtf_record_in_isoform_df.shape[0]} annotated genes record for ref gtf.')
    if gtop_gtf_path is not None:
        merge_df=pd.concat([annot_genes_gtf_record_in_isoform_df,novel_genes_gtf_record_in_novel_isoform_df,
                            isoform_gtf.gtf_df],ignore_index=True) #ref_gtf.gtf_df,
        out_gtf=GTFParser()
        out_gtf.create_from_arr(merge_df)
        # sort merged GTF
        logger.info(f'GTOP: start sorting gtf')
        out_gtf.sort()
        logger.info(f'GTOP: start writing gtf')
        out_gtf.write_file(gtop_gtf_path,gtop_gtf_comment)
        out_gtf.gtf_summary()
        logger.info(f"GTOP: Writing GTF: {gtop_gtf_path}")

def __merge_cds_to_enhanced_gtf(cds_path,enhanced_gtf_path,out_path,enhanced_with_CDS_comment):
    cds=GTFParser(cds_path).load_df()
    gtf=GTFParser(enhanced_gtf_path).load_df()
    novel_iso=[g for g in gtf.get_transcript_gene_id_map().keys() if not g.startswith('ENS')]
    logger.info(f'load {len(novel_iso)} novel isoforms')
    cds_df=cds.gtf_df
    logger.info(f'start merge')
    cds_exon_df=cds_df.loc[(cds_df['feature']=='CDS') & (cds_df['transcript_id'].isin(novel_iso)), :]
    trans_attr=gtf.get_transcript_id_attr_map()
    cds_exon_df['attributes']=cds_exon_df['transcript_id'].apply(lambda x: trans_attr[x])
    df=pd.concat([gtf.gtf_df,cds_exon_df],ignore_index=True)
    final_gtf=GTFParser().create_from_arr(df)
    logger.info(f'start write')
    final_gtf.write_file(out_path,f'{enhanced_with_CDS_comment}')


class GTFAnnotation:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file

    def conduct_gene_length(self):
        gtf_file = self.gtf_file
        types = {'gene': '-l', 'transcript': '-r'}
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        for tk, tv in types.items():
            cmd = f'python {CURRENT_DIR}/gtftools.py {tv} {wdir}/{f[:-4]}.{tk}_length {wdir}/{f}'
            print(cmd)
            os.system(cmd)

    def make_gene_transcript_tab(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_file = f'{wdir}/{f[:-4]}.gene_transcript_map.txt'
        seen = set()
        with open(gtf_file) as gtf, open(out_file, "w") as out:
            out.write("transcript_id\tchr\tstart\tend\tgene_id\tgene_name\ttranscript_type\tgene_type\n")
            for line in gtf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                attr_field = fields[8]
                attr_dict = {}
                for attr in attr_field.strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        continue
                tid = attr_dict.get("transcript_id")
                gid = attr_dict.get("gene_id")
                gname = attr_dict.get("gene_name")
                transcript_type = attr_dict.get("transcript_type")
                gene_type = attr_dict.get("gene_type")
                chr = fields[0]
                start = fields[3]
                end = fields[4]
                if tid and gid and tid not in seen:
                    out.write(f"{tid}\t{chr}\t{start}\t{end}\t{gid}\t{gname}\t{transcript_type}\t{gene_type}\n")
                    seen.add(tid)

    def make_gene_annot_tab(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_file = f'{wdir}/{f[:-4]}.gene_annot.txt'
        seen = set()
        with open(gtf_file) as gtf, open(out_file, "w") as out:
            for line in gtf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                if fields[2].strip() != 'gene':
                    continue
                attr_field = fields[8]
                attr_dict = {}
                for attr in attr_field.strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        continue
                gid = attr_dict.get("gene_id")
                gname = attr_dict.get("gene_name")
                gtype = attr_dict.get("gene_type")
                chr = fields[0]
                st = fields[3]
                ed = fields[4]
                chain = fields[6]
                out.write(f"{gname}\t{gtype}\t{gid}\t{chr}\t{st}\t{ed}\t{chain}\n")
        logger.info(f'complete {gtf_file}')

    def make_autosome_transcript_tab(self):
        keep_autosomes_only = True
        protein_coding_only = False
        # wdir = f'/media/dubai/home/xuechao/project/GMTiP-RNA/20251031/long_read/HPC/output/enhanced_gtf'
        # files=['GTOP_novel-GENCODE_v47.gtf','GTOP_all.gtf']
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_tab = f'{wdir}/{f[:-4]}.autosome_isoform.csv'
        gtf = pd.read_csv(
            gtf_file,
            sep="\t",
            comment="#",
            header=None,
            names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        )
        transcripts = gtf[gtf["feature"].isin(["transcript", "mRNA"])]  #
        if protein_coding_only:
            transcripts = transcripts[transcripts["attribute"].str.contains('gene_biotype "protein_coding"')]
        transcripts["gene_id"] = transcripts["attribute"].str.extract('gene_id "([^"]+)"')
        transcripts["transcript_id"] = transcripts["attribute"].str.extract('transcript_id "([^"]+)"')

        def normalize_chr(c):
            if pd.isna(c):
                return c
            c = str(c).strip()
            if c.startswith("chr"):
                return c
            if re.match(r"^\d+$", c):
                return "chr" + c
            if c.lower() in ["m", "mt", "chrM", "chrMT", "mitochondria", "mitochondrial"]:
                return "chrM"
            # sex chromosomes
            if c.lower() in ["x"]:
                return "chrX"
            if c.lower() in ["y"]:
                return "chrY"
            return c

        transcripts["chromosome"] = transcripts["chr"].apply(normalize_chr)
        if keep_autosomes_only:
            autosomes = [f"chr{i}" for i in range(1, 23)]
            transcripts = transcripts[transcripts["chromosome"].isin(autosomes)]
        output_df = transcripts[["transcript_id", "gene_id", "chromosome"]].drop_duplicates()
        output_df.to_csv(out_tab, index=False)
        logger.info(f'save to {out_tab}')

    def index_gtf(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        gtf_kw = f.split('.')[0]
        sorted_gtf = f'{wdir}/{gtf_kw}.sorted.gtf'
        gz_gtf = f'{wdir}/{gtf_kw}.sorted.gtf.gz'
        logger.info(f"Sorting {gtf_file}...")
        subprocess.run(f"sort -k1,1 -k4,4n {gtf_file} > {sorted_gtf}", shell=True, check=True)
        logger.info(f"Compressing {sorted_gtf}...")
        subprocess.run(f"bgzip -c {sorted_gtf} > {gz_gtf}", shell=True, check=True)
        logger.info(f"Indexing {gz_gtf}...")
        subprocess.run(f"tabix -p gff {gz_gtf}", shell=True, check=True)
        logger.info("Done!")
        os.remove(sorted_gtf)
        logger.info(f"Removed {sorted_gtf}")
        pass

    def make_gtf_fa(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        f_prefix = f[:-4]
        cmd = f'{GFFREAD_BIN} {gtf_file} -g {REF_GENOME_FASTA} -w {wdir}/{f_prefix}.fa'
        logger.info(f"cmd: {cmd}")
        subprocess.run(cmd, shell=True)

    def annotate_gtf_task(self):
        # self.make_gtf_fa()
        # self.make_gene_annot_tab()
        # self.conduct_gene_length()
        # self.make_gene_transcript_tab()
        # self.make_autosome_transcript_tab()
        self.index_gtf()


def make_enhanced_gtf_and_annotation_task(LRS_out_dir):
    enhanced_gtf_dir=f'{LRS_out_dir}/enhanced_gtf'
    sqanti_dir=f'{LRS_out_dir}/sqanti3_clean'
    sqanti3_gtf=f'{sqanti_dir}/filter.clean.gtf'
    class_file=f'{sqanti_dir}/filter.clean.RulesFilter_result_classification.txt'
    cds_gtf=f'{sqanti_dir}/filter.clean.cds.gff3'
    ref_gtf=REF_GTF
    novel_gtf=f'{enhanced_gtf_dir}/GTOP_novel.gtf'
    novel_gtf_comment = f'# GTOP novel isoform annotation \n'
    gtop_gtf=f'{enhanced_gtf_dir}/GTOP.gtf'
    gtop_comment = f'# GTOP isoform annotation: GTOP novel + annotated \n'
    gtop_all_gtf=f'{enhanced_gtf_dir}/GTOP_all.gtf'
    gtop_all_comment = f'# GTOP isoform annotation: GTOP all \n'
    enhanced_gtf = f'{enhanced_gtf_dir}/GTOP_novel-GENCODE_v47.gtf'
    enhanced_comment = f'# Enhanced isoform annotation: GTOP novel + GENCODE v47 \n'
    enhanced_with_CDS_gtf = f'{enhanced_gtf_dir}/GTOP_novel-GENCODE_v47.with_CDS.gtf'
    enhanced_with_CDS_comment = f'# Enhanced isoform annotation: GTOP novel + GENCODE v47 with CDS\n'
    os.makedirs(os.path.dirname(enhanced_gtf), exist_ok=True)
    #
    # logger.info(f"Loading annotated GTF: {sqanti3_gtf}")
    # gtf_par=GTFParser(sqanti3_gtf).load_df()
    # # correct gtf generated by SQANTI3, because the start/end position of transcript in "-" chain is error.
    # gtf_par.correct_transcript_boundaries()
    # gtf_par.gtf_summary()
    #
    # logger.info(f"Loading reference GTF: {ref_gtf}")
    # ref_gtf_par=GTFParser(ref_gtf).load_df()
    # ref_gtf_par.gtf_summary()
    # logger.info(f"Merging GTF")
    # __make_enhanced_gtf(class_file,gtf_par,ref_gtf_par,enhanced_gtf,enhanced_comment,
    #                   novel_gtf ,novel_gtf_comment,
    #                   gtop_gtf,gtop_comment)
    # # generate GTOP raw isoform with known genes.
    # __make_GTOP_gtf(class_file,gtf_par,ref_gtf_par, gtop_all_gtf,gtop_all_comment)
    # logger.info("✅ Done. Novel isoforms GTF written to:", gtop_all_gtf)
    # # merge CDS to enhanced gtf.
    # __merge_cds_to_enhanced_gtf(cds_gtf,enhanced_gtf,enhanced_with_CDS_gtf,enhanced_with_CDS_comment)
    # logger.info("✅ Done. CDS written to:", enhanced_with_CDS_gtf)
    # annotate GTF
    cus_gtf=['GTOP_novel-GENCODE_v47.gtf','GTOP_all.gtf','GTOP.gtf']
    gtf_files=[]
    for g in cus_gtf:
        gtf_files.append(f'{enhanced_gtf_dir}/{g}')
    gtf_files.append(REF_GTF)
    for gtf in gtf_files:
        logger.info(f"start to annotate {gtf}")
        GTFAnnotation(gtf).annotate_gtf_task()


import logging

def extract_transcripts_gtf(gtf_path, transcript_ids, out_path):
    ids = set(transcript_ids)
    all_tids = set()
    kept_tids = set()
    with open(gtf_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            tid = None
            for item in parts[8].split(";"):
                item = item.strip()
                if item.startswith("transcript_id"):
                    tid = item.split(" ", 1)[1].strip('"')
                    break
            if tid is None:
                continue
            all_tids.add(tid)
            if tid in ids:
                kept_tids.add(tid)
                fout.write(line)
    logger.info(f"Total transcripts in GTF: {len(all_tids)}")
    logger.info(f"Requested transcripts: {len(ids)}")
    logger.info(f"Successfully extracted transcripts: {len(kept_tids)}")


def fix_transcript_coords(gtf_in, gtf_out):
    coords = {}
    records = []

    with open(gtf_in) as f:
        for line in f:
            if line.startswith("#"):
                records.append((None, line))
                continue
            p = line.rstrip("\n").split("\t")
            tid = None
            for x in p[8].split(";"):
                x = x.strip()
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                    break
            records.append((tid, p))
            if p[2] == "exon" and tid is not None:
                coords.setdefault(tid, []).extend((int(p[3]), int(p[4])))

    with open(gtf_out, "w") as f:
        for tid, rec in records:
            if tid is None:
                f.write(rec)
            else:
                if rec[2] == "transcript" and tid in coords:
                    rec[3], rec[4] = map(str, (min(coords[tid]), max(coords[tid])))
                f.write("\t".join(rec) + "\n")


def transfer_gene_id(ref_gtf, in_gtf, out_gtf):
    logger = logging.getLogger(__name__)
    tx2gene = {}
    genes = set()

    with open(ref_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            tid = gid = None
            for x in p[8].split(";"):
                x = x.strip()
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                elif x.startswith("gene_id"):
                    gid = x.split(" ", 1)[1].strip('"')
            if tid and gid and tid not in tx2gene:
                tx2gene[tid] = gid
                genes.add(gid)

    logger.info(f"Transcripts found in reference GTF: {len(tx2gene)}")
    logger.info(f"Genes found in reference GTF: {len(genes)}")

    updated = set()
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                fout.write(line)
                continue

            tid = None
            attrs = []
            old_gene = None

            for x in p[8].split(";"):
                x = x.strip()
                if not x:
                    continue
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                    attrs.append(x)
                elif x.startswith("gene_id"):
                    old_gene = x
                else:
                    attrs.append(x)

            if tid in tx2gene:
                attrs.insert(0, f'gene_id "{tx2gene[tid]}"')
                updated.add(tid)
            elif old_gene is not None:
                attrs.insert(0, old_gene)

            p[8] = "; ".join(attrs) + ";"
            p[1] = 'GTOP'
            fout.write("\t".join(p) + "\n")

    logger.info(f"Transcripts updated or filled: {len(updated)}")

import re

def replace_gtf_ids(in_gtf, out_gtf):
    """
    Replace transcript_id and gene_id in a GTF file.

    Rules:
    - transcript_id: transcript_number -> GTOPTXXXXXXXX
    - gene_id:
        1) transcript_number -> same as transcript_id
        2) LOC_number -> GTOPGXXXXXXXX
    """

    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"

    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue

            attr = fields[8]

            # ---- transcript_id ----
            m_tx = re.search(r'transcript_id "([^"]+)"', attr)
            tx_new = None
            if m_tx:
                tx_id = m_tx.group(1)
                if tx_id.startswith("transcript_"):
                    num = tx_id.split("_", 1)[1]
                    tx_new = format_id("GTOPT", num)
                    attr = re.sub(
                        r'transcript_id "[^"]+"',
                        f'transcript_id "{tx_new}"',
                        attr
                    )

            # ---- gene_id ----
            m_gene = re.search(r'gene_id "([^"]+)"', attr)
            if m_gene:
                gene_id = m_gene.group(1)

                if gene_id.startswith("transcript_") and tx_new is not None:
                    # rule 1: same as transcript_id
                    attr = re.sub(
                        r'gene_id "[^"]+"',
                        f'gene_id "{tx_new}"',
                        attr
                    )

                elif gene_id.startswith("LOC_"):
                    # rule 2: LOC_number -> GTOPGXXXXXXXX
                    num = gene_id.split("_", 1)[1]
                    gene_new = format_id("GTOPG", num)
                    attr = re.sub(
                        r'gene_id "[^"]+"',
                        f'gene_id "{gene_new}"',
                        attr
                    )

            fields[8] = attr
            fout.write("\t".join(fields) + "\n")

def replace_sqanti3_annot_txid(in_annot, out_annot, tool_support_path):
    os.makedirs(os.path.dirname(out_annot),exist_ok=True)
    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"
    df=pd.read_csv(in_annot, sep='\t')
    # add tool support info
    tool_df=pd.read_csv(tool_support_path, sep='\t', index_col=0)
    tool_df["support_tools"] = tool_df.apply(lambda row: ",".join([tool for tool in tool_df.columns if row[tool] == 1]), axis=1)
    df['support_tools']=df['isoform'].apply(lambda x:tool_df.loc[x,'support_tools'])
    df['isoform']=df['isoform'].apply(lambda x:format_id('GTOPT',x.split('_')[1]))
    df.to_csv(out_annot, sep='\t', index=False)

def replace_faa_ids(in_faa, out_faa):
    os.makedirs(os.path.dirname(out_faa),exist_ok=True)
    def replace_tid(line):
        return re.sub(
            r'transcript_(\d+)',
            lambda m: f'GTOPT{int(m.group(1)):09d}',
            line
        )
    with open(in_faa, 'r') as fin, open(out_faa, "w") as fout:
        lines = []
        for line in fin:
            val=line.rstrip()
            if val.startswith(">"):
                val=replace_tid(val)
            lines.append(val+'\n')
        fout.writelines(lines)
    logger.info(f'save faa to {out_faa}')


def replace_flair_quant_ids(in_quant, clean_quant):
    df=pd.read_csv(in_quant, sep='\t',index_col=0)
    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"
    df.index=df.index.map(lambda x: format_id('GTOPT',x.split('_')[1]))
    os.makedirs(os.path.dirname(clean_quant),exist_ok=True)
    df.to_csv(clean_quant, sep='\t')
    pass


def assign_new_gene_loci(LRS_out_dir):
    tag='.custom'
    sqanti3_annot_path=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.RulesFilter_result_classification.txt'
    in_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.gtf'
    in_faa=f'{LRS_out_dir}/sqanti3_merged/filtered.faa'
    in_cds_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.cds.gff3'
    fixed_in_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.fixed.gtf'
    isoform_in_novel_gene_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.fixed.isoform_in_novel_gene.gtf'
    novel_gene_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.fixed.buildLoci_novel_gene.gtf'
    out_gtf=f'{LRS_out_dir}/sqanti3_merged/filtered{tag}.final.gtf'
    in_quant=f'{LRS_out_dir}/flair_quant/combined/raw_gtf.transcript.count.flair.tsv'
    #
    # # fix gtf, because transcript coor from SQANTI3 is wrong.
    # logger.info(f'start fix gtf')
    # fix_transcript_coords(in_gtf, fixed_in_gtf)
    # # get known gene id
    # ga=GeneAnnotation().load_df(REF_GTF[:-4],remove_ver=False)
    # gids=ga.get_gene_id_name_map().keys()
    # txids=ga.get_isoform_gene_map().keys()
    # logger.info(f'load {len(gids)} genes from GENCODE')
    # # extract isoforms in novel gene
    # sq_df = pd.read_csv(sqanti3_annot_path,sep='\t',index_col=0)
    # sq_df = sq_df.loc[sq_df['filter_result']=='Isoform',:]
    # print(sq_df.shape)
    # sq_df=sq_df.loc[sq_df['associated_transcript']=='novel', :]
    # isoform_ING_ids=sq_df.loc[~sq_df['associated_gene'].isin(gids), :].index.tolist()
    # logger.info(f'load {len(isoform_ING_ids)} isoforms in novel gene')
    # # save to tmp gtf
    # extract_transcripts_gtf(fixed_in_gtf,isoform_ING_ids,isoform_in_novel_gene_gtf)
    # # build gene loci
    # buildLoci_cmd='perl /lustre/home/cxue/software/buildLoci/buildLoci.pl'
    # bedtools_cmd='/lustre/home/cxue/software/bedtools2/bedtools'
    # cmd=f'{bedtools_cmd} intersect -s -wao -a {isoform_in_novel_gene_gtf} -b {isoform_in_novel_gene_gtf} | {buildLoci_cmd} - > {novel_gene_gtf}'
    # print(cmd)
    # os.system(cmd)
    # # assign novel gene id to gtf
    # transfer_gene_id(novel_gene_gtf,fixed_in_gtf,out_gtf)
    #
    # # rename gene id/transcript id
    tool_support_path=f'{LRS_out_dir}/raw_isoform/raw_merged.meta.tsv'
    clean_sqanti3_dir=f'{LRS_out_dir}/sqanti3_clean'
    os.makedirs(clean_sqanti3_dir,exist_ok=True)
    clean_gtf=f'{clean_sqanti3_dir}/filter.clean.gtf'
    clean_cds_gtf=f'{clean_sqanti3_dir}/filter.clean.cds.gff3'
    clean_annot=f'{clean_sqanti3_dir}/filter.clean.RulesFilter_result_classification.txt'
    clean_faa=f'{clean_sqanti3_dir}/filter.clean.faa'
    clean_quant=f'{clean_sqanti3_dir}/filter.clean.flair_quant.tsv'
    # replace_gtf_ids(out_gtf,clean_gtf)
    # replace_gtf_ids(in_cds_gtf,clean_cds_gtf)
    # replace_sqanti3_annot_txid(sqanti3_annot_path,clean_annot,tool_support_path)
    # replace_faa_ids(in_faa,clean_faa)
    replace_flair_quant_ids(in_quant, clean_quant)
    pass


if __name__ == '__main__':
    # ARGS = sys.argv[1:]
    # STEP = ARGS[0]
    # CONF_CSV = ARGS[1]
    # ROOT_OUTPUT_DIR = ARGS[2]
    # LOG_NAME = ARGS[3]
    # N_TASK = int(ARGS[4])
    # NT_PER_TASK = int(ARGS[5])
    # # single task
    # if STEP == 'enhanced_gtf':
    #     make_enhanced_gtf_and_annotation_task(ROOT_OUTPUT_DIR)

    LRS_out_dir=f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged'
    assign_new_gene_loci(LRS_out_dir)
    # make_enhanced_gtf_and_annotation_task(LRS_out_dir)