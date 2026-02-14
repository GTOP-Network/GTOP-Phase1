# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/7 17:58
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import hashlib
import logging
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from joblib import Parallel, delayed

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

PROJ_DIR = os.environ.get("PROJECT_ROOT")
LR_meta_path=f'{PROJ_DIR}/raw_data/GMTiP/meta/RNA/LRS-RNA-177-sample-annot.xlsx'
Tissue_meta_path=f'{PROJ_DIR}/raw_data/GMTiP/meta/Clinical/tissue_code.csv'


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
        df.loc[:,'transcript_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('transcript_id',pd.NA))
        df.loc[:,'gene_id']=df['attributes'].map(lambda d: self.parse_attributes(d).get('gene_id',pd.NA))
        self.gtf_df = df
        self.gtf_summary()
        return self

    def write_file(self, out_file, comment_str):
        df = self.gtf_df
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
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
        for t in ['transcript_id', 'gene_id', 'transcript_name','gene_type','transcript_type']:
            df.loc[:,t]=df['attributes'].map(lambda d: self.parse_attributes(d).get(t,np.nan))
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

    def get_gene_info_map(self,gene_item='gene_id'):
        anno = self.gtf_df
        gene_map = {}
        anno_gene_mask = anno['feature'] == 'gene'
        for attr in anno.loc[anno_gene_mask, 'attributes']:
            d = self.parse_attributes(attr)
            if gene_item in d:
                gene_map[d[gene_item]] = d
        return gene_map

    def get_transcript_gene_id_map(self,gene_item='gene_id'):
        anno = self.gtf_df
        transcript_gene_map = {}
        anno_gene_mask = anno['feature'] == 'transcript'
        for attr in anno.loc[anno_gene_mask, 'attributes']:
            d = self.parse_attributes(attr)
            if gene_item in d and 'transcript_id' in d:
                transcript_gene_map[d['transcript_id']] = d[gene_item]
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


class GeneAnnotation:
    def __init__(self):
        self.df=None
        pass

    def load_df(self, gene_annot_prefix, remove_ver=True):
        path=f'{gene_annot_prefix}.gene_transcript_map.txt'
        df=pd.read_csv(path,sep='\t',dtype=str)
        if remove_ver:
            for gid in ['gene_id','transcript_id']:
                df[gid]=df[gid].apply(lambda x:x.split('.')[0] if x.startswith('ENS') else x)
        self.df=df
        return self
        pass
    def get_gene_id_name_map(self):
        return dict(zip(self.df['gene_id'],self.df['gene_name']))

    def get_gene_isoforms_map(self,gene_key='gene_id'):
        gene_isoform_map = self.df.groupby(gene_key)['transcript_id'].apply(list).to_dict()
        return gene_isoform_map

    def get_isoform_gene_map(self,gene_key='gene_id',no_version=False):
        self.df['transcript_id_no_version']=self.df['transcript_id'].map(lambda x:x.split('.')[0] if x.startswith('ENS') else x)
        if no_version:
            gene_isoform_map = dict(zip(self.df['transcript_id_no_version'], self.df[gene_key]))
        else:
            gene_isoform_map = dict(zip(self.df['transcript_id'], self.df[gene_key]))
        return gene_isoform_map

    def get_isoform_pos(self):
        self.df['pos_str']=self.df['chr']+':'+self.df['start']+'-'+self.df['end']
        return dict(zip(self.df['transcript_id'],self.df['pos_str']))

    def get_isoform_type_map(self):
        return dict(zip(self.df['transcript_id'],self.df['transcript_type']))

    def get_gene_type_map(self):
        return dict(zip(self.df['gene_id'],self.df['gene_type']))


class ExprCorr:
    def __init__(self):
        pass

    def align_expression_matrices(
            self,
            df1: pd.DataFrame,
            df2: pd.DataFrame
    ):
        common_genes = df1.index.intersection(df2.index)
        common_samples = df1.columns.intersection(df2.columns)
        if len(common_genes) == 0:
            raise ValueError("No overlapping genes between df1 and df2")

        df1 = df1.loc[common_genes, common_samples]
        df2 = df2.loc[common_genes, common_samples]
        print(df1.shape)
        return df1, df2

    def spearman_1d(self, x: np.ndarray, y: np.ndarray, min_value=0.1) -> float:
        mask = (x > min_value) & (y > min_value)
        x_filtered = x[mask]
        y_filtered = y[mask]
        x_transformed = np.log2(x_filtered + 1)
        y_transformed = np.log2(y_filtered + 1)
        rx = rankdata(x_transformed, method="average")
        ry = rankdata(y_transformed, method="average")
        rx -= rx.mean()
        ry -= ry.mean()
        numerator = rx @ ry
        denominator = np.sqrt((rx @ rx) * (ry @ ry))
        return numerator / denominator if denominator != 0 else 0.0


    def spearman_by_sample_fast(
            self,
            df1: pd.DataFrame,
            df2: pd.DataFrame,
            n_jobs: int = 32
    ) -> pd.DataFrame:
        df1, df2 = self.align_expression_matrices(df1, df2)
        X = df1.values
        Y = df2.values
        sample_ids = df1.columns.to_numpy()
        spearman_vals = Parallel(
            n_jobs=n_jobs,
            backend="loky",
            batch_size=1
        )(
            delayed(self.spearman_1d)(X[:, i], Y[:, i])
            for i in range(X.shape[1])
        )
        return pd.DataFrame({
            "sample_id": sample_ids,
            "spearman_r": spearman_vals
        })

class GTOP_LRS_Meta:
    def __init__(self):
        self.LR_meta_path=LR_meta_path
        self.meta_map = None
        self.avail_keys = None
        self.df = None
        self.load_meta()

    def load_meta(self):
        rna_meta_df = pd.read_excel(self.LR_meta_path, dtype=str)
        self.df=rna_meta_df
        keys=['RIN','Age','Sex','Tissue','Batch']
        self.avail_keys=keys
        meta_map={}
        for k in keys:
            meta_map[k]=dict(zip(rna_meta_df['Sample_ID'],rna_meta_df[k]))
        self.meta_map=meta_map

    def get_meta_by_sample_id(self,sample_id,key='RIN'):
        arr=sample_id.split('-')
        sid=sample_id
        if len(arr)>2:
            sid='-'.join(arr[1:3])
        return self.meta_map[key][sid]

class GTOP_Tissue_Meta:
    def __init__(self):
        self.df = None
        self.tissue_name_meta_map = None
        self.tissue_code_meta_map = None
        self.tissue_meta_path=Tissue_meta_path
        self.load_meta()

    def load_meta(self):
        meta_df = pd.read_csv(self.tissue_meta_path, dtype=str)
        self.df=meta_df
        tissue_code_meta_map={}
        for k in meta_df.columns:
            tissue_code_meta_map[k]=dict(zip(meta_df['Tissue_Code'],meta_df[k]))
        tissue_name_meta_map={}
        for k in meta_df.columns:
            tissue_name_meta_map[k]=dict(zip(meta_df['Tissue'],meta_df[k]))
        self.tissue_code_meta_map=tissue_code_meta_map
        self.tissue_name_meta_map=tissue_name_meta_map

    def get_meta_by_tissue_code(self,tissue_code,key='Tissue'):
        return self.tissue_code_meta_map[key][tissue_code]

    def get_meta_by_tissue_name(self,tissue_name,key='Tissue_Color_Code'):
        return self.tissue_name_meta_map[key][tissue_name]


if __name__ == '__main__':
    logger.info('info')
    logger.warning('warning')
    logger.error('error')
    logger.critical('critical')
