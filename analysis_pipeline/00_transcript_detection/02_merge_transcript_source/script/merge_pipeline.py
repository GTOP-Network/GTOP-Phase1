# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/29 21:19
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import logging
import os

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


REF_DIR='/lustre/home/cxue/raw_data/GMTiP/ref/LRS'
REF_GENOME_FASTA = f"{REF_DIR}/genome.fa"

## common para
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_cerberus_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate cerberus'
gffread_PATH='/lustre/home/cxue/software/gffread-0.12.7.Linux_x86_64/gffread'

PROJ_DIR = os.environ.get("PROJECT_ROOT")
main_dir = f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery'
merge_main_dir=f'{main_dir}/merged'
gtf_dir=f'{merge_main_dir}/raw_isoform'
gtf_conf_path = f'{gtf_dir}/gtfs_path.csv'

raw_merged_gtf = f'{gtf_dir}/raw_merged.gtf'
raw_merged_meta_tsv=f'{gtf_dir}/raw_merged.meta.tsv'
tool2_filter_gtf = f'{gtf_dir}/tool2_filter.gtf'
tool2_filter_gtf_chr_dir = f'{merge_main_dir}/chr_based'


def make_gtf_conf():
    # tools_dir=f'{PROJ_DIR}/project/GMTiP-RNA/20251031/long_read/HPC/output_lrw/isoform_discovery'
    gtfs = {
        'bambu': f'{main_dir}/bambu/raw_isoform/bambu.gtf',
        'isoseq': f'{main_dir}/isoseq/raw_isoform/isoseq.gtf',
        'flair': f'{main_dir}/flair/raw_isoform/flair.gtf',
    }
    data=[]
    for k,v in gtfs.items():
        data.append([k, v])
    df=pd.DataFrame(data,columns=['tool_name','gtf_path'])
    os.makedirs(os.path.dirname(gtf_conf_path), exist_ok=True)
    df.to_csv(gtf_conf_path, index=False)
    pass

def merge_gtf_by_intron_chain():
    cmd=(f'{LOAD_cerberus_ENVS_CMD} '
         f'&& python {CURRENT_DIR}/merge_gtf_by_ic.py {gtf_conf_path} {raw_merged_gtf}')
    print(cmd)
    os.system(cmd)
    pass


def stat_isoform_tools_support():
    def parse_attributes(attr_str: str) -> dict:
        """
        Parse GFF3 or GTF attributes column into dict
        """
        attrs = {}
        for item in attr_str.strip().strip(";").split(";"):
            if not item:
                continue
            item = item.strip()
            if "=" in item:  # GFF3
                k, v = item.split("=", 1)
            else:  # GTF
                k, v = item.split(" ", 1)
                v = v.strip('"')
            attrs[k.strip()] = v.strip()
        return attrs

    import polars as pl

    def gtf_tools_to_polars_table(
            gff_file: str,
            feature_type: str = "transcript",
            transcript_key: str = "transcript_id",  # GTF: transcript_id | GFF3: ID
            tools_key: str = "tools",
    ):
        """
        Convert GFF/GTF with tools field to transcript × tool one-hot table
        """

        rows = []

        with open(gff_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) != 9:
                    continue

                feature = cols[2]
                if feature != feature_type:
                    continue

                attr_str = cols[8]
                attrs = parse_attributes(attr_str)

                if transcript_key not in attrs:
                    continue

                tid = attrs[transcript_key]

                tools = attrs.get(tools_key, "")
                if not tools:
                    continue

                tools = [t for t in tools.split(",") if t]

                for tool in tools:
                    rows.append(
                        {
                            "transcript_id": tid,
                            "tool": tool,
                            "value": 1,
                        }
                    )

        if not rows:
            raise ValueError("No valid transcript entries found.")

        df_long = pl.DataFrame(rows)
        df_wide = (
            df_long
            .pivot(
                values="value",
                index="transcript_id",
                columns="tool",
                aggregate_function="max",
            )
            .fill_null(0)
            .sort("transcript_id")
        )

        return df_wide

    # run

    df = gtf_tools_to_polars_table(
        raw_merged_gtf,
        feature_type="transcript",
        transcript_key="transcript_id",
        tools_key="tools",
    )

    df.write_csv(raw_merged_meta_tsv, separator="\t")


def plot_stat():
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv(
        raw_merged_meta_tsv,
        sep="\t"
    )

    tool_cols = [c for c in df.columns if c != "transcript_id"]
    df["n_tools"] = df[tool_cols].sum(axis=1)

    tool_counts = df[tool_cols].sum(axis=0)
    plt.figure()
    bars =plt.bar(tool_counts.index, tool_counts.values)
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{int(height)}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    plt.xlabel("Tool")
    plt.ylabel("Number of supported transcripts")
    plt.title("Number of transcripts supported by each tool")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


    max_tools = df["n_tools"].max()

    support_levels = list(range(1, max_tools + 1))
    counts = [
        (df["n_tools"] >= k).sum()
        for k in support_levels
    ]

    plt.figure()
    plt.plot(support_levels, counts, marker="o")
    plt.xlabel("At least N supporting tools")
    plt.ylabel("Number of transcripts")
    plt.title("Cumulative distribution of transcript support")
    for x, y in zip(support_levels, counts):
        plt.text(
            x,
            y,
            f"{y}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    plt.tight_layout()
    plt.show()


def prefilter_by_n_tools(n_tools=2):
    def parse_gtf_attributes(attr_str):
        attrs = {}
        for item in attr_str.strip().strip(";").split(";"):
            if not item:
                continue
            key, val = item.strip().split(" ", 1)
            attrs[key] = val.strip('"')
        return attrs

    def format_gtf_attributes(attrs):
        return "; ".join([f'{k} "{v}"' for k, v in attrs.items()]) + ";"

    df = pd.read_csv(
        raw_merged_meta_tsv,
        sep="\t"
    )
    tool_cols = [c for c in df.columns if c != "transcript_id"]
    df["n_tools"] = df[tool_cols].sum(axis=1)
    supported_tids = set(
        df.loc[df["n_tools"] >= n_tools, "transcript_id"]
    )
    print(f"Transcripts supported by ≥{n_tools} tools: {len(supported_tids)}")
    pass
    tmp_tool2_filter_gtf=f'{tool2_filter_gtf}.tmp'
    normal_chrs=[f'chr{x}' for x in [i for i in range(1,23)]+['X','Y','M']]
    with open(raw_merged_gtf) as fin, open(tmp_tool2_filter_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            if cols[0] not in normal_chrs:
                continue
            attr_str = cols[8]
            attrs = parse_gtf_attributes(attr_str)

            tid = attrs.get("transcript_id")
            if tid is None:
                continue
            if tid not in supported_tids:
                continue
            attrs["gene_id"] = tid
            cols[8] = format_gtf_attributes(attrs)
            fout.write("\t".join(cols) + "\n")

    # Convert input GFF/GTF file to a standardized GTF format using gffread
    cmd = f'{gffread_PATH} {tmp_tool2_filter_gtf} -T -o {tool2_filter_gtf} && rm {tmp_tool2_filter_gtf}'
    print(cmd)
    os.system(cmd)
    print(f"Convert to a standardized GTF format using gffread")


def split_gtf_chr():
    input_gtf = tool2_filter_gtf
    outdir = tool2_filter_gtf_chr_dir
    os.makedirs(outdir, exist_ok=True)
    valid = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    handles = {}
    def get_handle(name):
        if name not in handles:
            odir=f'{outdir}/{name}/gtf'
            os.makedirs(odir, exist_ok=True)
            handles[name] = open(f'{odir}/merged.gtf', "w")
        return handles[name]

    with open(input_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chr_ = line.split("\t", 1)[0]
            if chr_ in valid:
                get_handle(chr_).write(line)
            # else:
            #     get_handle("other").write(line)
    for h in handles.values():
        h.close()


def get_novel_gene_isoform():
    ## get GTF of isoform without SQANTI3 associated genes.


    pass


def add_gene_info():
    buildLoci_cmd='perl /lustre/home/cxue/software/buildLoci/buildLoci.pl'
    bedtools_cmd='/lustre/home/cxue/software/bedtools2/bedtools'
    # in_gtf=output_2tool_gtf
    # out_gtf=output_2tool_with_gene_gtf
    # os.makedirs(os.path.dirname(output_gff), exist_ok=True)
    # cmd=f'{bedtools_cmd} intersect -s -wao -a {in_gtf} -b {in_gtf} | {buildLoci_cmd} - > {out_gtf}'
    # print(cmd)

def merge_novel_gene_to_gtf():
    pass


def make_fa():
    cmd=f'{gffread_PATH} {tool2_filter_gtf} -g {REF_GENOME_FASTA} -w {tool2_filter_gtf[:-4]}.fa'
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':
    make_gtf_conf()
    merge_gtf_by_intron_chain()
    stat_isoform_tools_support()
    prefilter_by_n_tools(n_tools=2)
    split_gtf_chr()
    make_fa()

    pass
    # plot_stat()


