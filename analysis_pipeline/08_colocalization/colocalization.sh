
main(){
    finemapping_for_GWAS
    GWAS_QTL_pair
    finemapping_for_QTL
    coloc_for_GWAS_QTL
}


finemapping_for_GWAS(){

    mkdir -p slurm/fmGWAS/log
    # grep -v "Height" input/EAS_addGWAS_all_sentinal.356.txt | sed "s/\"//g" | tail -n+2 | while read line;
    cat input/d.txt | while read line
    do
        GWAS_name=`echo ${line} | cut -f 1 -d " "`
        CHR_name=`echo ${line} | cut -f 2 -d " "`
        sentinel_pos=`echo ${line} | cut -f 3 -d " "`
        sentinel_rsid=`echo ${line} | cut -f 5 -d " "`
        task=${GWAS_name}_${sentinel_rsid}
        echo $task
        echo "#!/bin/bash

#SBATCH --job-name=$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=slurm/fmGWAS/log/${task}.err
#SBATCH --output=slurm/fmGWAS/log/${task}.out

Rscript scripts/a1.finemapping_for_GWAS.R ${GWAS_name} ${CHR_name} ${sentinel_pos} ${sentinel_rsid}

" > slurm/fmGWAS/${task}.slurm
        sbatch slurm/fmGWAS/${task}.slurm
    
        process=$(squeue|wc -l)
        while((process >= 200))
        do
            echo "Current jobs' count is larger than 200"
            echo "Wait for another 60 s"
            sleep 60
            process=$(squeue|wc -l)
        done
    done
}


GWAS_QTL_pair(){

    mkdir -p slurm/gwasqtl/log
    
    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"
    for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL

    # TISSUE="Adipose_Visceral_Omentum,Adrenal_Gland,Liver,Muscle_Skeletal,Pancreas,Skin_Not_Sun_Exposed_Suprapubic,Spleen,Whole_Blood"
    # for xQTL_type in GTExv8_eQTL  GTExv8_juQTL
    
    # TISSUE="Whole_Blood"
    # for xQTL_type in JCTF_eQTL JCTF_juQTL    

    # TISSUE="LCL"
    # for xQTL_type in MAGE_eQTL MAGE_juQTL
    do
        for tissue in ${TISSUE//,/ }
        do
            task=${xQTL_type}_${tissue}
            echo $task
    echo "#!/bin/bash

#SBATCH --job-name=$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=slurm/gwasqtl/log/${task}.err
#SBATCH --output=slurm/gwasqtl/log/${task}.out

Rscript scripts/a2.finemapping_prepare_genes.R ${xQTL_type} ${tissue}

" > slurm/gwasqtl/${task}.slurm
            
            sbatch slurm/gwasqtl/${task}.slurm

            process=$(squeue|wc -l)
            while((process >= 100))
            do
                echo "Current jobs' count is larger than 10"
                echo "Wait for another 60 s"
                sleep 60
                process=$(squeue|wc -l)
            done
        done
    done
}


# cat input/finemapping_QTL/SNV_eQTL/gwas_qtl_pair_* | grep -v "xQTL_type" | cut -f 1-3 | sort -u > input/list_finemapping_gene.txt

finemapping_for_QTL(){

    tasknum=task11
    mkdir -p slurm/fmqtl/${tasknum}/log slurm/fmqtl/${tasknum}/task

    split -l 200 -d input/list_finemapping_gene.txt slurm/fmqtl/${tasknum}/task/task_

    for f in `ls slurm/fmqtl/${tasknum}/task/task_* `
    do
        task=`basename $f`
        echo $task
        echo "#!/bin/bash

#SBATCH --job-name=${tasknum}_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --error=slurm/fmqtl/${tasknum}/log/${task}.err
#SBATCH --output=slurm/fmqtl/${tasknum}/log/${task}.out

cat $f | while read line
do
    QTLTYPE=\`echo \$line | awk '{print \$1}'\`
    TISSUE=\`echo \$line | awk '{print \$2}'\`
    EGENE=\`echo \$line | awk '{print \$3}'\`
    
    Rscript scripts/a3.finemapping_for_QTL.R \${QTLTYPE} \${TISSUE} \${EGENE}

done

"   >   slurm/fmqtl/${tasknum}/${task}.slurm

        sbatch slurm/fmqtl/${tasknum}/${task}.slurm
    done
}


# cat input/finemapping_QTL/SNV_eQTL/gwas_qtl_pair_* | grep -v "xQTL_type" > input/list_for_coloc.txt
# wc -l input/list_for_coloc.txt

coloc_for_GWAS_QTL(){

    tasknum=task5
    mkdir -p slurm/coloc/${tasknum}/log slurm/coloc/${tasknum}/task slurm/coloc/${tasknum}/output

    split -l 3200 input/list_for_coloc.txt slurm/coloc/${tasknum}/task/task_

    for f in `ls slurm/coloc/${tasknum}/task/task_* `
    do
        task=`basename $f`
        echo $task
        echo "#!/bin/bash

#SBATCH --job-name=coloc_${tasknum}_$task
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --error=slurm/coloc/${tasknum}/log/${task}.err
#SBATCH --output=slurm/coloc/${tasknum}/log/${task}.out

Rscript scripts/a4.coloc.R $f 32 slurm/coloc/${tasknum}/output/${task}.out

"   >   slurm/coloc/${tasknum}/${task}.slurm

        sbatch slurm/coloc/${tasknum}/${task}.slurm

    done
}


main
