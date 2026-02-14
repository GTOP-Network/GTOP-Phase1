
QTL_type=SNV_eQTL

dir_permutation=output/${QTL_type}/permutation
dir_nominal=output/${QTL_type}/nominal
dir_output_Strong=output/${QTL_type}/Strong
dir_output_Random=output/${QTL_type}/Random

./bin/MashR/1-1Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_output_Strong}
./bin/MashR/1-2Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_nominal} ${dir_output_Strong}
./bin/MashR/1-3Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_output_Strong}
./bin/MashR/2-1Prepare-MashR-RandomPairs.sh ${dir_nominal} ${dir_output_Random}

mkdir -p slurm/log
rm slurm/${QTL_type}_*
perm_files=(`find ${dir_permutation} -name "*.txt"`)
echo ${#perm_files[@]}
nFile=${#perm_files[@]}
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/.txt/})

for i in `seq 0 $[nFile-1]`
do
    name=${NAMEs[i]}
    nominal_files=(`find ${dir_nominal} -name "${name}.txt.gz"`)
    echo "/bin/sh

#SBATCH --job-name=${QTL_type}_${name}
#SBATCH --partition=Compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --error=slurm/log/${QTL_type}_${name}.err
#SBATCH --output=slurm/log/${QTL_type}_${name}.out

# if [[ ! -f ${dir_output_Random}/${name}.nominal_pairs.extracted_pairs.txt.gz ]];then
echo ${dir_nominal}/${name}.txt.gz > ${dir_output_Random}/${name}.nominal_files.txt
python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/extract_pairs_tjy.py ${dir_output_Random}/${name}.nominal_files.txt ${dir_output_Random}/nominal_pairs.combined_signifpairs.txt.gz ${name}.nominal_pairs -o ${dir_output_Random}
# rm -f ${dir_output_Random}/${name}.nominal_files.txt

"   > slurm/${QTL_type}_${name}.slurm
    sed -i "s/\/bin\/sh/\#\!\/bin\/sh/" slurm/${QTL_type}_${name}.slurm
done

ls slurm/${QTL_type}_*.slurm | parallel -j 4 bash {}

ls output/${QTL_type}/Random/*.nominal_pairs.extracted_pairs.txt.gz | wc -l

./bin/MashR/2-3Prepare-MashR-RandomPairs.sh ${dir_output_Random} 

subset_size=1000000
Rscript ./bin/MashR/MashR-random_subset.R ${dir_output_Random} ${subset_size}

ll -h output/${QTL_type}/Random/*.MashR_input.txt.gz
ll -h output/${QTL_type}/Random/MashR.random_subset_1000000.RDS

Rscript ./bin/MashR/run_MashR.R ${dir_output_Strong}/strong_pairs.MashR_input.txt.gz ${dir_output_Random}/MashR.random_subset_${subset_size}.RDS 0 ./output/${QTL_type}/top_pairs


## finemapping result
QTL_type=SNV_eQTL

dir_finemapping=output/${QTL_type}/finemapping
dir_nominal=output/${QTL_type}/nominal
dir_output_finemapping=output/${QTL_type}/Strong2

./bin/MashR/3-1Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_output_finemapping}
./bin/MashR/3-2Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_nominal} ${dir_output_finemapping}
./bin/MashR/3-3Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_output_finemapping}
