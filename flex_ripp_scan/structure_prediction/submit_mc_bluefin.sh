
#Lets pass these to bluefin
#for item in `seq 1 9`; do ls 6PQG/${item}_* | wc -l; done
#for item in `seq 10 99`; do ls 6PQG/${item}* | wc -l; done
#1AJ1_top.txt take the top peptides from here

mkdir -p temp_pdb_names

BASE=$1
for item in `seq 1 9`; do 
count=`cat ${BASE}_top.txt | grep ${BASE}/${item}_ | wc -l`; 

if [ $count -ge 1 ]; then 
#echo $count
cat ${BASE}_top.txt | grep ${BASE}/${item}_ > temp_pdb_names/${BASE}_pdbs_${item}.txt 

#split inital file
split -l 50 temp_pdb_names/${BASE}_pdbs_${item}.txt temp_pdb_names/${BASE}_pdbs_${item}_

i=1
for file in temp_pdb_names/${BASE}_pdbs_${item}_*; do
    if [[ ! "$file" =~ \.txt$ ]]; then #don't want to count files that were previously created by the script
        mv $file temp_pdb_names/${BASE}_pdbs_${item}_${i}.txt
        sbatch sample_mc_list.sh ${BASE}_mc.options temp_pdb_names/${BASE}_pdbs_${item}_${i}.txt
        ((i++))
    fi  
done

fi

#sbatch sample_mc_list.sh ${BASE}_mc.options temp_pdb_names/${BASE}_pdbs_${item}.txt

done

for item in `seq 11 99`; do 
count=`cat ${BASE}_top.txt | grep ${BASE}/${item} | wc -l`; 

if [ $count -ge 1 ]; then 
#echo $count
cat ${BASE}_top.txt | grep ${BASE}/${item} > temp_pdb_names/${BASE}_pdbs_${item}.txt 

#split inital file
split -l 50 temp_pdb_names/${BASE}_pdbs_${item}.txt temp_pdb_names/${BASE}_pdbs_${item}_

i=1
for file in temp_pdb_names/${BASE}_pdbs_${item}_*; do
    if [[ ! "$file" =~ \.txt$ ]]; then #don't want to count files that were previously created by the script
        mv $file temp_pdb_names/${BASE}_pdbs_${item}_${i}.txt
        sbatch sample_mc_list.sh ${BASE}_mc.options temp_pdb_names/${BASE}_pdbs_${item}_${i}.txt
        ((i++))
    fi  
done

fi

#sbatch sample_mc_list.sh ${BASE}_mc.options temp_pdb_names/${BASE}_pdbs_${item}.txt

done

