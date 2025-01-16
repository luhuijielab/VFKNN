#!/bin/bash
samplename=$1
project=$2
echo ${samplename}
echo ${project}
workpath=$3
echo ${workpath}
#
mkdir -p temp/pathogen/${samplename}_kraphlan
time1=$(date "+%Y-%m-%d %H:%M:%S")

echo "For Stage 1"
########################################################################################
##  提取kraken,默认不结合kaken2
# awk -v FS='\t' -v OFS='\t' '{print $2,$3}' ./temp/kraken2/${samplename}.output | sed -E 's/\t.+\(taxid /\t/g;s/\)//g' |sed "1 i readid\tkrataxid" > temp/pathogen/${samplename}_kraphlan/${samplename}_kra.reatax # 取出readid/krataxid列


##  提取meatphlan
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' ./temp/metaphlan/${samplename}.bowtie2out.txt | sed -E 's/__.+\t.+SGB/\tSGB/g' |sed "1 i readid\tSGBid" > temp/pathogen/${samplename}_kraphlan/${samplename}_phlan_reatax.temp; # readid/SGBid列
grep "t__" ./temp/metaphlan/${samplename}_tavg_g_profiled_metagenome.txt | cut -f1,2 | sed -E "s/.+\|t__//g" | sort -k 1 > temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.tmp; #  SGBid/phlan_NCBI_taxid列 
tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.tmp | cut -f2 | grep -E -o "([0-9]+\|)$" | sed "s/|//g" > temp/pathogen/${samplename}_kraphlan/temp.tmp # 
paste temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.tmp temp/pathogen/${samplename}_kraphlan/temp.tmp | cut -f1,3 | sed "1 i SGBid\tphlantaxid" > temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.txt # SGBid/phlantaxid
rm temp/pathogen/${samplename}_kraphlan/temp.tmp temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.tmp
#
Rscript ${workpath}/${project}/combin_files.r -a temp/pathogen/${samplename}_kraphlan/${samplename}_phlan_reatax.temp -b temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.txt -c SGBid -o temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan_1.tsv > temp/pathogen/${samplename}_kraphlan/log.log 2>&1; # readid SGBid phlantaxid

## 生成${samplename}_kraphlan.tsv
tail -n+2 temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan_1.tsv | sed "1 i readid\tSGBid\ttaxid" > temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv;

## taxid list
tail -n+2 temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv | cut -f 3 | sed "/^#/d;s/NA//g" | LC_ALL=C sort -k 1 | uniq > temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst # taxid
sed -i '/^0$/d' temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst #删掉0


##  结合kra/phlan

# Rscript ${workpath}/${project}/combin_files.r -a temp/pathogen/${samplename}_kraphlan/${samplename}_kra.reatax -b temp/pathogen/${samplename}_kraphlan/${samplename}_phlan.reatax -c readid -o temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan_1.tsv > log.log 2>&1; # readid  krataxid  SGBid phlantaxid
# tail -n+2 temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan_1.tsv | awk -v FS='\t' -v OFS='\t' '{if($4=="NA") print $1,$2,$3,$4,$2; else print $1,$2,$3,$4,$4}' | sed "1 i readid\tkrataxid\tSGBid\tphlantaxid\ttaxid" > temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv; # kra综合phlan的taxid  # readid  krataxid  SGBid phlantaxid  taxid
# rm temp/pathogen/${samplename}_kraphlan/${samplename}_phlan_reatax.temp temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan_1.tsv
# ## taxid list
# awk -v FS='\t' -v OFS='\t' '{if($2=="NA") print $1,$2,$3,$4,$5}' temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv > temp/pathogen/${samplename}_kraphlan/special.txt
# tail -n+2 temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv | cut -f 5 | sed "/^#/d;s/NA//g" | LC_ALL=C sort -k 1 | uniq > temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst # taxid
# sed -i '/^0$/d' temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst #删掉0

# ## taxonomy 依据 kra taxonkit
# # kraken 的rank_id_name文件  用taxonkit
# taxonkit lineage temp/pathogen/${samplename}_kraphlan/${samplename}_krataxid.lst -r | awk '$2!=""' > temp/pathogen/${samplename}_kraphlan/${samplename}temp
# cat temp/pathogen/${samplename}_kraphlan/${samplename}temp | taxonkit reformat -P -t | cut -f1,3,4,5 | sed "s/;/\|/g" | sed "1 i krataxid\tkra_NCBI_taxonomy_minlv\tkra_NCBI_taxonomy\tkra_NCBI_taxid" > temp/pathogen/${samplename}_kraphlan/${samplename}krataxidlineage.tsv
# # metaphlan 的rank_id_name文件  用taxonkit
# tail -n+2 temp/pathogen/${samplename}_kraphlan/${samplename}_SGB_taxid.txt | cut -f2 | taxonkit lineage -r | awk '$2!=""' > temp/pathogen/${samplename}_kraphlan/${samplename}temp
# cat temp/pathogen/${samplename}_kraphlan/${samplename}temp | taxonkit reformat -P -t | cut -f1,3,4,5 | sed "s/;/\|/g" | sed -E 's/$/|/g' | sed "1 i phlantaxid\tphlan_NCBI_taxonomy_minlv\tphlan_NCBI_taxonomy\tphlan_NCBI_taxid" > temp/pathogen/${samplename}_kraphlan/${samplename}phlantaxidlineage.tsv
# # taxid 的rank_id_name文件  用taxonkit
# taxonkit lineage temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst -r | awk '$2!=""' > temp/pathogen/${samplename}_kraphlan/${samplename}temp
# cat temp/pathogen/${samplename}_kraphlan/${samplename}temp | taxonkit reformat -P -t | cut -f1,3,4,5 | sed "s/;/\|/g" | sed "1 i taxid\tNCBI_taxonomy_minlv\tNCBI_taxonomy\tNCBI_taxid" > temp/pathogen/${samplename}_kraphlan/${samplename}taxidlineage.tsv

# rm temp/pathogen/${samplename}_kraphlan/${samplename}temp 
#
time2=$(date "+%Y-%m-%d %H:%M:%S")
echo "For Stage 2"
########################################################################################
## 依据taxid提取序列
# 获取总reads数
nreads=$(sed -n 'x;$p'  ./temp/metaphlan/${samplename}.bowtie2out.txt | cut -f2)
all=$(awk "BEGIN{print $(echo "${nreads}*0.00001" | bc -l) + 0}")
# echo "nreads: ${nreads}"

# 按照taxid提取id
mkdir -p temp/pathogen/${samplename}_kraphlan/extractseqid temp/pathogen/${samplename}_kraphlan/extractseq
tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst | cut -f1 |rush -j 20 "grep -w '{1}$' temp/pathogen/${samplename}_kraphlan/${samplename}_kraphlan.tsv | cut -f1 > temp/pathogen/${samplename}_kraphlan/extractseqid/{1}_id.txt"

# 重构taxid.lst
cp temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst temp/pathogen/${samplename}_kraphlan/${samplename}_taxid_backup.lst
for i in `tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst | cut -f1 `;do
     num=$(awk "BEGIN{print $(sed -n '$=' temp/pathogen/${samplename}_kraphlan/extractseqid/${i}_id.txt)*2 + 0}")
     # declare -i num
     if [ `echo "$num < $all"|bc` -eq 1 ]
     then
          # print $num
          sed -i "/^${i}$/d" temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst
     fi
done
# 提取seq
tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst | cut -f1 |rush -j 24 "seqtk subseq temp/concat/${samplename}.fq temp/pathogen/${samplename}_kraphlan/extractseqid/{1}_id.txt | grep -A1 '@' | sed -e '/^--/d;s/@/>/g' >> temp/pathogen/${samplename}_kraphlan/extractseq/{1}_out.fa"

time3=$(date "+%Y-%m-%d %H:%M:%S")
echo "For Stage 3"
#######################################################################################
## reads比对vfdb
mkdir -p temp/pathogen/${samplename}_kraphlan/etractseqblastvf/
tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_taxid.lst | cut -f1 | rush -j 20 'blastn -query temp/pathogen/{samplename}_kraphlan/extractseq/{1}_out.fa -db /media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/vfdb/VFDB_B_nt -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -num_threads 6 -out temp/pathogen/{samplename}_kraphlan/etractseqblastvf/{1}_vfblastout.txt' -v samplename=${samplename}

##  计算TPM、cell信息
mkdir -p temp/pathogen/${samplename}_kraphlan/16S_cells_metadata/ temp/pathogen/${samplename}_kraphlan/vftpm/ temp/pathogen/${samplename}_kraphlan/temp/
cd ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/etractseqblastvf
find . -name "*" -type f -size +1c | sed 's/_vfblastout.txt//g' | sed 's/.\///g' > ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/${samplename}_vftaxid.txt
cd ${workpath}/${project}/

# 借鉴自ARG-OAP
for i in `tail -n+1 temp/pathogen/${samplename}_kraphlan/${samplename}_vftaxid.txt|cut -f 1`;do
     # echo "for:${i}"
     python /media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/scripts/py/16S_cells_abundance_stage1.py -i ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/extractseq/${i}_out.fa -o ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/ --sampname ${i} -t 10  #step1   16S_cells_metadata.tsv
     python /media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC/scripts/py/16S_cells_abundance_stage2.py -i ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/etractseqblastvf/${i}_vfblastout.txt -o ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/ --sampname ${i} -t 10 -m ${workpath}/${project}/temp/pathogen/${samplename}_kraphlan/16S_cells_metadata/${i}_16S_cells_metadata.tsv --dbtype nucl  #step2  _vf_tpm.tsv
     # echo "done\n"
done
time4=$(date "+%Y-%m-%d %H:%M:%S")
echo $time1
echo $time2
echo $time3
echo $time4
## 移动文件
mkdir -p temp/pathogen/${samplename}
mv temp/pathogen/${samplename}_kraphlan/* temp/pathogen/${samplename}