mkdir -p /${output}/rma6_collection

for i in $(awk -F'\t' '{print $1}' ~/eager/screening_data/mapping/sediment_streptoccocus_equinus_mapping.tsv)
    do
    i_dir=$(echo $i | rev | cut -d'_' -f2,3,4,5 | rev)
    cp ${input}/single_samples_$i_dir/results/metagenomic_classification/malt/$i.unmapped.rma6 ${output}/rma6_collection
done

rma-tabuliser -d ${output}/rma6_collection -t ${threads} -n
## clean copied rma6 files
mv ${output}/rma6_collection/count_table.tsv ${output}/count_table.tsv
rm -r ${output}/rma6_collection