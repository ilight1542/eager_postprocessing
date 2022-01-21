/ptmp/iclight/eager/mapping/dairy_screening/yersinia_pestis/results/deduplication
conda activate basic_genomics

find . -name "*.bam" > files

for i in $(seq 1 $(wc -l <files))
    do
    current_file=$(awk -v i="$i" 'NR == i {print}' files)
    echo $current_file >> output_read_counts.txt
    samtools idxstats $current_file | cut -f 1,3 >> output_read_counts.txt
done

for i in $(seq 1 $(wc -l <files))
    do
    current_file=$(awk -v i="$i" 'NR == i {print}' files)
    if [[ $i == 1 ]]
    then
        echo $current_file | sed 's#/#\t#g' | cut -f 3 | awk '{print "\t\t"$0}' >> output_read_counts.txt.${i}
        samtools idxstats $current_file | cut -f 1,3  >> output_read_counts.txt.${i}
    else
        echo $current_file | sed 's#/#\t#g' | cut -f 3 >> output_read_counts.txt.${i}
        samtools idxstats $current_file | cut -f 3  >> output_read_counts.txt.${i}
    fi
done

paste -d '\t' output_read_counts.txt.* | head -n+5 > ${output_name}.mapped_reads