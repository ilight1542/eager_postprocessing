## run inside the directory containing all the loop iterations (eg subdirs of 1,2,3,..100,101,..)

function malt_extract {
    output=$1
    input_dir=$2

    mkdir -p ${output}/ancient/coverage/
    mkdir -p ${output}/ancient/damageMismatch/
    mkdir -p ${output}/ancient/editDistance/
    mkdir -p ${output}/ancient/filterInformation/
    mkdir -p ${output}/ancient/readDist/
    mkdir -p ${output}/ancient/percentIdentity/
    mkdir -p ${output}/default/coverage/
    mkdir -p ${output}/default/damageMismatch/
    mkdir -p ${output}/default/editDistance/
    mkdir -p ${output}/default/filterInformation/
    mkdir -p ${output}/default/readDist/
    mkdir -p ${output}/default/percentIdentity/

    anc_runsum=ancient/RunSummary.txt
    def_runsum=default/RunSummary.txt
    anc_totalcount=ancient/TotalCount.txt
    def_totalcount=default/TotalCount.txt

    # get first col for starting files and move unique non-shared files 
    for i in ${input_dir}/*/results/maltextract/results/
        do
        cp ${i}/ancient/coverage/*.txt ${output}/ancient/coverage/
        cp ${i}/ancient/damageMismatch/*.txt ${output}/ancient/damageMismatch/
        cp ${i}/ancient/editDistance/*.txt ${output}/ancient/editDistance/
        cp ${i}/ancient/filterInformation/*.txt ${output}/ancient/filterInformation/
        cp ${i}/ancient/percentIdentity/*.txt ${output}/ancient/percentIdentity/
        cp ${i}/ancient/readDist/*.txt ${output}/ancient/readDist/
        cp ${i}/default/coverage/*.txt ${output}/default/coverage/
        cp ${i}/default/damageMismatch/*.txt ${output}/default/damageMismatch/
        cp ${i}/default/editDistance/*.txt ${output}/default/editDistance/
        cp ${i}/default/filterInformation/*.txt ${output}/default/filterInformation/
        cp ${i}/default/percentIdentity/*.txt ${output}/default/percentIdentity/
        cp ${i}/default/readDist/*.txt ${output}/default/readDist/
    done

    python3 ~/bin/merge_summary_maltextract.py -i ${input_dir}/*/results/maltextract/results/${anc_runsum} -o RunSummary.txt
    mv RunSummary.txt ${output}/ancient/
    python3 ~/bin/merge_summary_maltextract.py -i ${input_dir}/*/results/maltextract/results/${def_runsum} -o RunSummary.txt
    mv RunSummary.txt ${output}/default/
    python3 ~/bin/merge_summary_maltextract.py -i ${input_dir}/*/results/maltextract/results/${anc_totalcount} -o TotalCount.txt
    mv TotalCount.txt ${output}/ancient/
    python3 ~/bin/merge_summary_maltextract.py -i ${input_dir}/*/results/maltextract/results/${def_totalcount} -o TotalCount.txt
    mv TotalCount.txt ${output}/default/

    ## running malt extract
    conda activate basic_genomics
    if [ $# -eq 3 ];
	then
        postprocessing.AMPS.r -r ${output} -n ${3} -t 8 -m def_anc
    else
        postprocessing.AMPS.r -r ${output} -n ~/bin/HOPS/Resources/default_list.txt -t 8 -m def_anc
    fi
}
