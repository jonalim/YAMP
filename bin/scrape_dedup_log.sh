# MultiQC doesn't have a module for clumpify yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

# Log some figures about sequences passing de-duplication
totR=$(grep "Reads In:" $1 | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
remR=$(grep "Duplicates Found:" $1 | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
survivedR=$(($totR - $remR))
percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
percentage=`printf "%.2f" $percentage`
time=$(grep "Total time:" $1 | cut -d: -f 2 | cut -f 2 | sed 's/s\./s/g')

samplename=$(head -n 1 $1 | grep -o "\sin=[^\S\.]*" | sed 's\[[:blank:]]in=\\' )

#printf "\tinput\tduplicated\tsurviving\ttime\n"
printf "%s\t%s\t%s\t%s\t%s\n" "${samplename}" "${totR}" "${remR}" "${survivedR}" "${time}"

# Dump to YAML
# echo "id: 'deduplication'"
# echo "section_name: 'YAMP Deduplication'" 
# echo "section_href: 'https://github.com/alesssia/yamp'" 
# echo "plot_type: 'html'" 
# echo "description: 'This information is collected at run time from the software output.'" 
# echo "data: |" 
# echo "    <dl class="dl-horizontal">" 
# echo "        <dt>Input:</dt><dd>"$totR"</dd>" 
# echo "        <dt>Duplicated Found:</dt><dd>"$remR"</dd>" 
# echo "        <dt>Surviving:</dt><dd>"$survivedR" ("$percentage"% of input)</dd>" 
# echo "        <dt>Total time:</dt><dd>"$time"</dd>" 
# echo "    </dl>" 


