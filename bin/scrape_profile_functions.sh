# MultiQC doesn't have a module for human yet. As a consequence, I
# had to create a YAML file with all the info I need via a bash script

totR=$(gunzip < $3 | grep "^@" | wc -l)
tot_species_prescreening=$(grep "Total species selected from prescreen:" $2 | cut -d: -f 2 | sed 's/ //g')
selected_species_explain=$(grep "Selected species explain" $2 | cut -d" " -f 4- | sed "s/% of predicted community composition/% of QC'd reads/g" )
tot_reads_aligned_nucleotide=$(grep -zo "Total bugs from nucleotide alignment.*Total gene families from nucleotide alignment" $2 | grep -Eoa "[0-9]+ hits" | grep -Eo "[0-9]+" | paste -sd+ | bc | sed 's/ //g')
unaligned_reads_nucleotide=$(grep "Unaligned reads after nucleotide alignment:" $2 | cut -d: -f 2 | sed 's/ //g')
tot_reads_aligned_both=$(grep -zo "Total bugs after translated alignment.*Total gene families after translated alignment" $2 | grep -Eoa "[0-9]+ hits" | grep -Eo "[0-9]+" | paste -sd+ | bc | sed 's/ //g')
unaligned_reads_translated=$(grep "Unaligned reads after translated alignment:" $2 | cut -d: -f 2 | sed 's/ //g')
tot_gene_family=$(grep "Total gene families after translated alignment:" $2 | cut -d: -f 2 | sed 's/ //g')

printf "\tinput\tspecies\texplain\taligned_nucleotide\tunaligned_nucleotide_percent\taligned_both\tunaligned_both\tgene_families\n"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $1 ${totR} ${tot_species_prescreening} ${selected_species_explain} ${tot_reads_aligned_nucleotide} ${unaligned_reads_nucleotide} ${tot_reads_aligned_both} ${unaligned_reads_translated} ${tot_gene_family}
# # Dump to YAML (header)
# echo "id: 'humann'"
# echo "section_name: 'HUMAnN'" 
# echo "section_href: 'https://github.com/alesssia/yamp'" 
# echo "plot_type: 'html'" 
# echo "description: 'This information is collected at run time from the software output.'" 
# echo "data: |" 
# echo "    <dl class="dl-horizontal">" 
# echo  "        <dt>Input after QC</dt><dd>"${totR}"</dd>" 
# echo  "        <dt>Selected from prescreen</dt><dd>"${tot_species_prescreeing}" species</dd>" 
# echo  "        <dt>Selected species explain</dt><dd>"${selected_species_explain}"</dd>" 
# echo  "        <dt>Unaligned</dt><dd>"${unaligned_reads_nucleotide}" QC'd reads unaligned after nucleotide alignment</dd>" 
# echo  "        <dt>Unaligned</dt><dd>"${unaligned_reads_translated}" QC'd reads unaligned after translated alignment</dd>" 
# echo  "        <dt>Total gene families</dt><dd>"${tot_gene_family}" (after translated alignment)</dd>" 
# echo  "        <dt>Complete log</dt><dd>"$1"_HUMAnN.log</dd>" 
# echo "    </dl>" 

