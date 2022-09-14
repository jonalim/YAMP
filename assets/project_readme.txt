DIRECTORY STRUCTURE

SAMPLE DIRECTORIES

ILLUMINA_DATA - contains raw reads and Megablast results.

ANALYSIS
- 01_fastqc - Quality reports on raw reads.
- 02_dedup_log.txt - Log file from Clumpify.
- 03_syndecontam_log.txt - Log file from BBDuk.
- 04_trim+hostremove_log.txt - Log file from Kneaddata.
- 05_QCd.fq.gz - Quality-controlled reads.
- 05_qc_stats.txt - Summary of stats extracted from tool logs.
- 05_fastqc - Quality reports on QC'd reads.

- 06_*.metaphlan_bugs_list.tsv - Sample relative abundances calculated by 
	Metaphlan. (https://huttenhower.sph.harvard.edu/metaphlan)
- 06_*.biom - Sample relative abundances calculated by Metaphlan.
	(https://huttenhower.sph.harvard.edu/metaphlan)
- 06_*.alpha-diversity.tsv - Alpha-diversity metrics calculated using QIIME2.
- 07_*.HUMAnN.log - HUMAnN log file.
- 07_*_genefamilies.tsv - Abundance-RPKs of gene families.
	(https://github.com/biobakery/humann#1-gene-families-file)
- 07_*_pathabundance.tsv - Sample pathway abundance in RPK.
	(https://github.com/biobakery/humann#2-pathway-abundance-file)
- 07_*_pathcoverage.tsv - Sample pathway coverage.
	(https://github.com/biobakery/humann#3-pathway-coverage-file)
- 07_profiling_stats.txt - Summary of stats extracted from profiling output.
- 07_humann_intermediate - Files listing the gene family assignments of individual 
	reads. Not included unless requested for specific samples.
	(https://github.com/biobakery/humann#4-intermediate-temp-output-files)

PROJECT DIRECTORY

08_merged_metaphlan_abundance_table.txt - percent relative abundance, at all 
	taxonomic levels.
	(https://github.com/biobakery/biobakery/wiki/metaphlan3#merge-outputs)
08_merged_metaphlan_abundance_table.species.txt - percent relative abundance of
	all identified species.
08_all_genefamilies.rpk.tsv - gene family counts from all samples (reads per
	kilobase). (https://github.com/biobakery/humann#1-gene-families-file)
08_all_genefamilies.copm.unstratified.tsv - Unstratified gene family counts 
	(copies per million).
08_all_genefamilies.copm.stratified.tsv - Taxonomically stratified gene family
	abundances (copies per million).
08_all_pathabundance.rpk.tsv - pathway abundances from all samples.
	(https://github.com/biobakery/humann#2-pathway-abundance-file)
08_all_pathabundance.copm.unstratified.tsv - Unstratified pathway abundances
	(copies per million).
08_all_pathabundance.copm.stratified.tsv - Taxonomically stratified pathway
	abundances (copies per million).
08_all_pathcoverage.tsv - pathway coverage scores from all samples.
	(https://github.com/biobakery/humann#3-pathway-coverage-file)
08_all_pathcoverage.unstratified.tsv - Unstratified pathway coverage values.
08_all_pathcoverage.stratified.tsv - Taxonomically stratified pathway coverage
	values.

09_metaphlan.species.top50.hclust2.png - Heatmap of top 50 species, sorted by
	sum of relative abundance.
09_all_genefamilies.copm.stratified.top50.hclust2.png - Heatmap of top 50
	taxonomically stratified gene families.
09_all_pathabundance.copm.unstratified.top50.hclust2.png - Abundance heatmap of 
	top 50 most abundant pathways.
09_all_pathcoverage.unstratified.no-unintegrated.top50.hclust2.png - Heatmap of coverage values
	of top 50 most covered pathways. UNMAPPED and UNINTEGRATED values have been
	omitted as they are meaningless. 

10_qc_stats.txt - Project-level table of sample read throughput.
10_profiling_stats.txt - Project-level table summarizing taxonomic and 
	functional profiling.
10_*_multiqc_report.html - MultiQC report summarizing throughput and profiling.
10_*_multiqc_report_data/ - Additional data curated by MultiQC.


Where indicated that data are in copies per million, depth-normalization was 
applied as described: https://github.com/biobakery/humann#humann_renorm_table

