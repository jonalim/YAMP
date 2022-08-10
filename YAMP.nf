#!/usr/bin/env nextflow

/**
MD Genomics Metagenomic Analysis Pipeline
Copyright (C) 2022 Jonathan Lim.

This program incorporates a modified version of Yet Another Metagenomic Pipeline (YAMP).
Copyright (C) 2017-2021    Dr Alessia Visconti

This script is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this script. If not, see <http://www.gnu.org/licenses/>.
*/

def versionMessage()
{
    log.info"""

    MD Genomics Metagenomic Analysis Pipeline - Version: ${workflow.manifest.version}
    """.stripIndent()
}

def helpMessage()
{
    log.info"""

MD Genomics Metagenomic Analysis Pipeline - Version: ${workflow.manifest.version}
Copyright (C) 2022 Jonathan Lim.

This program incorporates a modified version of Yet Another Metagenomic Pipeline (YAMP).
Copyright (C) 2017-2021    Dr Alessia Visconti

This pipeline is distributed in the hope that it will be useful
but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.

  Usage:
  nextflow run YAMP.nf --reads1 R1 --reads2 R2 --prefix prefix --outdir path [options]

  Mandatory arguments:
    --reads1   R1      Forward (if paired-end) OR all reads (if single-end) file path
    [--reads2] R2      Reverse reads file path (only if paired-end library layout)
    --prefix   prefix  Prefix used to name the result files
    --outdir   path    Output directory (will be outdir/prefix/)

  Profiles:
    --profile base,conda,sge    Profile for all MD Genomics analysis.
	--profile test,conda,sge	Profile for validation runs

  Other options:
  Parameters for removing synthetic contaminants:
    --qin                 <33|64> Input quality offset
    --artefacts           path    FASTA file with artefacts
    --phix174ill          path    FASTA file with phix174_ill

  Parameters for adapter/quality trimming:
    --phred               value   regions with average quality BELOW this will be trimmed
    --minlength           value   reads shorter than this after trimming will be discarded
    --adapters            path    FASTA file with adapters

  Parameters for decontamination:
    --foreign_genome_ref  path    folder for for contaminant (pan)genome (directory of bowtie2 indexes)

  MetaPhlAn parameters for taxa profiling:
    --metaphlan_databases path    folder for the MetaPhlAn database
    --bt2options          value   BowTie2 options

  HUMANn parameters for functional profiling:
    --chocophlan          path    folder for the ChocoPhlAn database
    --uniref              path      folder for the UniRef database
  
  --only_samples		  path		Changes the way in which this pipeline finds samples to process.
									The pipeline will only process the samples listed in the file
									provided as the value of this option. This file should be a TSV with
									each sample name at the start of each line, without a header.

  --combine_multiqc					Gather existing summary stats in the project directory and combine them with new
									data files to create a unified MultiQC report. Any new data for samples that
									were previously processed will override old data. NB: Running this pipeline
									twice in the same directory will ALWAYS overwrite any existing files, including the
									summary_stats.txt file. Before running the pipeline in a directory containing
									existing data, back up the summary_stats.txt file.

MD Genomics Metagenomic Analysis Pipeline supports FASTQ and compressed FASTQ files.
"""
}

/**
Prints version when asked for
*/
if (params.version) {
    versionMessage()
    exit 0
}

/**
Prints help when asked for
*/

if (params.help) {
    helpMessage()
    exit 0
}

/**
STEP 0.

Checks input parameters and (if it does not exists) creates the directory
where the results will be stored (aka working directory).
Initialises the log file.

The working directory is named after the prefix and located in the outdir
folder. The log file, that will save summary statistics, execution time,
and warnings generated during the pipeline execution, will be saved in the
working directory as "prefix.log".
*/

//Checking user-defined parameters
if (params.qin != 33 && params.qin != 64) {
    exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)"
}

//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
// if (!params.skip_preprocess && !params.singleEnd && (params.reads2 == "null") ) {
//     exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
// }

// Given --indir, this script knows how to find the inputs for any entrypoint
if (params.indir == "null") {
    exit 1, "Please set the --indir parameter"
}

//Creates working dir
// workingpath = params.outdir + "/" + params.prefix
// workingdir = file(workingpath)
// if( !workingdir.exists() ) {
//     if( !workingdir.mkdirs() )     {
//         exit 1, "Cannot create working directory: $workingpath"
//     }
// }

// Header log info
log.info """---------------------------------------------
MD Genomics Metagenomic Analysis Pipeline
---------------------------------------------

Analysis introspection:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date()
summary['Pipeline Name'] = 'MD Genomics Metagenomic Analysis Pipeline'
summary['Pipeline Version'] = workflow.manifest.version

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume

summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version")
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name

summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container

if (!params.skip_preprocess)
{
    if (workflow.containerEngine == 'singularity') {
        summary['FastQC'] = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
    } else if (workflow.containerEngine == 'docker') {
        summary['FastQC'] = "quay.io/biocontainers/fastqc:0.11.9--0"
    } else {
        summary['FastQC'] = "No container information"
    }
}

if (! params.skip_profile_taxa || ! params.skip_profile_function)
{
    if (workflow.containerEngine == 'singularity') {
        summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
        summary['qiime'] = "qiime2/core:2020.8"
    } else if (workflow.containerEngine == 'docker') {
        summary['biobakery'] = "biobakery/workflows:3.0.0.a.6.metaphlanv3.0.7"
        summary['qiime'] = "qiime2/core:2020.8"
    } else {
        summary['biobakery'] = "No container information"
        summary['qiime'] = "No container information"
    }
}

if (workflow.containerEngine == 'singularity') {
    summary['MultiQC'] = "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
} else if (workflow.containerEngine == 'docker') {
    summary['MultiQC'] = "quay.io/biocontainers/multiqc:1.9--py_1"
} else {
    summary['MultiQC'] = "No container information"
}

if(workflow.profile == 'awsbatch'){
    summary['AWS Region'] = params.awsregion
    summary['AWS Queue'] = params.awsqueue
}

//General
summary['Running parameters'] = workflow.commandLine
summary['Layout'] = 'Paired-End'
// summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'

if (!params.skip_preprocess)
{
    // if (!params.singleEnd)
    // {
    //     summary['QC mate pairs separately'] = params.qc_matepairs
    // }
    summary['Performing de-duplication'] = params.dedup
    summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'

    //remove_synthetic_contaminants
    summary['Synthetic contaminants:'] = ""
    summary['Artefacts'] = params.artefacts
    summary['Phix174ill'] = params.phix174ill
    summary['K-mer size'] = "31"

    //Trimming
    summary['Adapter/quality trimming:'] = ""
    summary['Adapters'] = params.adapters
    summary['Trimmomatic parameters:'] = "ILLUMINACLIP:[adapters]:2:30:10:8:True SLIDINGWINDOW:4:$params.phred MINLEN:$params.minlength"

    //Decontamination
    summary['Host decontamination:'] = ""
    if (params.foreign_genome_ref != "") {
        summary['Contaminant (pan)genome'] = params.foreign_genome_ref + " (indexed)"
    } else if (    params.foreign_genome_ref == "") {
        summary['Contaminant (pan)genome'] = params.foreign_genome
    }
    summary['Decontamination parameters'] = "Kneaddata defaults"
}

if (!params.skip_profile_taxa || !params.skip_profile_function)
{
    //BowTie2 databases for metaphlan
    summary['MetaPhlAn parameters:'] = ""
    summary['MetaPhlAn database'] = params.metaphlan_databases
    summary['Bowtie2 options'] = params.bt2options

    // ChocoPhlAn and UniRef databases
    summary['HUMAnN parameters:'] = ""
    summary['Chocophlan database'] = params.chocophlan
    summary['Uniref database'] = params.uniref
}

//Folders
summary['Folders:'] = ""
summary['Input directory'] = params.indir
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Launching dir'] = workflow.launchDir

log.info summary.collect { k,v -> "${k.padRight(27)}: $v" }.join("\n")
log.info ""

// def sample_pubdir(name) {
//     folder = new File(params.outdir + "/" + name + "/ANALYSIS/")
//     if( !folder.exists() ) {
//         folder.mkdirs()
//     }
//     return folder.getPath()
// }

/**
    Prepare workflow introspection

    This process adds the workflow introspection (also printed at runtime) in the logs
    This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'workflow-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'Pipeline Workflow Summary'
    section_href: 'https://github.com/jonalim/YAMP'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd>$v</dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/**
    Gets software version.

    This process ensures that software version are included in the logs.
*/

process get_software_versions {
    label 'biobakery'

    //Starting the biobakery container. I need to run metaphlan and Humann to get
    //their version number (due to the fact that they live in the same container)
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

    output:
    file("software_versions_mqc.yaml") into software_versions_yaml

    script:
    //I am using a multi-containers scenarios, supporting docker and singularity
    //with the software at a specific version (the same for all platforms). Therefore, I
    //will simply parse the version from there. Perhaps overkill, but who cares?
    //This is not true for the biobakery suite (metaphlan/humann) which extract the
    //information at runtime from the actual commands (see comment above)
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt

    echo $params.docker_container_fastqc | cut -d: -f 2 > v_fastqc.txt
    echo $params.docker_container_bbmap | cut -d: -f 2 > v_bbmap.txt

    metaphlan --version > v_metaphlan.txt
    humann --version > v_humann.txt
    echo $params.docker_container_qiime2 | cut -d: -f 2 > v_qiime.txt

    echo $params.docker_container_multiqc | cut -d: -f 2 > v_multiqc.txt

    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

// software_versions_yaml = software_versions_yaml.first()

// Defines channels for foreign_genome file

// index_foreign_genome_log = Channel.empty()
// foreign_genome = Channel.empty()

if (params.foreign_genome_ref == "") {
	// FIXME write code to download the bowtie2 index to a work directory
	//Stage boilerplate log when the contaminant (pan)genome is indexed
	// index_foreign_genome_log = Channel.from(file("$baseDir/assets/foreign_genome_indexing_mqc.yaml"))
	// foreign_genome = file( "${params.foreign_genome}", type: "file", checkIfExists: true )
} else {
	//When the indexed contaminant (pan)genome is already available, its path should be pushed in
	// the correct channel
	ref_foreign_genome = Channel.value(file(params.foreign_genome_ref, checkIfExists: true ))
}


if(params.only_samples) {
    // Load samples from this manifest into a glob pattern
    glob_pattern = "${params.indir}/{"
    file(params.only_samples, checkIfExists: true)
        .readLines()
        .each { it ->
            samplename = (it =~ /^(\S+)/)[0][1] // FIXME doens't work when last line of file is blank
            if(file("${params.indir}/${samplename}/ILLUMINA_DATA/*_R{1,2}.fastq.gz").size() != 2) {
                raise Exception("Couldn't find raw files for ${samplename}. Glob pattern: ${params.indir}/${samplename}/ILLUMINA_DATA/*_R{1,2}.fastq.gz.")
            }
            glob_pattern += samplename + ","
        }
    glob_pattern = StringUtils.chop(glob_pattern)
    glob_pattern += "}/ILLUMINA_DATA/*_R{1,2}.fastq.gz"
} else  {
    glob_pattern = "${params.indir}/*/ILLUMINA_DATA/*_R{1,2}.fastq.gz"
}
Channel.fromFilePairs(glob_pattern, size: -1) { file -> file.getParent().getParent().getName() }
        .set { to_combine_reads }

// Defines channels for resources file
artefacts = file(params.artefacts, type: "file", checkIfExists: true )
phix174ill = file(params.phix174ill, type: "file", checkIfExists: true)
adapters = file(params.adapters, type: "file", checkIfExists: true)

process merge_paired_end_reads {
	tag "$name"

	input:
    tuple val(name), path(reads) from to_combine_reads

	output:
	tuple val(name), path("${name}.fq.gz") into raw_reads_to_qa, raw_reads_to_prepro
	val(name) into all_sample_names

	script:
	"""
	 #Sets the maximum memory to 4/5 of the value requested in the config file
    maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') * 4 / 5))G\"

    # concat read pairs
    reformat.sh -Xmx\"\$maxmem\" in1="${reads[0]}" in2="${reads[1]}" out="${name}.fq.gz"
	"""
}

// ------------------------------------------------------------------------------
//    QUALITY CONTROL
// ------------------------------------------------------------------------------

/**
    Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

    This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
    but not otherwise (identical reads will identify natural duplicates).

    Quality control - STEP 2. A decontamination of synthetic sequences and artefacts
    is performed.

    Quality control - STEP 3. Trimming of low quality bases and of adapter sequences.
    Short reads are discarded.

    If dealing with paired-end reads, when either forward or reverse of a paired-read
    are discarded, the surviving read is saved on a file of singleton reads.

    Quality control - STEP 4. Decontamination. Removes external organisms' contamination,
    using given genomes.
*/

process preprocess {
    tag "$name"

    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "*.fq.gz", overwrite: true
    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "01_dedup_log.txt", overwrite: true
    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "02_syndecontam_log.txt", overwrite: true
    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "03_trim+hostremove_log.txt", overwrite: true

    input:
	tuple val(name), path(reads) from raw_reads_to_prepro
    file(artefacts) from artefacts
    file(phix174ill) from phix174ill
    file(adapters) from adapters
    path(ref_foreign_genome) from ref_foreign_genome

    output:
    tuple val(name), path("01_dedup_log.txt") into dedup_log
    tuple val(name), path("02_syndecontam_log.txt") into synthetic_contaminants_log
    tuple val(name), path("03_trim+hostremove_log.txt") into trim_decontam_log
    tuple val(name), path("${name}.QCd.fq.gz") into qcd_reads_to_qa, to_profile_taxa

    // parameters params.singleEnd, params.qc_matepairs, params.dedup, params.mode come into effect here

    when:
    !params.skip_preprocess

    script:
    // // When paired-end are used, decontamination is carried on independently on paired reads
    // // and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads
    // // and once on the singleton ones, merging the results on a single output file
    // def decontam_in = params.singleEnd ? "in=\"${name}_trimmed.fq.gz\"" :  "in1=\"${name}_trimmed_R1.fq.gz\",\"${name}_trimmed_singletons.fq.gz\" in2=\"${name}_trimmed_R2.fq.gz\",null"
    // def outputu = "\"decontam/${name}.fq.gz\""
    // def outputm = "\"${name}_contamination.fq.gz\""

    """
    #Sets the maximum memory to 4/5 of the value requested in the config file
    maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') * 4 / 5))G\"

    # dedup
    clumpify.sh -Xmx\"\$maxmem\" in=${name}.fq.gz out=\"${name}.dedup.fq.gz\" qin=$params.qin dedupe \\
        reorder subs=0 threads=${task.cpus} &> 01_dedup_log.txt

    # MultiQC doesn't have a module for clumpify yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    #mkdir dedup
    #bash scrape_dedup_log.sh > dedup/${name}.yaml
    # bash scrape_dedup_log.sh 01_dedup_log.txt ${name} > 01_dedup_log.txt # convert to MQC table format

    # remove synthetic contaminants
    bbduk.sh -Xmx\"\$maxmem\" in=\"${name}.dedup.fq.gz\" \\
        out=\"${name}.no_synthetic_contaminants.fq.gz\" k=31 ref=$phix174ill,$artefacts \\
        qin=$params.qin ordered=t threads=${task.cpus} ow &> 02_syndecontam_log.txt
    rm \"${name}.dedup.fq.gz\"
    # rm \"${name}_dedup_R2.fq.gz\"

    # MultiQC doesn't have a module for bbduk yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    #mkdir syndecontam
    #bash scrape_remove_synthetic_contaminants_log.sh > syndecontam/${name}.yaml

    # trim
    kneaddata -i \"${name}.no_synthetic_contaminants.fq.gz\" \\
        -o ./ --output-prefix=\"${name}\" --remove-intermediate-output --log 03_trim+hostremove_log.txt \\
        --trimmomatic-options="ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:4:$params.phred MINLEN:$params.minlength" \\
        -db $ref_foreign_genome \\
        --max-memory \"\$maxmem\" --threads=${task.cpus}
    rm \"${name}.no_synthetic_contaminants.fq.gz\" # input to kneaddata, no longer needed
    gzip \"${name}.fastq\" --stdout > \"${name}.QCd.fq.gz\" && rm \"${name}.fastq\" # gzip final file
    foreign_genome_fn=\$( basename $ref_foreign_genome/*.rev.1.bt2 )
    foreign_name=\${foreign_genome_fn%.rev.1.bt2}
    gzip \"${name}_\${foreign_name}_bowtie2_contam.fastq\" --stdout > \"${name}.contamination.fq.gz\" \\
        && rm \"${name}_\${foreign_name}_bowtie2_contam.fastq\"
    """
}


// if(params.skip_preprocess) {
// 	raw_reads_to_profile_taxa.mix(qcd_reads_to_profile_taxa)
// 		.set{to_profile_taxa}
// }

// ------------------------------------------------------------------------------
//    QUALITY ASSESSMENT
// ------------------------------------------------------------------------------

process quality_assessment {

    tag "$name"

    //Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_fastqc
    } else {
        container params.docker_container_fastqc
    }

    publishDir "${params.outdir}/${name}/ANALYSIS/fastqc", mode: 'link', overwrite: true

    input:
    tuple val(name), file(reads: 'raw/*'), file(qcd_reads: "qcd/${name}.fq.gz") from raw_reads_to_qa.join(qcd_reads_to_qa, failOnDuplicate: true, failOnMismatch: false)

    output:
    tuple(val(name), path("raw/*_fastqc.zip")) into fastqc_raw_log
    tuple(val(name), path("qcd/*_fastqc.zip")) optional true into fastqc_qcd_log

    when:
    !params.skip_quality_assessment

    script:
    """
    fastqc -q $reads
    fastqc -q $qcd_reads || true
    """
}

// ------------------------------------------------------------------------------
//  COMMUNITY CHARACTERISATION
// ------------------------------------------------------------------------------

/**
    Community Characterisation - STEP 1. Performs taxonomic binning and estimates the
    microbial relative abundances using MetaPhlAn and its databases of clade-specific markers.
*/

bowtie2_metaphlan_databases = file( params.metaphlan_databases, type: 'dir', checkIfExists: true )

process profile_taxa {
    label 'biobakery'
    tag "$name"

    //Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "*.{biom,tsv}", overwrite: true
    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "${name}.metaphlan_bugs_list.tsv", overwrite: true

    input:
    tuple val(name), file(reads) from to_profile_taxa
    file(bowtie2db) from bowtie2_metaphlan_databases

    output:
    tuple val(name), path("*.biom") into to_alpha_diversity
    file("*.metaphlan_bugs_list.tsv") into to_collect_taxonomic_profiles, nSamples
    tuple val(name), path(reads), path("*.metaphlan_bugs_list.tsv") into to_profile_function
    tuple val(name), path("${name}.metaphlan_bugs_list.tsv") into profile_taxa_log

    when:
    !params.skip_profile_taxa

    script:
    """
    #If a file with the same name is already present, Metaphlan2 used to crash, leaving this here just
    # in case
    rm -rf ${name}_bt2out.txt

    metaphlan --input_type fastq --tmp_dir=. --biom ${name}.biom --bowtie2out=${name}_bt2out.txt \\
        --bowtie2db $bowtie2db --bt2_ps ${params.bt2options} --add_viruses --sample_id ${name} \\
        --nproc ${task.cpus} $reads ${name}.metaphlan_bugs_list.tsv &> profile_taxa_mqc.txt

    # MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
    # had to create a YAML file with all the info I need via a bash script
    bash scrape_profile_taxa_log.sh ${name}.metaphlan_bugs_list.tsv > ${name}.yaml
    """
}

if (params.skip_profile_taxa) {
    Channel.fromPath("${params.indir}/*/ANALYSIS/*.metaphlan_bugs_list.tsv")
        .set { to_collect_taxonomic_profiles; nSamples }
}

nSamples = nSamples.collect().count()

process collect_taxonomic_profiles {
    label 'biobakery'

    publishDir "${params.outdir}", mode: 'link', pattern: "merged_metaphlan_abundance_table.txt", overwrite: true
    publishDir "${params.outdir}", mode: 'link', pattern: "merged_metaphlan_abundance_table.species.txt", overwrite: true

    input:
    file '*' from to_collect_taxonomic_profiles.collect()

    output:
    file('merged_metaphlan_abundance_table.species.txt') into to_hclust_taxonomic_profiles

    when:
    !params.skip_profile_taxa

    script:
    """
    # merge_metaphlan_tables.py will put the whole filename (without .tsv) as column names
    for i in *.metaphlan_bugs_list.tsv; do
        mv \$i \${i/.metaphlan_bugs_list/}
    done

    merge_metaphlan_tables.py *.tsv > merged_metaphlan_abundance_table.txt
    head -n 2 merged_metaphlan_abundance_table.txt > merged_metaphlan_abundance_table.species.txt
    grep "s__" merged_metaphlan_abundance_table.txt >> merged_metaphlan_abundance_table.species.txt
    """
}

process hclust_taxonomic_profiles {
    label 'hclust'

    publishDir "${params.outdir}", mode: 'link', overwrite: true

    input:
    path(table_species) from to_hclust_taxonomic_profiles
    val nSamples

    output:
    path '*.png'

    when:
    !params.skip_profile_taxa

    script:
    """
    aspect=\$( echo "$nSamples*0.03" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
    aspect=\$( [[ "\$aspect" -gt 9 ]] && echo "9" || echo "\$aspect" ) # bounded to 1 and 9
    aspect=\$( [[ "\$aspect" -lt 1 ]] && echo "1" || echo "\$aspect" )

    hclust2.py -i ${table_species} \\
        -o metaphlan.species.top50.hclust2.png --skip_rows 0 --ftop 50 \\
        --f_dist_f cosine --s_dist_f braycurtis --slinkage complete --no_fclustering \\
        --sqrt_scale \\
        --flabel_size 4 --no_slabels --max_flabel_len 100 \\
        --title "% Relative abundance" --dpi 300 --cell_aspect_ratio \$aspect
    """
}

process alpha_diversity {

    tag "$name"

    //Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_qiime2
    } else {
        container params.docker_container_qiime2
    }

    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "*.{tsv}", overwrite: true

    input:
    tuple val(name), file(metaphlan_bug_list) from to_alpha_diversity

    output:
    file "${name}.alpha-diversity.tsv" into alpha_diversity_log

    when:
    !params.skip_profile_taxa

    script:
    """
    #It checks if the profiling was successful, that is if identifies at least three species
    n=\$(grep -o s__ $metaphlan_bug_list | wc -l  | cut -d\" \" -f 1)
    if (( n <= 3 )); then
        #The file should be created in order to be returned
        touch ${name}.alpha-diversity.tsv
    else
        echo $name > ${name}.alpha-diversity.tsv
        qiime tools import --input-path $metaphlan_bug_list --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ${name}_abundance_table.qza
        for alpha in berger_parker_d brillouin_d dominance enspie esty_ci fisher_alpha gini_index goods_coverage heip_e kempton_taylor_q lladser_pe margalef mcintosh_d mcintosh_e menhinick michaelis_menten_fit observed_features osd pielou_e robbins shannon simpson simpson_e singles strong
        do
            qiime diversity alpha --i-table ${name}_abundance_table.qza --p-metric \$alpha --output-dir \$alpha
            qiime tools export --input-path \$alpha/alpha_diversity.qza --output-path \${alpha}
            if [[ -e "\${alpha}/alpha-diversity.tsv" ]]; then
                value=\$(sed -n '2p' \${alpha}/alpha-diversity.tsv | cut -f 2)
            else
                value="#N/A"
            fi
            echo -e  \$alpha'\t'\$value >> ${name}.alpha-diversity.tsv
        done
    fi
    """
}

/**
    Community Characterisation - STEP 2. Performs the functional annotation using HUMAnN.
*/

// Defines channels for bowtie2_metaphlan_databases file
chocophlan_databases = file( params.chocophlan, type: 'dir', checkIfExists: true )
uniref_databases = file( params.uniref, type: 'dir', checkIfExists: true )

process profile_function {
    label 'biobakery'
    tag "$name"

    //Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link',
        pattern: "*.{tsv,log}", overwrite: true
    publishDir "${params.outdir}/${name}/ANALYSIS/humann_intermediate",
        mode: 'link', pattern: "${name}_humann_temp/*",
        saveAs: {
            filename -> f = new File(filename)
            f.name
        }, overwrite: true

    input:
    tuple val(name), file(reads), file(metaphlan_bug_list) from to_profile_function
    file(chocophlan) from chocophlan_databases
    file(uniref) from uniref_databases

    output:
    tuple val(name), path("*.HUMAnN.log") into profile_functions_log
    path "${name}_humann_temp/*"
    file "*_genefamilies.tsv" into humann_gene_families
    file "*_pathcoverage.tsv" into humann_path_coverage
    file "*_pathabundance.tsv" into humann_path_abundance

    when:
    !params.skip_profile_function

    script:
    """
    #HUMAnN will uses the list of species detected by the profile_taxa process
	
	if [[ -n "$metaphlan_bug_list" ]]; then
		humann --input $reads --output . --output-basename ${name} --taxonomic-profile $metaphlan_bug_list \
			--nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads \
			${task.cpus} --memory-use maximum &> ${name}.HUMAnN.log
	else
		humann --input $reads --output . --output-basename ${name} \
			--nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads \
			${task.cpus} --memory-use maximum &> ${name}.HUMAnN.log
	fi
    
    rm ${name}_humann_temp/*.bt2

    shopt -s nullglob
    for file in ${name}_humann_temp/*.sam; do bgzip \$file --compress-level 9; done
    for file in ${name}_humann_temp/*.ffn; do bgzip \$file --compress-level 9; done
    for file in ${name}_humann_temp/*.tsv; do bgzip \$file --compress-level 9; done
    for file in ${name}_humann_temp/*.fa; do bgzip \$file --compress-level 9; done
     """
}

process collect_functional_profiles {
    label 'biobakery'

    publishDir "${params.outdir}", mode: 'link', overwrite: true,
        pattern: "{all_genefamilies*.tsv,all_pathcoverage*.tsv,all_pathabundance*.tsv}"

    //Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

    input:
    file('*') from humann_gene_families.collect()
    file('*') from humann_path_coverage.collect()
    file('*') from humann_path_abundance.collect()

    output:
    path('*.tsv') into collected_functional_profiles
    path('all_genefamilies.copm.stratified.tsv') into gene_fams_to_hclust
    path('all_pathabundance.copm.unstratified.no-unintegrated.tsv') into pathabunds_to_hclust
    path('all_pathcoverage.unstratified.tsv') into pathcovs_to_hclust

    when:
    !params.skip_profile_function

    script:
 /* groovylint-disable-next-line UnnecessaryGString */
    """
    humann_join_tables -i ./ -o all_genefamilies.rpk.tsv --file_name genefamilies
    humann_join_tables -i ./ -o all_pathcoverage.temp.tsv --file_name pathcoverage
    humann_join_tables -i ./ -o all_pathabundance.rpk.tsv --file_name pathabundance

    # UNMAPPED and UNINTEGRATED coverage values are meaningless
    grep -vP "(?:^UNMAPPED)|(?:^UNINTEGRATED)" all_pathcoverage.temp.tsv > all_pathcoverage.tsv
    rm -f all_pathcoverage.temp.tsv

    humann_renorm_table -i all_genefamilies.rpk.tsv -o all_genefamilies.copm.tsv --units cpm

    humann_split_stratified_table --input all_genefamilies.copm.tsv --output ./
    rm all_genefamilies.copm.tsv
    mv all_genefamilies.copm_stratified.tsv all_genefamilies.copm.stratified.tsv
    mv all_genefamilies.copm_unstratified.tsv all_genefamilies.copm.unstratified.tsv

    humann_renorm_table -i all_pathabundance.rpk.tsv -o all_pathabundance.copm.tsv --units cpm

    humann_split_stratified_table --input all_pathabundance.copm.tsv --output ./
    rm all_pathabundance.copm.tsv
    mv all_pathabundance.copm_stratified.tsv all_pathabundance.copm.stratified.tsv
    mv all_pathabundance.copm_unstratified.tsv all_pathabundance.copm.unstratified.tsv

    humann_split_stratified_table --input all_pathcoverage.tsv --output ./
    mv all_pathcoverage_stratified.tsv all_pathcoverage.stratified.tsv
    mv all_pathcoverage_unstratified.tsv all_pathcoverage.unstratified.tsv

    grep -vP "(?:^UNMAPPED)|(?:^UNINTEGRATED)" all_pathabundance.copm.unstratified.tsv > all_pathabundance.copm.unstratified.no-unintegrated.tsv
    """
}

process hclust_functional_profiles {
    label 'hclust'

    publishDir "${params.outdir}", mode: 'link', overwrite: true

    input:
    path gene_fam_table from gene_fams_to_hclust
    path path_abund_table from pathabunds_to_hclust
    path path_cov_table from pathcovs_to_hclust
    val nSamples

    output:
    path '*.png'

	when:
    !params.skip_profile_function

    script:
    """
    if [[ \$(awk '{print NF}' ${gene_fam_table} | head -n 2 | tail -n 1) -gt 3 ]]; then
        cluster_columns=""
    else
        cluster_columns="--no_fclustering"
    fi

    aspect=\$( echo "$nSamples*0.03" | bc -l | awk '{printf("%d\\n",\$1 + 0.5)}' )
    aspect=\$( [[ "\$aspect" -gt 9 ]] && echo "9" || echo "\$aspect" ) # bounded to 1 and 9
    aspect=\$( [[ "\$aspect" -lt 1 ]] && echo "1" || echo "\$aspect" )

    file="$gene_fam_table"
    hclust2.py -i ${gene_fam_table} -o \${file%.*}.top50.hclust2.png \\
        --skip_rows 0 --ftop 50 \\
        \$cluster_columns --f_dist_f cosine --s_dist_f braycurtis --slinkage complete \\
        --dpi 300 --sqrt_scale \\
        --no_slabels --flabel_size 4 --colorbar_font_size 4 --max_flabel_len 100 \\
        --title "Gene family abundance (copies per million)" --cell_aspect_ratio \$aspect

    file="$path_cov_table"
    hclust2.py -i ${path_cov_table} -o \${file%.*}.top50.hclust2.png \\
         --skip_rows 0 --ftop 50 \\
        \$cluster_columns --f_dist_f cosine --s_dist_f cityblock --slinkage complete \\
        --dpi 300 --sqrt_scale \\
        --no_slabels --flabel_size 4 --colorbar_font_size 4 --max_flabel_len 100 \\
        --title "Pathway coverage" --cell_aspect_ratio \$aspect

    file="$path_abund_table"
    hclust2.py -i ${path_abund_table} -o \${file%.*}.top50.hclust2.png --skip_rows 0 --ftop 50 \\
        \$cluster_columns --f_dist_f cosine --s_dist_f canberra --slinkage complete \\
        --dpi 300 --sqrt_scale \\
        --no_slabels --flabel_size 4 --colorbar_font_size 4 --max_flabel_len 100 \\
        --title "Pathway abundance (copies per million)" --cell_aspect_ratio \$aspect
    """
}

// if(params.skip_preprocess) {
// 	all_sample_names.map{ it -> [it, "", "", ""] }
// 		.set{ log_tuple_1 }
// } else {
// 	dedup_log.join(synthetic_contaminants_log, remainder: true )
// 		.join( trim_decontam_log, remainder: true )
// 		.set(log_tuple_1)
// }
dedup_log.join(synthetic_contaminants_log, remainder: true )
	.join( trim_decontam_log, remainder: true )
	.join( profile_taxa_log, remainder: true )
	.join( profile_functions_log, remainder: true )
    .set{ all_logs }

// This should happen RIGHT after finishing HUMANN. Collect the individual tool logs into one
// tool-agnostic data file that is forward-compatible to any future pipeline runs.
process sample_multiqc_logging {

    tag "$name"

    publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "04_*_stats.txt", overwrite: true

    input:
    tuple val(name), path(dedup_log), path(syndecontam_log), path(trim_log), path(metaphlan_log), path(humann_log) from all_logs

    output:
    path("04_qc_stats.txt") into sample_genstats
	path("04_profiling_stats.txt") into sample_profilestats

    script:
    """
    printf "Sample\\traw\\tdeduped\\tsynDecontam\\ttrimmed\\thost\\tqcd\\n" > 04_qc_stats.txt

	if [[ -n "$dedup_log" ]]; then
		totR=\$(grep "Reads In:" "$dedup_log" | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
		remR=\$(grep "Duplicates Found:" "$dedup_log" | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
		deduped=\$((\$totR - \$remR))
	else
		totR=""
		deduped=""
	fi
    printf "%s\\t%s\\t%s" "\\"$name\\"" "\${totR}" "\${deduped}" >> 04_qc_stats.txt

	if [[ -n "$syndecontam_log" ]]; then
		contam=\$(grep "Contaminants:" "$syndecontam_log" | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
		syndecontam=\$((\$deduped-\$contam))
	else
		syndecontam=""
	fi
    sedstr="s/(\\"$name\\"\\t.*)\$/\\1\\t\${syndecontam}/"
    sed -i -E "\${sedstr}" 04_qc_stats.txt

	if [[ -n "$trim_log" ]]; then
		filtered=\$(grep -oP "TrimmomaticSE: Started with arguments:.*Dropped: [0-9]*" "$trim_log" | grep -oP "[0-9]*\$" )
		trimmed=\$((\$syndecontam-\$filtered))

		final=\$(grep -P "READ COUNT: final single" "$trim_log" | grep -oP "[0-9\\.]*\$" | grep -oP "^[0-9]*" )
    	host=\$((\$trimmed-\$final))
	else
		trimmed=""
		host=""
		final=""
	fi
    sedstr="s/(\\"$name\\"\\t.*)\$/\\1\\t\${trimmed}\\t\${host}\\t\${final}\\n/"
    sed -i -E "\${sedstr}" 04_qc_stats.txt

	printf "Sample\\tnSpecies\\tgene_families_nucleotide\\tcontributing_nucleotide" > 04_profiling_stats.txt
    printf "\\tgene_families_final\\tn_contributing_final\\t" >> 04_profiling_stats.txt
    printf "perc_contributing_final\\tunmapped\\n" >> 04_profiling_stats.txt
	if [[ -n "$metaphlan_log" ]]; then
    	tot_species=\$(grep "s__" "$metaphlan_log" | wc -l)
	else
		tot_species=""
	fi
    printf "%s\\t%s\\n" "\\"$name\\"" \${tot_species} >> 04_profiling_stats.txt

	if [[ -n "$humann_log" ]]; then
    	gene_fams_nucleotide=\$(grep "Total gene families from nucleotide alignment:" "$humann_log" | cut -d: -f 2 | sed 's/ //g')
    	unaligned_nucleotide=\$(grep -o "Unaligned reads after nucleotide alignment: [0-9\\.]*" "$humann_log" | grep -o "[0-9\\.]*")
    	contributing_nucleotide=\$(echo "100 - \$unaligned_nucleotide" | bc)

		unaligned_final=\$(grep "Unaligned reads after translated alignment:" "$humann_log" | grep -o "[0-9\\.]*")
		contributing_final=\$(echo "100 - \$unaligned_final" | bc)
		nContributing_final=\$(echo "\$contributing_final / 100 * \$final" | bc -l | awk '{printf("%d\\n",\$1 + 0.5)}' )
		unmapped=\$((\$final-\$nContributing_final))
		gene_fams_final=\$(grep "Total gene families after translated alignment:" "$humann_log" | cut -d: -f 2 | sed 's/ //g')
	else
		gene_fams_nucleotide=""
		contributing_nucleotide=""
		gene_fams_final=""
		nContributing_final=""
		contributing_final=""
		unmapped=""
	fi
    sedstr="s/(\\"$name\\"\\t.*)\$/\\1\\t\${gene_fams_nucleotide}\\t\${contributing_nucleotide}\\t\${gene_fams_final}\\t\${nContributing_final}\\t\${contributing_final}\\t\${unmapped}/"
    sed -i -E "\${sedstr}" 04_profiling_stats.txt
    """
}

// ------------------------------------------------------------------------------
//    MULTIQC LOGGING
// ------------------------------------------------------------------------------


/**
    Generate Logs.

    Logs generate at each analysis step are collected and processed with MultiQC
*/

// These skip-parameters ...how do they interact with --only-samples?
// When ready...fix this. please. skip means to change the DSL flow so that it's as if the process wasn't done at all.
// if(params.skip_preprocess) {
//     Channel.fromPath("${params.indir}/*/ANALYSIS/01_dedup_log.txt")
//         .set { dedup_log }
//     Channel.fromPath("${params.indir}/*/ANALYSIS/02_syndecontam_log.txt")
//         .set { synthetic_contaminants_log }
//     Channel.fromPath("${params.indir}/*/ANALYSIS/03_trim+decontam_log.txt")
//         .set { trimming_log }
// }
// if(params.skip_quality_assessment) {
//     // find multiQC inputs from files
//     Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/raw/*")
//         .set { fastqc_raw_log }
//     Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/qcd/*")
//         .set { fastqc_qcd_log }
// }
// if(params.skip_profile_taxa) {
//     Channel.fromPath("${params.indir}/*/ANALYSIS/*.metaphlan_bugs_list.tsv")
//         .set { profile_taxa_log }
//     Channel.fromPath("${params.indir}/*/ANALYSIS/*.alpha-diversity.tsv")
//         .set {alpha_diversity_log}
//     Channel.fromPath("${params.indir}/*/ANALYSIS/*.HUMAnN.log")
//         .set { profile_functions_log }
// }


// Get sample summary stats from any other samples in the project directory
// Any way to ensure we only have 1 file per sample (in case --combine_multiqc and --only_samples both used and there was an overlap?)
archived_genstats = Channel.empty()
archived_profilestats = Channel.empty()
if(params.combine_multiqc) {
	// get only_samples
	// get globbed files
	// prioritize data files from THIS run of the pipeline, if an archived data file also exists for the same sample
	Channel.fromPath("${params.indir}/qc_stats.txt")
		.set { archived_genstats }
	Channel.fromPath("${params.indir}/profiling_stats.txt")
		.set { archived_profilestats }
	Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/raw/*_fastqc.zip")
        .map { it -> [file(it).getParent().getParent().getParent().getParent().getName(), it] }
        .join(fastqc_raw_log, remainder: true)
		.map { it -> [it[2] != null ? it[2] : it[1]] }
		.set { fastqc_raw_log }
	Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/qcd/*_fastqc.zip")
        .map { it -> [file(it).getParent().getParent().getParent().getParent().getName(), it] }
        .join(fastqc_qcd_log, remainder: true)
		.map { it -> [it[2] != null ? it[2] : it[1]] }
		.set { fastqc_qcd_log }
} else {
	// Drop the keys, they will not be needed or accepted by the next process
	fastqc_raw_log.map {it -> it[1]}
		.set { fastqc_raw_log }
	fastqc_qcd_log.map {it -> it[1]}
		.set { fastqc_qcd_log }
}

process log {

    publishDir "${params.outdir}", mode: 'link', overwrite: true

    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_multiqc
    } else {
        container params.docker_container_multiqc
    }

    stageInMode "copy"

    input:
    file multiqc_config from file(params.multiqc_config, type: 'file', checkIfExists: true )
    path "yamp.css" from file(params.multiqc_css, type: 'file', checkIfExists: true )
    path "replace_sample_names.txt" from file(params.multiqc_replace_names, type: 'file', checkIfExists: true )
    file workflow_summary from create_workflow_summary(summary)
    file "software_versions_mqc.yaml" from software_versions_yaml.collect()
    path "fastqc_raw/*" from fastqc_raw_log.collect().ifEmpty([])
    path "fastqc_QCd/*" from fastqc_qcd_log.collect().ifEmpty([])
    path("sample_stats/genstats_*.txt") from sample_genstats.collect().ifEmpty([])
	path("sample_stats/profiling_stats_*.txt") from sample_profilestats.collect().ifEmpty([])
	path("old_qc_stats.txt") from archived_genstats.ifEmpty([])
	path("old_profiling_stats.txt") from archived_profilestats.ifEmpty([])

    output:
	path("qc_stats.txt")
	path("profiling_stats.txt")
    path "multiqc_report.html"
    path "multiqc_data/*"

    script:
    """
    # collect stats files into one project-level table
    statsfiles=( sample_stats/genstats*.txt )
    nFiles=\${#statsfiles[@]}
	for (( i=0; i<\$nFiles; i++ ))
	do
		if [[ \$i -eq 0 ]]; then
			cat \${statsfiles[\$i]} > qc_stats.txt
		else
			head -2 \${statsfiles[\$i]} | tail -1 >> qc_stats.txt
		fi
	done

    statsfiles=( sample_stats/profiling_stats*.txt )
    nFiles=\${#statsfiles[@]}
	for (( i=0; i<\$nFiles; i++ ))
	do
		if [[ \$i -eq 0 ]]; then
			cat \${statsfiles[\$i]} > profiling_stats.txt
		else
			head -2 \${statsfiles[\$i]} | tail -1 >> profiling_stats.txt
		fi
	done

	# combine with old stats, if called for
    if [[ -e "old_qc_stats.txt" ]]; then
		python -s \$( which combine_tsvs.py ) "old_qc_stats.txt" "qc_stats.txt"
    fi
    if [[ -e "old_profiling_stats.txt" ]]; then
		python -s \$( which combine_tsvs.py ) "old_profiling_stats.txt" "profiling_stats.txt"
    fi

    # extract columns verbatim for a barplot
    printf "Sample\\tHost\\tMicrobial\\tUnmapped\\n" > host_microbial_composition.txt
    {
        read # skip first line
        while read sample raw deduped synDecontam trimmed host qcd
        do
            printf "%s\\t%s\\n" "\$sample" "\$host" >> host_microbial_composition.txt
        done
    } < qc_stats.txt
    {
        read
        while read sample nSpecies gene_families_nucleotide contributing_nucleotide gene_families_final n_contributing_final perc_contributing_final unmapped
        do
            sedstr="s/(\${sample}\\t.*)\$/\\1\\t\${n_contributing_final}\\t\${unmapped}/"
            sed -i -E "\${sedstr}" host_microbial_composition.txt
        done
    } < profiling_stats.txt

    multiqc --config $multiqc_config . -f --custom-css-file yamp.css
    """
}
