#!/usr/bin/env nextflow

/**
MD Genomics Metagenomic Analysis Pipeline
Copyright (C) 2022 Jonathan Lim.

This program incorporates a modified version of Yet Another Metagenomic Pipeline (YAMP).
Copyright (C) 2017-2021	Dr Alessia Visconti

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
Copyright (C) 2017-2021	Dr Alessia Visconti

This pipeline is distributed in the hope that it will be useful
but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.

  Usage:
  nextflow run YAMP.nf --reads1 R1 --reads2 R2 --prefix prefix --outdir path [options]

  Mandatory arguments:
    --reads1   R1      Forward (if paired-end) OR all reads (if single-end) file path
    [--reads2] R2      Reverse reads file path (only if paired-end library layout)
    --prefix   prefix  Prefix used to name the result files
    --outdir   path    Output directory (will be outdir/prefix/)

  Main options:
    --mode       <QC|characterisation|complete>
    --singleEnd  <true|false>   whether the layout is single-end
	--qc_matepairs  <true|false> whether to wait until after QC to combine mate pairs
    --dedup      <true|false>   whether to perform de-duplication

  Profiles:
	--profile msl	Profile for all MSL analysis.
						Dedup: false
						Decontaminate: true
						Adapters: /local/projects/drasko/thazen/scripts/illumina.adapters.current.fasta.txt
						Phix174ill:  /local/projects/drasko/thazen/scripts/phiX174.fas

  Other options:
  BBduk parameters for removing synthetic contaminants and trimming:
    --qin                 <33|64> Input quality offset
    --kcontaminants       value   kmer length used for identifying contaminants
    --phred               value   regions with average quality BELOW this will be trimmed
    --minlength           value   reads shorter than this after trimming will be discarded
    --mink                value   shorter kmer at read tips to look for
    --hdist               value   maximum Hamming distance for ref kmer
    --artefacts           path    FASTA file with artefacts
    --phix174ill          path    FASTA file with phix174_ill
    --adapters            path    FASTA file with adapters

  BBwrap parameters for decontamination:
    --foreign_genome      path    FASTA file for contaminant (pan)genome
    --foreign_genome_ref  path    folder for for contaminant (pan)genome (pre indexed)
    --mind                value   approximate minimum alignment identity to look for
    --maxindel            value   longest indel to look for
    --bwr                 value   restrict alignment band to this

  MetaPhlAn parameters for taxa profiling:
    --metaphlan_databases path    folder for the MetaPhlAn database
    --bt2options          value   BowTie2 options

  HUMANn parameters for functional profiling:
    --chocophlan          path    folder for the ChocoPhlAn database
    --uniref              path	  folder for the UniRef database

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
// 	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
// }

// Given --indir, this script knows how to find the inputs for any entrypoint
if (params.indir == "null") {
	exit 1, "Please set the --indir parameter"
}

//Creates working dir
// workingpath = params.outdir + "/" + params.prefix
// workingdir = file(workingpath)
// if( !workingdir.exists() ) {
// 	if( !workingdir.mkdirs() ) 	{
// 		exit 1, "Cannot create working directory: $workingpath"
// 	}
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
	   	summary['BBmap'] = "https://depot.galaxyproject.org/singularity/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
	} else if (workflow.containerEngine == 'docker') {
    	summary['BBmap'] = "quay.io/biocontainers/bbmap:38.87--h1296035_0"
		summary['FastQC'] = "quay.io/biocontainers/fastqc:0.11.9--0"
	} else {
		summary['BBmap'] = "No container information"
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
summary['Prefix'] = params.prefix
summary['Layout'] = 'Paired-End'
// summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'

if (!params.skip_preprocess)
{
	// if (!params.singleEnd)
	// {
	// 	summary['QC mate pairs separately'] = params.qc_matepairs
	// }
	summary['Performing de-duplication'] = params.dedup

	//remove_synthetic_contaminants
	summary['Synthetic contaminants:'] = ""
	summary['Artefacts'] = params.artefacts
	summary['Phix174ill'] = params.phix174ill

	//Trimming
	summary['Adapters'] = params.adapters
	summary['Trimming parameters:'] = ""
	summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'
	summary['Min phred score'] = params.phred
	summary['Min length'] = params.minlength
	summary['kmer lenght'] = params.kcontaminants
	summary['Shorter kmer'] = params.mink
	summary['Max Hamming distance'] = params.hdist

	//Decontamination
	summary['Decontamination parameters:'] = ""
	if (params.foreign_genome_ref != "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome_ref + " (indexed)"
	} else if (	params.foreign_genome_ref == "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome
	}
	summary['Min alignment identity'] = params.mind
	summary['Max indel length'] = params.maxindel
	summary['Max alignment band'] = params.bwr
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
// 	folder = new File(params.outdir + "/" + name + "/ANALYSIS/")
// 	if( !folder.exists() ) {
// 		folder.mkdirs()
// 	}
// 	return folder.getPath()
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
foreign_genome = file( "${params.foreign_genome}", type: "file", checkIfExists: true )

//Stage boilerplate log when the contaminant (pan)genome is indexed
if (!params.skip_preprocess && params.foreign_genome_ref == "") {
	index_foreign_genome_log = Channel.from(file("$baseDir/assets/foreign_genome_indexing_mqc.yaml"))
} else {
	index_foreign_genome_log = Channel.empty()
}

process index_foreign_genome {

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_bbmap
    } else {
        container params.docker_container_bbmap
    }

	input:
	file(foreign_genome) from foreign_genome

	output:
	path("ref/", type: 'dir') into ref_foreign_genome

	when:
	!params.skip_preprocess && params.foreign_genome_ref == ""

	script:
	"""
	#Sets the maximum memory to the value requested in the config file
	#maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') ))G\"

	# This step will have a boilerplate log because the information saved by bbmap are not relevant
	bbmap.sh -Xmx\"\$maxmem\" ref=$foreign_genome &> foreign_genome_index_mqc.txt
	"""
}

//When the indexed contaminant (pan)genome is already available, its path should be pushed in the correct channel
if (params.foreign_genome_ref != "") {
	ref_foreign_genome = file(params.foreign_genome_ref, checkIfExists: true )
}

/**
	Creates a set of channels for input read files.
	- read_files_fastqc is used for the first QC assessment (on the raw reads)
	- read_files_dedup  is used for the deduplication step (which is optional and may skip to trimming)
	- read_files_trim   is used for the decontamination from synthetic contaminants (used only if
	  deduplication is not run)
*/

if (params.singleEnd) {
	// FIXME UNSUPPORTED
	// Channel
	// .from([[params.prefix, [file(params.reads1)]]])
	// .into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
} else {
	if(params.qc_matepairs) {
		Channel.fromFilePairs("${params.indir}/*/ILLUMINA_DATA/*_R{1,2}.fastq.gz", size: -1) { file -> file.getParent().getParent().getName() }
			.into {read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants; read_files_log }
		// Channel
		// .fromFilePairs("${params.indir}/*_R{1,2}.fastq.gz", checkIfExists: true)
		// //.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
		// .into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants; read_files_log }
	} else {
		Channel.fromFilePairs("${params.indir}/*/ILLUMINA_DATA/*_R{1,2}.fastq.gz", size: -1) { file -> file.getParent().getParent().getName() }
			.set { to_combine_reads }
		// Channel.fromFilePairs("${params.indir}/*_R{1,2}.fastq.gz", checkIfExists: true)
		// .set {to_combine_reads}
	}
}

// Defines channels for resources file
artefacts = file(params.artefacts, type: "file", checkIfExists: true )
phix174ill = file(params.phix174ill, type: "file", checkIfExists: true)
adapters = file(params.adapters, type: "file", checkIfExists: true)

process preprocess {
	tag "$name"

	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "decontam/*.fq.gz", saveAs: {filename -> file(filename).getName()}, overwrite: true
	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "01_dedup_log.txt", overwrite: true
	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "02_syndecontam_log.txt", overwrite: true
	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "03_trim_log.txt", overwrite: true
	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "04_decontam_log.txt", overwrite: true

	input:
	tuple val(name), file(reads) from to_combine_reads
	file(artefacts) from artefacts
	file(phix174ill) from phix174ill
	file(adapters) from adapters
	path(ref_foreign_genome) from ref_foreign_genome

	output:
	tuple val(name), path("${name}.fq.gz") into read_files_fastqc, read_files_log
	file "01_dedup_log.txt" into dedup_log
	file "02_syndecontam_log.txt" into synthetic_contaminants_log
	file "03_trim_log.txt" into trimming_log
	tuple val(name), path("decontam/*.QCd.fq.gz") into qcd_reads
	tuple val(name), path("decontam/*.QCd.fq.gz") into to_profile_taxa_decontaminated
	file "04_decontam_log.txt" into decontaminate_log

	// parameters params.singleEnd, params.qc_matepairs, params.dedup, params.mode come into effect here

	script:
	// myList = [1776, -1, 33, 99, 0, 928734928763]
	// def input = (params.singleEnd || ! params.qc_matepairs) ? ["in=": "\"${reads[0]}\""] :  ["in1=": \"${reads[0]}\"", "in2=":"\"${reads[1]}\""]
	// def deduped = (params.singleEnd || ! params.qc_matepairs) ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""
	// def syndecontam_in = (params.singleEnd || ! params.qc_matepairs) ? "in=\"${name}_dedup.fq.gz\"" :  "in1=\"${name}_dedup_R1.fq.gz\" in2=\"${name}_dedup_R2.fq.gz\""
	// def syndecontam_out = (params.singleEnd || ! params.qc_matepairs) ? "out=\"${name}_no_synthetic_contaminants.fq.gz\"" :  "out=\"${name}_no_synthetic_contaminants_R1.fq.gz\" out2=\"${name}_no_synthetic_contaminants_R2.fq.gz\""
	// def trim_in = (params.singleEnd || ! params.qc_matepairs) ? "in=\"${name}_no_synthetic_contaminants.fq.gz\"" :  "in1=\"${name}_no_synthetic_contaminants_R1.fq.gz\" in2=\"${name}_no_synthetic_contaminants_R2.fq.gz\""
	// def trim_out = params.singleEnd ? "out=\"${name}_trimmed.fq.gz\"" :  "out=\"${name}_trimmed_R1.fq.gz\" out2=\"${name}_trimmed_R2.fq.gz\" outs=\"${name}_trimmed_singletons.fq.gz\""
	// // When paired-end are used, decontamination is carried on independently on paired reads
	// // and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads
	// // and once on the singleton ones, merging the results on a single output file
	// def decontam_in = params.singleEnd ? "in=\"${name}_trimmed.fq.gz\"" :  "in1=\"${name}_trimmed_R1.fq.gz\",\"${name}_trimmed_singletons.fq.gz\" in2=\"${name}_trimmed_R2.fq.gz\",null"
	// def outputu = "\"decontam/${name}.fq.gz\""
	// def outputm = "\"${name}_contamination.fq.gz\""

	"""
	#Sets the maximum memory to 4/5 of the value requested in the config file
	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') * 4 / 5))G\"

	# concat read pairs
	reformat.sh -Xmx\"\$maxmem\" in1="${reads[0]}" in2="${reads[1]}" out="${name}.fq.gz"

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
	bbduk.sh -Xmx\"\$maxmem\" in=\"${name}.no_synthetic_contaminants.fq.gz\" \\
		out=\"${name}.trimmed.fq.gz\" ktrim=r k=$params.kcontaminants mink=$params.mink \\
		hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters \\
		qin=$params.qin ordered=t threads=${task.cpus} tbo tpe ow &> 03_trim_log.txt
	rm \"${name}.no_synthetic_contaminants.fq.gz\"
	# rm \"${name}_no_synthetic_contaminants_R2.fq.gz\"

	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	#mkdir trim
	#bash scrape_trimming_log.sh > trim/${name}.yaml

	# decontaminate
	mkdir decontam

	bbwrap.sh -Xmx\"\$maxmem\" mapper=bbmap append=t in=\"${name}.trimmed.fq.gz\" \\
		outu=\"decontam/${name}.QCd.fq.gz\" outm=\"decontam/${name}.contamination.fq.gz\" \\
		minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl \\
		trimq=$params.phred path="./" qin=$params.qin threads=${task.cpus} untrim quickmatch fast \\
		ordered=t ow &> 04_decontam_log.txt
	rm \"${name}.trimmed.fq.gz\"
	#rm \"${name}_trimmed_singletons.fq.gz\"
	# rm \"${name}_trimmed_R2.fq.gz\"

	# MultiQC doesn't have a module for bbwrap yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script

	"""
}

// ------------------------------------------------------------------------------
//	QUALITY CONTROL
// ------------------------------------------------------------------------------

/**
	Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

	This step is OPTIONAL. De-duplication should be carried on iff you are
    using PCR amplification (in this case identical reads are technical artefacts)
	but not otherwise (identical reads will identify natural duplicates).
*/

// process dedup {

//     tag "$name"

// 	//Enable multicontainer settings
//     if (workflow.containerEngine == 'singularity') {
//         container params.singularity_container_bbmap
//     } else {
//         container params.docker_container_bbmap
//     }

// 	input:
// 	tuple val(name), file(reads) from read_files_dedup

// 	output:
// 	tuple val(name), path("${name}_dedup*.fq.gz") into to_synthetic_contaminants
// 	file "*.yaml" into dedup_log

// 	// when:
// 	// !params.skip_preprocess && params.dedup

// 	script:
// 	// This is to deal with single and paired end reads
// 	def input = (params.singleEnd || ! params.qc_matepairs) ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
// 	def output = (params.singleEnd || ! params.qc_matepairs) ? "out=\"${name}_dedup.fq.gz\"" :  "out1=\"${name}_dedup_R1.fq.gz\" out2=\"${name}_dedup_R2.fq.gz\""

// 	"""
// 	#Sets the maximum memory to 4/5 of the value requested in the config file
// 	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') / 5 * 4))G\"
// 	#maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
// 	echo \"$reads\"
//     clumpify.sh -Xmx\"\$maxmem\" $input $output qin=$params.qin dedupe subs=0 threads=${task.cpus} &> dedup_mqc.txt

// 	# MultiQC doesn't have a module for clumpify yet. As a consequence, I
// 	# had to create a YAML file with all the info I need via a bash script
// 	bash scrape_dedup_log.sh > ${name}.yaml
// 	"""
// }

// /**
// 	Quality control - STEP 2. A decontamination of synthetic sequences and artefacts
// 	is performed.
// */

// //When the de-duplication is not done, the raw file should be pushed in the correct channel
// //FIXME: make this also optional?
// if (!params.dedup & !params.skip_preprocess) {
// 	to_synthetic_contaminants = read_files_synthetic_contaminants
// 	dedup_log = Channel.from(file("$baseDir/assets/no_dedup.yaml"))
// }


// process remove_synthetic_contaminants {

// 	tag "$name"

// 	//Enable multicontainer settings
//     if (workflow.containerEngine == 'singularity') {
//         container params.singularity_container_bbmap
//     } else {
//         container params.docker_container_bbmap
//     }

// 	input:
// 	tuple val(name), file(reads) from to_synthetic_contaminants
// 	file(artefacts) from artefacts
// 	file(phix174ill) from phix174ill

// 	output:
// 	tuple val(name), path("${name}_no_synthetic_contaminants*.fq.gz") into to_trim
// 	file "*.yaml" into synthetic_contaminants_log

// 	when:
// 	!params.skip_preprocess

//    	script:
// 	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
// 	def output = params.singleEnd ? "out=\"${name}_no_synthetic_contaminants.fq.gz\"" :  "out=\"${name}_no_synthetic_contaminants_R1.fq.gz\" out2=\"${name}_no_synthetic_contaminants_R2.fq.gz\""
// 	"""
// 	#Sets the maximum memory to the value requested in the config file
// 	#maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
// 	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') / 5 * 4))G\"
// 	bbduk.sh -Xmx\"\$maxmem\" $input $output k=31 ref=$phix174ill,$artefacts qin=$params.qin threads=${task.cpus} ow &> synthetic_contaminants_mqc.txt

// 	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
// 	# had to create a YAML file with all the info I need via a bash script
// 	bash scrape_remove_synthetic_contaminants_log.sh > ${name}.yaml
// 	"""
// }


/**
	Quality control - STEP 3. Trimming of low quality bases and of adapter sequences.
	Short reads are discarded.

	If dealing with paired-end reads, when either forward or reverse of a paired-read
	are discarded, the surviving read is saved on a file of singleton reads.
*/


// process trim {

// 	tag "$name"

// 	//Enable multicontainer settings
//     if (workflow.containerEngine == 'singularity') {
//         container params.singularity_container_bbmap
//     } else {
//         container params.docker_container_bbmap
//     }

// 	input:
// 	tuple val(name), file(reads) from to_trim
// 	file(adapters) from adapters
// 	output:
// 	tuple val(name), path("${name}_trimmed*.fq.gz") into to_decontaminate
// 	file "*.yaml" into trimming_log

// 	when:
// 	!params.skip_preprocess

//    	script:
// 	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
// 	def output = params.singleEnd ? "out=\"${name}_trimmed.fq.gz\"" :  "out=\"${name}_trimmed_R1.fq.gz\" out2=\"${name}_trimmed_R2.fq.gz\" outs=\"${name}_trimmed_singletons.fq.gz\""
// 	"""
// 	#Sets the maximum memory to the value requested in the config file
// 	#maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
// 	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') / 5 * 4))G\"

// 	bbduk.sh -Xmx\"\$maxmem\" $input $output ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe ow &> trimming_mqc.txt

// 	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
// 	# had to create a YAML file with all the info I need via a bash script
// 	bash scrape_trimming_log.sh > ${name}.yaml
// 	"""
// }


/**
	Quality control - STEP 4. Decontamination. Removes external organisms' contamination,
	using given genomes.

	When an indexed contaminant (pan)genome is not provided, the index_foreign_genome process is run
	before the decontamination process. This process require the FASTA file of the contaminant (pan)genome.
*/


// process decontaminate {

//     tag "$name"

// 	//Enable multicontainer settings
//     if (workflow.containerEngine == 'singularity') {
//         container params.singularity_container_bbmap
//     } else {
//         container params.docker_container_bbmap
//     }

// 	publishDir "${params.outdir}/${name}", mode: 'link', pattern: "*QCd.fq.gz"

// 	input:
// 	tuple val(name), file(reads) from to_decontaminate
// 	path(ref_foreign_genome) from ref_foreign_genome

// 	output:
// 	tuple val(name), path("decontam/*.fq.gz") into qcd_reads
// 	tuple val(name), path("decontam/*.fq.gz") into to_profile_taxa_decontaminated
// 	file "*.yaml" into decontaminate_log

// 	when:
// 	!params.skip_preprocess

// 	script:
// 	// When paired-end are used, decontamination is carried on independently on paired reads
// 	// and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads
// 	// and once on the singleton ones, merging the results on a single output file
// 	def input = params.singleEnd ? "in=\"${reads[0]}\"" :  "in1=\"${reads[0]}\",\"${reads[2]}\" in2=\"${reads[1]}\",null"
// 	def outputu = "\"decontam/${name}.fq.gz\""
// 	def outputm = "\"${name}_contamination.fq.gz\""
// 	"""
// 	mkdir decontam

// 	#Sets the maximum memory to the value requested in the config file
// 	#maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
// 	maxmem=\"\$((\$(echo ${task.memory} | sed 's/ GB//g') / 5 * 4))G\"

// 	bbwrap.sh -Xmx\"\$maxmem\"  mapper=bbmap append=t $input outu=$outputu outm=$outputm minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path="./" qin=$params.qin threads=${task.cpus} untrim quickmatch fast ow &> decontamination_mqc.txt

// 	# MultiQC doesn't have a module for bbwrap yet. As a consequence, I
// 	# had to create a YAML file with all the info I need via a bash script
// 	bash scrape_decontamination_log.sh $outputu $outputm > ${name}.yaml
// 	"""
// }


// ------------------------------------------------------------------------------
//	QUALITY ASSESSMENT
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
    tuple val(name), file(reads: 'raw/*'), file(qcd_reads: "qcd/${name}.fq.gz") from read_files_fastqc.join(qcd_reads, failOnDuplicate: true, failOnMismatch: true)

    output:
    path "raw/*_fastqc.{zip,html}" into fastqc_raw_log
	path "qcd/*_fastqc.{zip,html}" into fastqc_qcd_log

	when:
	!params.skip_preprocess

    script:
    """
    fastqc -q $reads
	fastqc -q $qcd_reads
    """
		// for f in {raw,qcd}/*; do
	// 	base=\${f##*/}   #=> "foo.cpp" (basepath)
	// 	dir=\${f%\${base}}
	// 	ext=\${f##*.}
	// 	mv "\${f}" "\${dir}/${name}_fastqc.\${ext}"
	// done
}

// ------------------------------------------------------------------------------
//  COMMUNITY CHARACTERISATION
// ------------------------------------------------------------------------------

// The user will specify the clean file either as a single clean file (that is the YAMP
// default behaviour), or as two files (forward/reverse). ]
// In the former case, the user will set singleEnd = true and only one file will be
// selected and used directly for taxa and community profiling.
// In the latter case, the user will set singleEnd = false and provide two files, that will
// be merged before feeding the relevant channels for profiling.
if (params.skip_preprocess && params.singleEnd) {
	// FIXME NOT SUPPORTED
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { reads_profile_taxa }

	//Initialise empty channels
	reads_merge_paired_end_cleaned = Channel.empty()
	//Init as value channel so it can be reused by log process
	merge_paired_end_cleaned_log = Channel.value([])
} else if (params.skip_preprocess && (!params.skip_profile_taxa || !params.skip_profile_function) && !params.singleEnd) {
	// FIXME NOT SUPPORTED
	Channel
	.from([[params.prefix, [file(params.reads1), file(params.reads2)]]])
	.set { reads_merge_paired_end_cleaned }

	//Stage boilerplate log
	merge_paired_end_cleaned_log = Channel.from(file("$baseDir/assets/merge_paired_end_cleaned_mqc.yaml"))

	//Initialise empty channels
	reads_profile_taxa = Channel.empty()
} else {
	reads_merge_paired_end_cleaned = Channel.empty()
	reads_profile_taxa = Channel.empty()
}

process merge_paired_end_cleaned {

	tag "$name"

	input:
	tuple val(name), file(reads) from reads_merge_paired_end_cleaned

	output:
	tuple val(name), path("*.fq.gz") into to_profile_taxa_merged

	when:
	(params.skip_preprocess && (!params.skip_profile_taxa || !params.skip_profile_function) && !params.singleEnd)

   	script:
	"""
	# This step will have no logging because the information are not relevant
	# I will simply use a boilerplate YAML to record that this has happened
	# If the files were not compressed, they will be at this stage
	if (file ${reads[0]} | grep -q compressed ) ; then
	    cat ${reads[0]} ${reads[1]} > ${name}.fq.gz
	else
		cat ${reads[0]} ${reads[1]} | gzip > ${name}.fq.gz
	fi
	"""
}

/**
	Community Characterisation - STEP 1. Performs taxonomic binning and estimates the
	microbial relative abundances using MetaPhlAn and its databases of clade-specific markers.
*/


// Defines channels for bowtie2_metaphlan_databases file
// log.info params.metaphlan_databases
// exit 1, params.metaphlan_databases
// bowtie2_metaphlan_databases = Channel.value( params.metaphlan_databases )
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
	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'link', pattern: "${name}_metaphlan_bugs_list.tsv", overwrite: true

	input:
	tuple val(name), file(reads) from to_profile_taxa_decontaminated.mix(to_profile_taxa_merged).mix(reads_profile_taxa)
	file(bowtie2db) from bowtie2_metaphlan_databases

	output:
	tuple val(name), path("*.biom") into to_alpha_diversity
	file("*_metaphlan_bugs_list.tsv") into to_collect_taxonomic_profiles
	tuple val(name), path(reads), path("*_metaphlan_bugs_list.tsv") into to_profile_function
	path "${name}_metaphlan_bugs_list.tsv" into profile_taxa_log

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	#If a file with the same name is already present, Metaphlan2 used to crash, leaving this here just
	# in case
	rm -rf ${name}_bt2out.txt

	metaphlan --input_type fastq --tmp_dir=. --biom ${name}.biom --bowtie2out=${name}_bt2out.txt \\
		--bowtie2db $bowtie2db --bt2_ps ${params.bt2options} --add_viruses --sample_id ${name} \\
		--nproc ${task.cpus} $reads ${name}_metaphlan_bugs_list.tsv &> profile_taxa_mqc.txt

	# MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_taxa_log.sh ${name}_metaphlan_bugs_list.tsv > ${name}.yaml
	"""
}

process collect_taxonomic_profiles {
	label 'biobakery'

	publishDir "${params.outdir}", mode: 'link', pattern: "merged_metaphlan_abundance_table.txt", overwrite: true

	input:
	file '*' from to_collect_taxonomic_profiles.collect()

	output:
	tuple file('merged_metaphlan_abundance_table.txt'), file('merged_metaphlan_abundance_table.species.txt') into to_hclust_taxonomic_profiles

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	# merge_metaphlan_tables.py will put the whole filename (without .tsv) as column names
	for i in *_metaphlan_bugs_list.tsv; do
		mv \$i \${i/_metaphlan_bugs_list/}
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
	tuple path(table_stratified), path(table_species) from to_hclust_taxonomic_profiles

	output:
	path '*.png'

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	#sed -i 's/\\t/ /g' ${table_stratified}

	hclust2.py -i ${table_stratified} -o metaphlan.top50.hclust2.png --skip_rows 0 \\
		--ftop 50 --f_dist_f cosine --s_dist_f braycurtis --cell_aspect_ratio 9 -s --fperc 99 \\
		--flabel_size 4 --legend_file metaphlan.hclust2.legend.png --max_flabel_len 100 \\
		--metadata_height 0.075 --minv 0.01 --no_slabels --dpi 300 --slinkage complete --no_fclustering

	# sed -i 's/\\t/ /g' ${table_species}

	hclust2.py -i ${table_species} -o metaphlan.species.top50.hclust2.png --skip_rows 0 \\
		--ftop 50 --f_dist_f cosine --s_dist_f braycurtis --cell_aspect_ratio 9 -s --fperc 99 \\
		--flabel_size 4 --legend_file metaphlan.species.hclust2.legend.png --max_flabel_len 100 \\
		--metadata_height 0.075 --minv 0.01 --no_slabels --dpi 300 --slinkage complete --no_fclustering
	"""
}

/**
	Community Characterisation - STEP 2. Performs the functional annotation using HUMAnN.
*/

// Defines channels for bowtie2_metaphlan_databases file
chocophlan_databases = file( params.chocophlan, type: 'dir', checkIfExists: true )
uniref_databases = file( params.uniref, type: 'dir', checkIfExists: true )

Channel
    .fromPath('reads/*')
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple()
    .set{ groups_ch }

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
	// FIXME the metaphlan bug list is A but the reads are B!!!
	tuple val(name), file(reads), file(metaphlan_bug_list) from to_profile_function
	file(chocophlan) from chocophlan_databases
	file(uniref) from uniref_databases

    output:
	file "*.HUMAnN.log" into profile_functions_log
	path "${name}_humann_temp/*"
	file "*_genefamilies.tsv" into humann_gene_families
	file "*_pathcoverage.tsv" into humann_path_coverage
	file "*_pathabundance.tsv" into humann_path_abundance

	when:
	!params.skip_profile_taxa && !params.skip_profile_function

	script:
	"""
	#HUMAnN will uses the list of species detected by the profile_taxa process

	humann --input $reads --output . --output-basename ${name} --taxonomic-profile $metaphlan_bug_list \
		--nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads \
		${task.cpus} --memory-use maximum &> ${name}.HUMAnN.log
	sleep 60
	bgzip ${name}_humann_temp/*.sam --compress-level 9 || true
	bgzip ${name}_humann_temp/*.tsv --compress-level 9 || true
	bgzip ${name}_humann_temp/*.fa --compress-level 9 || true
 	"""
}

process collect_functional_profiles {
	label 'biobakery'

	publishDir "${params.outdir}", mode: 'link', overwrite: true
	// "{all_genefamilies.tsv,all_pathcoverage.tsv,all_pathabundance.tsv,all_genefamilies-cpm.tsv,all_pathcoverage-cpm.tsv,all_pathabundance-cpm.tsv}"

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
	path('all_genefamilies-relab.tsv') into gene_fams_to_hclust
	path('all_pathabundance.tsv') into pathabunds_to_hclust
	path('all_pathcoverage.tsv') into pathcovs_to_hclust

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	humann_join_tables -i ./ -o all_genefamilies.tsv --file_name genefamilies
	humann_join_tables -i ./ -o all_pathcoverage.tsv --file_name pathcoverage
	humann_join_tables -i ./ -o all_pathabundance.tsv --file_name pathabundance

	humann_renorm_table -i all_genefamilies.tsv -o all_genefamilies-copm.tsv --units cpm
	humann_renorm_table -i all_genefamilies.tsv -o all_genefamilies-relab.tsv --units relab
	"""
}


process hclust_functional_profiles {
	label 'hclust'

	publishDir "${params.outdir}", mode: 'link', overwrite: true

	input:
	path gene_fam_table from gene_fams_to_hclust
	path path_abund_table from pathabunds_to_hclust
	path path_cov_table from pathcovs_to_hclust

	output:
	path '*.png'

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	# sed -i 's/\\t/ /g' ${gene_fam_table}
	# sed -i 's/\\t/ /g' ${path_abund_table}
	# sed -i 's/\\t/ /g' ${path_cov_table}

	if [[ \$(awk '{print NF}' ${gene_fam_table} | head -n 2 | tail -n 1) -gt 3 ]]; then
		cluster_columns=""
	else
		cluster_columns="--no_fclustering"
	fi

	hclust2.py -i ${gene_fam_table} -o all_genefamilies.top50.hclust2.png --skip_rows 0 \\
		--ftop 50 --f_dist_f cosine --s_dist_f braycurtis --cell_aspect_ratio 9 -s --fperc 99 \\
		--flabel_size 4 --legend_file all_genefamilies.top50.legend.png --max_flabel_len 100 \\
		--metadata_height 0.075 --minv 0.01 --no_slabels --dpi 300 --slinkage complete \$cluster_columns

	hclust2.py -i ${path_cov_table} -o all_pathcoverage.top50.hclust2.png --skip_rows 0 \\
		--ftop 50 --f_dist_f cosine --s_dist_f cityblock --cell_aspect_ratio 9 -s --fperc 99 \\
		--flabel_size 4 --legend_file all_pathcoverage.top50.legend.png --max_flabel_len 100 \\
		--metadata_height 0.075 --minv 0.01 --no_slabels --dpi 300 --slinkage complete \$cluster_columns

	hclust2.py -i ${path_abund_table} -o all_pathabundance.top50.hclust2.png \\
		--skip_rows 0 --ftop 50 --f_dist_f cosine --s_dist_f canberra --cell_aspect_ratio 9 -s \\
		--fperc 99 --flabel_size 4 --legend_file all_pathabundance.top50.legend.png \\
		--max_flabel_len 100 --metadata_height 0.075 --minv 0.01 --no_slabels --dpi 300 \\
		--slinkage complete \$cluster_columns
	"""
}


/**
	Community Characterisation - STEP 3. Evaluates several alpha-diversity measures.

*/

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
	file "${name}_alpha-diversity.tsv" into alpha_diversity_log

	when:
	!params.skip_profile_taxa || !params.skip_profile_function

	script:
	"""
	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $metaphlan_bug_list | wc -l  | cut -d\" \" -f 1)
	if (( n <= 3 )); then
		#The file should be created in order to be returned
		touch ${name}_alpha-diversity.tsv
	else
		echo $name > ${name}_alpha-diversity.tsv
		qiime tools import --input-path $metaphlan_bug_list --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ${name}_abundance_table.qza
		for alpha in ace berger_parker_d brillouin_d chao1 chao1_ci dominance doubles enspie esty_ci fisher_alpha gini_index goods_coverage heip_e kempton_taylor_q lladser_pe margalef mcintosh_d mcintosh_e menhinick michaelis_menten_fit osd pielou_e robbins shannon simpson simpson_e singles strong
		do
			qiime diversity alpha --i-table ${name}_abundance_table.qza --p-metric \$alpha --output-dir \$alpha &> /dev/null
			qiime tools export --input-path \$alpha/alpha_diversity.qza --output-path \${alpha} &> /dev/null
			if [[ -e "\${alpha}/alpha-diversity.tsv" ]]; then
                value=\$(sed -n '2p' \${alpha}/alpha-diversity.tsv | cut -f 2)
            else
                value="#N/A"
            fi
		    echo -e  \$alpha'\t'\$value
		done >> ${name}_alpha-diversity.tsv
	fi
	"""
}


// ------------------------------------------------------------------------------
//	MULTIQC LOGGING
// ------------------------------------------------------------------------------


/**
	Generate Logs.

	Logs generate at each analysis step are collected and processed with MultiQC
*/
if(params.skip_preprocess) {
	Channel.fromPath("${params.indir}/*/ANALYSIS/01_dedup_log.txt")
		.set { dedup_log }
		
	Channel.fromPath("${params.indir}/*/ANALYSIS/02_syndecontam_log.txt")
		.set { synthetic_contaminants_log }
	Channel.fromPath("${params.indir}/*/ANALYSIS/03_trim_log.txt")
		.set { trimming_log }
}
if(params.skip_quality_assessment) {
	// find multiQC inputs from files
	Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/raw/*")
		.set { fastqc_raw_log }
	Channel.fromPath("${params.indir}/*/ANALYSIS/fastqc/qcd/*")
		.set { fastqc_qcd_log }
}
if(params.skip_profile_taxa) {
	Channel.fromPath("${params.indir}/*/ANALYSIS/*_metaphlan_bugs_list.tsv")
		.set { profile_taxa_log }
	Channel.fromPath("${params.indir}/*/ANALYSIS/*_alpha-diversity.tsv")
		.set {alpha_diversity_log}
	Channel.fromPath("${params.indir}/*/ANALYSIS/*.HUMAnN.log")
		.set { profile_functions_log }
}

process log {

	publishDir "${params.outdir}", mode: 'link', overwrite: true

    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_multiqc
    } else {
        container params.docker_container_multiqc
    }

	input:
	file multiqc_config from file(params.multiqc_config, type: 'file', checkIfExists: true )
	path "yamp.css" from file(params.multiqc_css, type: 'file', checkIfExists: true )
	path "replace_sample_names.txt" from file(params.multiqc_replace_names, type: 'file', checkIfExists: true )
	file workflow_summary from create_workflow_summary(summary)
	file "software_versions_mqc.yaml" from software_versions_yaml.collect()
	path "fastqc_raw/*" from fastqc_raw_log.collect().ifEmpty([])
	file "dedup/log*.txt" from dedup_log.collect().ifEmpty([])
	file "syncontam/log*.txt" from synthetic_contaminants_log.collect().ifEmpty([])
	file "trimming/log*.txt" from trimming_log.collect().ifEmpty([])
// 	file "foreign_genome_indexing_mqc.yaml" from index_foreign_genome_log.collect()
	path "fastqc_QCd/*" from fastqc_qcd_log.collect().ifEmpty([])
	path "metaphlan/*" from profile_taxa_log.collect().ifEmpty([])
		//atm we do nothing with the qiime output
	path "qiime/*" from alpha_diversity_log.collect().ifEmpty([])
	file "humann/*" from profile_functions_log.collect().ifEmpty([])

	output:
	path "*multiqc_report*.html"
	path "*multiqc_data*"

	script:
	"""
	printf "\\traw\\tduplicated\\tdeduped\\ttime\\n" > dedup_data.txt # or copy header from one of the files
	for f in dedup/log*.txt
		do
		totR=\$(grep "Reads In:" "\${f}" | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
		remR=\$(grep "Duplicates Found:" "\${f}" | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
		survivedR=\$((\$totR - \$remR))
		percentage=\$(echo \$survivedR \$totR | awk '{print \$1/\$2*100}' )
		percentage=`printf "%.2f" \$percentage`
		time=\$(grep "Total time:" "\${f}" | cut -d: -f 2 | cut -f 2 | sed 's/s\\./s/g')
		samplename=\$(head -n 1 "\${f}" | grep -o "\\sin=[^\\S\\.]*" | sed 's/[[:blank:]]in=//' )
		printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\${samplename}" "\${totR}" "\${remR}" "\${survivedR}" "\${time}" >> dedup_data.txt
	done

	printf "\\tsynDecontam\\ttrimmed\\n" > bbduk_data.txt
	for f in syncontam/log*.txt
	 	do
		samplename=\$(head -n 1 \${f} | grep -o "\\sin=[^\\S\\.]*" | sed 's/[[:blank:]]in=//' )
		totR=\$(grep "Input:" \${f} | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
		remR=\$(grep "Contaminants:" \${f} | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
		reads=\$((\$totR-\$remR))
	 	printf "%s\\t%s\\n" "\${samplename}" "\${reads}" >> bbduk_data.txt
	done

	for f in trimming/log*.txt
	 	do
		samplename=\$(head -n 1 \${f} | grep -o "\\sin=[^\\S\\.]*" | sed 's/[[:blank:]]in=//' )
		totR=\$(grep "Input:" \${f} | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
		remR=\$(grep "Total Removed:" \${f} | cut -d: -f 2 | cut -f 2 | cut -d" " -f 1 | sed 's/ //g')
		reads=\$((\$totR-\$remR))
	 	sedstr="s/(\${samplename}.*)\$/\\1\\t\${reads}/"
	 	sed -i -E "\${sedstr}" bbduk_data.txt
	done

	printf "\\tdecontam\\n" > decontam_data.txt
	for f in fastqc_QCd/*.html
		do
		reads=\$(grep -o "<tr><td>Total Sequences</td><td>[0-9]*</td></tr>" "\${f}" | grep -o "[0-9]*")
		filename=\${f##*/}
		samplename=\${filename%_fastqc.html}
		printf "%s\\t%s\\n" "\${samplename}" "\${reads}" >> decontam_data.txt
	done

	printf "\\tnSpecies\\tperc_explained\\n" > profile_taxa_data.txt
	for f in humann/*.log
		do
		filename=\${f##*/}
		ext=\${filename#*.}
		samplename=\${filename%.\$ext}
		tot_species_prescreening=\$(grep "Total species selected from prescreen:" "\${f}" | cut -d: -f 2 | sed 's/ //g')
		selected_species_explain=\$(grep "Selected species explain" "\${f}" | cut -d" " -f 4 | grep -o "[0-9\\.]*")
		printf "%s\\t%s\\t%s" \$samplename \${tot_species_prescreening} \${selected_species_explain} >> profile_taxa_data.txt
	done

	printf "\\taligned_nucleotide\\tunaligned_nucleotide_percent\\t" > profile_functions_data.txt
	printf "aligned_translated\\tperc_explained\\tgene_families\\n" >> profile_functions_data.txt
	for f in humann/*.log
		do
		filename=\${f##*/}
		ext=\${filename#*.}
		samplename=\${filename%.\$ext}
		tot_reads_aligned_nucleotide=\$(grep -zo "Total bugs from nucleotide alignment.*Total gene families from nucleotide alignment" "\${f}" | grep -Eoa "[0-9]+ hits" | grep -Eo "[0-9]+" | paste -sd+ | bc | sed 's/ //g')
		unaligned_reads_nucleotide=\$(grep -o "Unaligned reads after nucleotide alignment: [0-9\\.]*" "\${f}" | grep -o "[0-9\\.]*")
		tot_reads_aligned_both=\$(grep -zo "Total bugs after translated alignment.*Total gene families after translated alignment" "\${f}" | grep -Eoa "[0-9]+ hits" | grep -Eo "[0-9]+" | paste -sd+ | bc | sed 's/ //g')
		tot_reads_aligned_translated=\$((\$tot_reads_aligned_both - \$tot_reads_aligned_nucleotide))
		unaligned_reads_translated=\$(grep "Unaligned reads after translated alignment:" "\${f}" | grep -o "[0-9\\.]*")
		perc_explained=\$(echo "100 - \$unaligned_reads_translated" | bc)
		tot_gene_family=\$(grep "Total gene families after translated alignment:" "\${f}" | cut -d: -f 2 | sed 's/ //g')
		printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \$samplename \${tot_reads_aligned_nucleotide} \${unaligned_reads_nucleotide} \${tot_reads_aligned_translated} \${perc_explained} \${tot_gene_family} >> profile_functions_data.txt
	done

	multiqc --config $multiqc_config . -f --custom-css-file yamp.css
    """
}



