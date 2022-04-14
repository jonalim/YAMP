#!/usr/bin/env nextflow
	
/**
Yet Another Metagenomic Pipeline (YAMP)
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
	
For any bugs or problems found, please go to:
- https://github.com/alesssia/YAMP/issues
*/


def versionMessage() 
{
	log.info"""
	 
	YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: ${workflow.manifest.version} 
	""".stripIndent()
}

def helpMessage() 
{
	log.info"""

YET ANOTHER METAGENOMIC PIPELINE (YAMP) - Version: ${workflow.manifest.version} 

This pipeline is distributed in the hope that it will be useful
but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.

Please report comments and bugs at https://github.com/alesssia/YAMP/issues.
Check https://github.com/alesssia/YAMP for updates, and refer to
https://github.com/alesssia/YAMP/wiki for more details.

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

YAMP supports FASTQ and compressed FASTQ files.
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
if (params.mode != "QC" && params.mode != "characterisation" && params.mode != "complete") {
	exit 1, "Mode not available. Choose any of <QC, characterisation, complete>"
}	

if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
if (params.mode != "characterisation" && !params.singleEnd && (params.reads2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
}

//--reads1 and --reads2 can be omitted (and the default from the config file used instead) 
//only when mode is "characterisation". Obviously, --reads2 should be always omitted when the
//library layout is single.
if (params.mode != "characterisation" && ( (!params.singleEnd && (params.reads1 == "null" || params.reads2 == "null")) || (params.singleEnd && params.reads1 == "null")) ) {
	exit 1, "Please set the reads1 and/or reads2 parameters"
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
YET ANOTHER METAGENOMIC PIPELINE (YAMP) 
---------------------------------------------

Analysis introspection:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date() 
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'YAMP'
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

if (params.mode != "characterisation") 
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

if (params.mode != "QC")
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
summary['Running parameters'] = ""
summary['Prefix'] = params.prefix
summary['Running mode'] = params.mode
summary['Layout'] = params.singleEnd ? 'Single-End' : 'Paired-End'

if (params.mode != "characterisation") 
{
	if (!params.singleEnd) 
	{
		summary['QC mate pairs separately'] = params.qc_matepairs
	}
	summary['Performing de-duplication'] = params.dedup

	//remove_synthetic_contaminants 
	summary['Synthetic contaminants'] = ""
	summary['Artefacts'] = params.artefacts
	summary['Phix174ill'] = params.phix174ill

	//Trimming
	summary['Adapters'] = params.adapters
	summary['Trimming parameters'] = ""
	summary['Input quality offset'] = params.qin == 33 ? 'ASCII+33' : 'ASCII+64'
	summary['Min phred score'] = params.phred
	summary['Min length'] = params.minlength
	summary['kmer lenght'] = params.kcontaminants
	summary['Shorter kmer'] = params.mink 
	summary['Max Hamming distance'] = params.hdist 

	//Decontamination
	summary['Decontamination parameters'] = ""
	if (params.foreign_genome_ref != "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome_ref + " (indexed)"
	} else if (	params.foreign_genome_ref == "") {
		summary['Contaminant (pan)genome'] = params.foreign_genome
	}	
	summary['Min alignment identity'] = params.mind
	summary['Max indel length'] = params.maxindel
	summary['Max alignment band'] = params.bwr
}

if (params.mode != "QC")
{
    //BowTie2 databases for metaphlan
	summary['MetaPhlAn parameters'] = ""
    summary['MetaPhlAn database'] = params.metaphlan_databases
    summary['Bowtie2 options'] = params.bt2options
  
    // ChocoPhlAn and UniRef databases
	summary['HUMAnN parameters'] = ""
	summary['Chocophlan database'] = params.chocophlan
	summary['Uniref database'] = params.uniref
}

//Folders
summary['Folders'] = ""
summary['Input directory'] = params.reads
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
    section_name: 'YAMP Workflow Summary'
    section_href: 'https://github.com/alesssia/yamp'
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
if (params.mode != "characterisation" && params.foreign_genome_ref == "") {
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
	params.mode != "characterisation" && params.foreign_genome_ref == ""

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
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants }
} else {
	if(params.qc_matepairs) {
		Channel.fromFilePairs("${params.reads}/*/ILLUMINA_DATA/*_R{1,2}.fastq.gz", size: -1) { file -> file.getParent().getParent().getName() }
			.into {read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants; read_files_log }
		// Channel
		// .fromFilePairs("${params.reads}/*_R{1,2}.fastq.gz", checkIfExists: true)
		// //.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
		// .into { read_files_fastqc; read_files_dedup; read_files_synthetic_contaminants; read_files_log }
	} else {
		Channel.fromFilePairs("${params.reads}/*/ILLUMINA_DATA/*_R{1,2}.fastq.gz", size: -1) { file -> file.getParent().getParent().getName() }
			.set { to_combine_reads }
		// Channel.fromFilePairs("${params.reads}/*_R{1,2}.fastq.gz", checkIfExists: true)
		// .set {to_combine_reads}
	}
}

// Defines channels for resources file 
artefacts = file(params.artefacts, type: "file", checkIfExists: true )
phix174ill = file(params.phix174ill, type: "file", checkIfExists: true)
adapters = file(params.adapters, type: "file", checkIfExists: true)

process preprocess {
	tag "$name"

	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'copy', pattern: "decontam/*.fq.gz", saveAs: {filename -> file(filename).getName()}

	input:
	tuple val(name), file(reads) from to_combine_reads
	file(artefacts) from artefacts
	file(phix174ill) from phix174ill
	file(adapters) from adapters
	path(ref_foreign_genome) from ref_foreign_genome

	output:
	tuple val(name), path("${name}.fq.gz") into read_files_fastqc, read_files_log
	file "dedup/*.yaml" into dedup_log
	file "syndecontam/*.yaml" into synthetic_contaminants_log
	file "trim/*.yaml" into trimming_log
	tuple val(name), path("decontam/*_QCd.fq.gz") into qcd_reads
	tuple val(name), path("decontam/*_QCd.fq.gz") into to_profile_taxa_decontaminated
	file "decontam/*.yaml" into decontaminate_log

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
    clumpify.sh -Xmx\"\$maxmem\" in=${name}.fq.gz out=\"${name}_dedup.fq.gz\" qin=$params.qin dedupe reorder subs=0 threads=${task.cpus} &> dedup_mqc.txt

	# MultiQC doesn't have a module for clumpify yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	mkdir dedup
	bash scrape_dedup_log.sh > dedup/${name}.yaml

	# remove synthetic contaminants
	bbduk.sh -Xmx\"\$maxmem\" in=\"${name}_dedup.fq.gz\" out=\"${name}_no_synthetic_contaminants.fq.gz\" k=31 ref=$phix174ill,$artefacts qin=$params.qin ordered=t threads=${task.cpus} ow &> synthetic_contaminants_mqc.txt
	rm \"${name}_dedup.fq.gz\"
	# rm \"${name}_dedup_R2.fq.gz\"

	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	mkdir syndecontam
	bash scrape_remove_synthetic_contaminants_log.sh > syndecontam/${name}.yaml

	# trim
	bbduk.sh -Xmx\"\$maxmem\" in=\"${name}_no_synthetic_contaminants.fq.gz\" out=\"${name}_trimmed.fq.gz\" ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred  minlength=$params.minlength ref=$adapters qin=$params.qin ordered=t threads=${task.cpus} tbo tpe ow &> trimming_mqc.txt
	rm \"${name}_no_synthetic_contaminants.fq.gz\"
	# rm \"${name}_no_synthetic_contaminants_R2.fq.gz\"

	# MultiQC doesn't have a module for bbduk yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	mkdir trim
	bash scrape_trimming_log.sh > trim/${name}.yaml

	# decontaminate
	mkdir decontam

	bbwrap.sh -Xmx\"\$maxmem\" mapper=bbmap append=t in=\"${name}_trimmed.fq.gz\" outu=\"decontam/${name}_QCd.fq.gz\" outm=\"decontam/${name}_contamination.fq.gz\" minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred path="./" qin=$params.qin threads=${task.cpus} untrim quickmatch fast ordered=t ow &> decontamination_mqc.txt
	rm \"${name}_trimmed.fq.gz\"
	#rm \"${name}_trimmed_singletons.fq.gz\"
	# rm \"${name}_trimmed_R2.fq.gz\"

	# MultiQC doesn't have a module for bbwrap yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_decontamination_log.sh \"decontam/${name}_QCd.fq.gz\" \"decontam/${name}_contamination.fq.gz\" > decontam/${name}.yaml
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
// 	// params.mode != "characterisation" && params.dedup

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
// if (!params.dedup & params.mode != "characterisation") {
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
// 	params.mode != "characterisation"

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
// 	params.mode != "characterisation"

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

// 	publishDir "${params.outdir}/${name}", mode: 'copy', pattern: "*QCd.fq.gz"

// 	input:
// 	tuple val(name), file(reads) from to_decontaminate
// 	path(ref_foreign_genome) from ref_foreign_genome

// 	output:
// 	tuple val(name), path("decontam/*.fq.gz") into qcd_reads
// 	tuple val(name), path("decontam/*.fq.gz") into to_profile_taxa_decontaminated
// 	file "*.yaml" into decontaminate_log

// 	when:
// 	params.mode != "characterisation"

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
	
	publishDir "${params.outdir}/${name}/ANALYSIS/fastqc", mode: 'copy' //,

    input:
    tuple val(name), file(reads: 'raw/*'), file(qcd_reads: "qcd/${name}.fq.gz") from read_files_fastqc.join(qcd_reads, failOnDuplicate: true, failOnMismatch: true)
	
    output:
    path "raw/*_fastqc.{zip,html}" into fastqc_sample_log, fastqc_raw_project_log
	path "qcd/*_fastqc.{zip,html}" into fastqc_qcd_project_log
	
	when:
	params.mode != "characterisation"

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
if (params.mode == "characterisation" && params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1)]]])
	.into { reads_profile_taxa }
	
	//Initialise empty channels
	reads_merge_paired_end_cleaned = Channel.empty()
	//Init as value channel so it can be reused by log process
	merge_paired_end_cleaned_log = Channel.value([])
} else if (params.mode == "characterisation" && !params.singleEnd) {
	Channel
	.from([[params.prefix, [file(params.reads1), file(params.reads2)]]] )
	.set { reads_merge_paired_end_cleaned }
	
	//Stage boilerplate log
	merge_paired_end_cleaned_log = Channel.from(file("$baseDir/assets/merge_paired_end_cleaned_mqc.yaml"))
	
	//Initialise empty channels
	reads_profile_taxa = Channel.empty()
} else if (params.mode != "characterisation")
{
	//Initialise empty channels
	reads_merge_paired_end_cleaned = Channel.empty()
	reads_profile_taxa = Channel.empty()
	//Init as value channel so it can be reused by log process
	merge_paired_end_cleaned_log = Channel.value([])
}

process merge_paired_end_cleaned {

	tag "$name"
		
	input:
	tuple val(name), file(reads) from reads_merge_paired_end_cleaned
	
	output:
	tuple val(name), path("*.fq.gz") into to_profile_taxa_merged
	
	when:
	params.mode == "characterisation" && !params.singleEnd

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

    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'copy', pattern: "*.{biom,tsv}"
	
	input:
	tuple val(name), file(reads) from to_profile_taxa_decontaminated.mix(to_profile_taxa_merged).mix(reads_profile_taxa)
	file(bowtie2db) from bowtie2_metaphlan_databases

	output:
	tuple val(name), path("*.biom") into to_alpha_diversity
	file("*_metaphlan_bugs_list.tsv") into to_collect_taxonomic_profiles
	tuple val(name), path(reads), path("*_metaphlan_bugs_list.tsv") into to_profile_function
	tuple val(name), path("${name}.yaml") into profile_taxa_sample_log
	path "${name}.yaml" into profile_taxa_project_log
	
	when:
	params.mode != "QC"
	
	script:
	"""
	#If a file with the same name is already present, Metaphlan2 used to crash, leaving this here just in case
	rm -rf ${name}_bt2out.txt

	metaphlan --input_type fastq --tmp_dir=. --biom ${name}.biom --bowtie2out=${name}_bt2out.txt --bowtie2db $bowtie2db --bt2_ps ${params.bt2options} --add_viruses --sample_id ${name} --nproc ${task.cpus} $reads ${name}_metaphlan_bugs_list.tsv &> profile_taxa_mqc.txt
	
	# MultiQC doesn't have a module for Metaphlan yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_taxa_log.sh ${name}_metaphlan_bugs_list.tsv > ${name}.yaml
	"""
}

process collect_taxonomic_profiles {
	publishDir "${params.outdir}", mode: 'copy', pattern: "merged_metaphlan_abundance_table.txt"

	input:
	file '*' from to_collect_taxonomic_profiles.collect()

	output:
	file 'merged_metaphlan_abundance_table.txt' into to_hclust

	when:
	params.mode != "QC"

	script:
	"""
	# merge_metaphlan_tables.py will put the whole filename (without .tsv) as column names
	for i in *_metaphlan_bugs_list.tsv; do
		mv \$i \${i/_metaphlan_bugs_list/} 
	done

	merge_metaphlan_tables.py *.tsv > merged_metaphlan_abundance_table.txt
	"""
}

process hclust_taxonomic_profiles {
	publishDir "${params.outdir}", mode: 'copy'

	input:
	file 'table' from to_hclust

	output:
	path '*.png'

	when:
	params.mode != "QC"

	script:
	"""
	awk '{\$2=""; print \$0}' ${table} > table.txt

	hclust2.py \
		-i table.txt \
		-o metaphlan.hclust2.sqrt_scale.png \
		--sep '\\s+' \
		--skip_rows 0 \
		--ftop 50 \
		--f_dist_f correlation \
		--s_dist_f braycurtis \
		--cell_aspect_ratio 9 \
		-s --fperc 99 \
		--flabel_size 4 \
		--legend_file metaphlan.hclust2.sqrt_scale.legend.png \
		--max_flabel_len 100 \
		--metadata_height 0.075 \
		--minv 0.01 \
		--no_slabels \
		--dpi 300 \
		--slinkage complete \
		--no_fclustering # don't know why, but fclustering doesn't work. Might be bc only 2 samples?
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
	
    tag "$name"

	//Enable multicontainer settings
    if (workflow.containerEngine == 'singularity') {
        container params.singularity_container_biobakery
    } else {
        container params.docker_container_biobakery
    }

	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'copy', pattern: "*.{tsv,log}"
	publishDir "${params.outdir}/${name}/ANALYSIS/humann_intermediate", mode: 'move', pattern: "${name}_humann_temp/*", saveAs: {filename -> 
		f = new File(filename);
		f.getName()
	}
	
	input:
	// FIXME the metaphlan bug list is A but the reads are B!!!
	tuple val(name), file(reads), file(metaphlan_bug_list) from to_profile_function
	file(chocophlan) from chocophlan_databases
	file(uniref) from uniref_databases

    output:
	file "*_HUMAnN.log"
	path "${name}_humann_temp/*"
	file "*_genefamilies.tsv" into humann_gene_families
	file "*_pathcoverage.tsv" into humann_path_coverage
	file "*_pathabundance.tsv" into humann_path_abundance
	file "*.yaml" into profile_functions_log, profile_functions_project_log

	when:
	params.mode != "QC"

	script:
	"""
	#HUMAnN will uses the list of species detected by the profile_taxa process

	humann --input $reads --output . --output-basename ${name} --taxonomic-profile $metaphlan_bug_list --nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads ${task.cpus} --memory-use maximum &> ${name}_HUMAnN.log 
	bgzip ${name}_humann_temp/*.sam --compress-level 9
	bgzip ${name}_humann_temp/*.fa --compress-level 9
	bgzip ${name}_humann_temp/*.tsv --compress-level 9

	# MultiQC doesn't have a module for humann yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash scrape_profile_functions.sh ${name} ${name}_HUMAnN.log ${reads} > ${name}.yaml
 	"""
}

process collect_functional_profiles {
	publishDir "${params.outdir}", mode: 'copy'
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

	when:
	params.mode != "QC"

	script:
	"""
	humann_join_tables -i ./ -o all_genefamilies.tsv --file_name genefamilies
	humann_join_tables -i ./ -o all_pathcoverage.tsv --file_name pathcoverage
	humann_join_tables -i ./ -o all_pathabundance.tsv --file_name pathabundance

	humann_renorm_table -i all_genefamilies.tsv -o all_genefamilies-cpm.tsv --units cpm
	humann_renorm_table -i all_pathcoverage.tsv -o all_pathcoverage-cpm.tsv --units cpm
	humann_renorm_table -i all_pathabundance.tsv -o all_pathabundance-cpm.tsv --units cpm
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

	publishDir "${params.outdir}/${name}/ANALYSIS/", mode: 'copy', pattern: "*.{tsv}"
	
	input:
	tuple val(name), file(metaphlan_bug_list) from to_alpha_diversity
		
    output:
	file "${name}.tsv" into alpha_diversity_project_log
	file "${name}.yaml" into alpha_diversity_sample_log
	
	when:
	params.mode != "QC"

	script:
	"""
	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $metaphlan_bug_list | wc -l  | cut -d\" \" -f 1)
	if (( n <= 3 )); then
		#The file should be created in order to be returned
		touch ${name}.tsv 
	else
		echo $name > ${name}.tsv
		qiime tools import --input-path $metaphlan_bug_list --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ${name}_abundance_table.qza
		for alpha in ace berger_parker_d brillouin_d chao1 chao1_ci dominance doubles enspie esty_ci fisher_alpha gini_index goods_coverage heip_e kempton_taylor_q lladser_pe margalef mcintosh_d mcintosh_e menhinick michaelis_menten_fit osd pielou_e robbins shannon simpson simpson_e singles strong
		do
			qiime diversity alpha --i-table ${name}_abundance_table.qza --p-metric \$alpha --output-dir \$alpha &> /dev/null
			qiime tools export --input-path \$alpha/alpha_diversity.qza --output-path \${alpha} &> /dev/null
			value=\$(sed -n '2p' \${alpha}/alpha-diversity.tsv | cut -f 2)
		    echo -e  \$alpha'\t'\$value 
		done >> ${name}.tsv  
	fi

	# MultiQC doesn't have a module for qiime yet. As a consequence, I
	# had to create a YAML file with all the info I need via a bash script
	bash generate_alpha_diversity_log.sh \${n} > ${name}.yaml	
	"""
}


// ------------------------------------------------------------------------------   
//	MULTIQC LOGGING
// ------------------------------------------------------------------------------   


/**
	Generate Logs. 

	Logs generate at each analysis step are collected and processed with MultiQC 
*/

// Stage config files
// read_files_log.view{ x -> "read_files_log contains: $x" }
// fastqc_log.view{ x -> "fastqc_log contains: $x" }
// dedup_log.view{ x -> "dedup_log contains: $x" }
// synthetic_contaminants_log.view{ x -> "synthetic_contaminants_log contains: $x" }
// trimming_log.view{ x -> "trimming_log contains: $x" }

// process log {
	
// 	publishDir "${params.outdir}/${name}", mode: 'copy'

//     if (workflow.containerEngine == 'singularity') {
//         container params.singularity_container_multiqc
//     } else {
//         container params.docker_container_multiqc
//     }

// 	input:
// 	file multiqc_config from file(params.multiqc_config, type: 'file', checkIfExists: true )
// 	path "yamp.css" from file(params.multiqc_css, type: 'file', checkIfExists: true )
// 	tuple val(name), file(reads) from read_files_log
// 	file workflow_summary from create_workflow_summary(summary)
// 	file "software_versions_mqc.yaml" from software_versions_yaml.collect()
// 	path "fastqc/*" from fastqc_sample_log
// 	file "dedup_mqc.yaml" from dedup_log
// 	file "synthetic_contaminants_mqc.yaml" from synthetic_contaminants_log
// 	file "trimming_mqc.yaml" from trimming_log
// 	file "foreign_genome_indexing_mqc.yaml" from index_foreign_genome_log.collect()
// 	file "decontamination_mqc.yaml" from decontaminate_log
// 	file profile_taxa_mqc from profile_taxa_sample_log
// 	file "merge_paired_end_cleaned_mqc.yaml" from merge_paired_end_cleaned_log
// 	file "profile_functions_mqc.yaml" from profile_functions_log
// 	file "alpha_diversity_mqc.yaml" from alpha_diversity_sample_log

// 	output:
// 	path "*multiqc_report*.html" into multiqc_report
// 	path "*multiqc_data*"

// 	script:
// 	"""
// 	multiqc --config $multiqc_config . -f --custom-css-file yamp.css
// 	mv multiqc_report.html ${name}_multiqc_report_${params.mode}.html
// 	mv multiqc_data ${name}_multiqc_data_${params.mode}
// 	"""
// }


process project_log {
	
	publishDir "${params.outdir}", mode: 'copy'

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
	path "fastqc_raw/*" from fastqc_raw_project_log.collect()
	file "dedup/*" from dedup_log.collect()
	file "syncontam/*" from synthetic_contaminants_log.collect()
	file "trimming/*" from trimming_log.collect()
// 	file "foreign_genome_indexing_mqc.yaml" from index_foreign_genome_log.collect()
	file "decontam/*" from decontaminate_log.collect()
	path "fastqc_QCd/*" from fastqc_qcd_project_log.collect()
	path "metaphlan/*" from profile_taxa_project_log.collect()
		//atm we do nothing with the qiime output
	path "qiime/*" from alpha_diversity_project_log.collect()
	file "humann/*" from profile_functions_project_log.collect()

	output:
	path "*multiqc_report*.html"
	path "*multiqc_data*"

	script:
	"""
	printf "\\traw\\tdeduped\\n" > dedup_data.txt
	for f in dedup/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		raw=\$(grep "<dt>Input:</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g' | grep -o "^[0-9]*")
		surviving=\$(grep "<dt>Surviving:</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g' | grep -o "^[0-9]*")
		printf "%s\\t%s\\n" "\${samplename}" "\${raw}" "\${surviving}" >> dedup_data.txt
	done

	printf "\\tsynDecontam\\ttrimmed\\n" > bbduk_data.txt
	for f in syncontam/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		reads=\$(grep "<dt>Surviving:</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g' | grep -o "^[0-9]*")
		printf "%s\\t%s\\n" "\${samplename}" "\${reads}" >> bbduk_data.txt
	done

	for f in trimming/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		reads=\$(grep "<dt>Surviving:</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g' | grep -o "^[0-9]*")
		sedstr="s/(\${samplename}.*)\$/\\1\\t\${reads}/"
		sed -i -E "\${sedstr}" bbduk_data.txt
	done

	printf "\\tdecontam\\n" > decontam_data.txt
	for f in decontam/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		reads=\$(grep "<dt>Surviving:</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g' | grep -o "^[0-9]*")
		printf "%s\\t%s\\n" "\${samplename}" "\${reads}" >> decontam_data.txt
	done

	printf "\\tnSpecies\\n" > profile_taxa_data.txt
	for f in metaphlan/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		species=\$(grep "<dt>Species</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g')
		printf "%s\\t%s\\n" "\${samplename}" "\${species}" >> profile_taxa_data.txt
	done

	printf "\\tpercExplained\\tgeneFamilies\\n" > profile_functions_data.txt
	for f in humann/*.yaml
		do
		filename=\${f##*/}
		samplename=\${filename%.yaml}
		species=\$(grep "<dt>Input after QC</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's/<\\/dd>.*\$//g')
		explained=\$(grep "<dt>Selected species explain</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed "s/% of QC'd reads<\\/dd>.*\$//g")
		genes=\$(grep "<dt>Total gene families</dt><dd>" \${f} | sed 's/^.*<dd>//g' | sed 's| (after translated alignment)</dd>||g')
		printf "%s\\t%s\\t%s\\n" "\${samplename}" "\${explained}" "\${genes}" >> profile_functions_data.txt
	done

	multiqc --config $multiqc_config . -f --custom-css-file yamp.css
    """
	
	// """
	// multiqc --config $multiqc_config . -f --custom-css-file yamp.css
	// mv multiqc_report.html ${name}_multiqc_report_${params.mode}.html
	// mv multiqc_data ${name}_multiqc_data_${params.mode}
	// """
}



