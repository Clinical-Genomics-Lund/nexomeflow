#!/usr/bin/env nextflow

params.fasta = "/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
params.bed = "/data/bnf/ref/agilent_clinical_research_exome_v2/S30409818_Regions.nochr.bed"
index = file(params.csv)


genome_file = file(params.fasta)
regions_bed = file(params.bed)


Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
    .set { fastq }


process bwa_align {
    cpus 20
    
    input: 
	set id, file(read1), file(read2) from fastq

    output:
    set id, file("${id}_bwa.sort.bam")

    script:
	"""
    /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon bwa mem -M -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -t ${task.cpus} $genome_file $read1 $read2 | /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon util sort -r $genome_file -o ${id}_bwa.sort.bam -t ${task.cpus} --sam2bam -i -
    """

}
