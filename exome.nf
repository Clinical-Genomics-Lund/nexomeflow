#!/usr/bin/env nextflow

params.fasta = "/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
params.bed = "/data/bnf/ref/agilent_clinical_research_exome_v2/S30409818_Regions.nochr.bed"
index = file(params.csv)


genome_file = file(params.fasta)
regions_bed = file(params.bed)


Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
    .set { fastq }


process bwa_align {
    cpus 20
    // publishDir '/trannel/proj/cmd-pipe/test-files/bam', mode: 'copy', overwrite: 'false'

    input: 
	set group, id, file(read1), file(read2) from fastq

    output:
    set group, id, file("${id}_bwa.sort.bam") into bwa_bam

    script:
	"""
    echo "$read1+$read2" > ${id}_bwa.sort.bam
    """
    //    /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon bwa mem -M -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -t ${task.cpus} $genome_file $read1 $read2 | /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon util sort -r $genome_file -o ${id}_bwa.sort.bam -t ${task.cpus} --sam2bam -i -

}

process markdup {
    cpus 10
    input:
	set group, id, file(sorted_bam) from bwa_bam

    output:
	set group, val(id), file("${id}.markdup.bam") into bam_markdup
    
    script:
    """
    wc -l $sorted_bam > ${id}.markdup.bam
    """
    //    sambamba markdup --tmpdir /data/tmp -t ${task.cpus} $sorted_bam ${id}.markdup.bam

}

bam_markdup.into {
    bam_qc
    bam_marked
}


process post_align_qc {
    cpus 6

    input:
    set group, id, file(markdup_bam) from bam_qc
    output:
    set id, file("${id}.bwa.qc")


    script:
    """
    echo "QC DONE" > ${id}.bwa.qc
    """
  //  /data/bnf/scripts/postaln_qc.pl $markdup_bam $regions_bed ${id} ${task.cpus} $regions_bed $genome_file > ${id}.bwa.qc


}



process gvcf_template {
    cpus 1

    publishDir '/trannel/proj/cmd-pipe/test-files', mode: 'copy', overwrite: 'true'
    input:
    set group, id, file(bam) from bam_marked.groupTuple()

    output:
    file("${group}.concat.count") into results

    script:
    """
    echo "${bam}" >> ${group}.concat.count
    """



}