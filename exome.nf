#!/usr/bin/env nextflow

params.fasta = "/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
params.bed = "/data/bnf/ref/agilent_clinical_research_exome_v2/S30409818_Regions.nochr.bed"

OUTDIR = file(params.outdir)

// might move these to config
genome_file = file(params.fasta)
regions_bed = file(params.bed)


Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
//    .map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
    .set { fastq }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype) }
    .set { ped }



// alignment using sentieon-bwa
process bwa_align {
    cpus 20
    // publishDir '/trannel/proj/cmd-pipe/test-files/bam', mode: 'copy', overwrite: 'false'

    input: 
	//set group, id, file(read1), file(read2) from fastq
    set clarity_sample_id, id, assay, sex, diagnosis, phenotype, group, father, mother, clarity_pool_id, platform, file(read1), file(read2), analysis_dir from fastq

    output:
    set group, id, analysis_dir, file("${id}_bwa.sort.bam") into bwa_bam

    script:
	"""
    echo "$read1+$read2" > ${id}_bwa.sort.bam
    """
    //    /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon bwa mem -M -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -t ${task.cpus} $genome_file $read1 $read2 | /data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon util sort -r $genome_file -o ${id}_bwa.sort.bam -t ${task.cpus} --sam2bam -i -

}

// mark duplicates, will switch to sentieon
process markdup {
    cpus 10
    input:
	set group, id, analysis_dir, file(sorted_bam) from bwa_bam

    output:
	set group, val(id), analysis_dir, file("${id}.markdup.bam") into bam_markdup
    
    script:
    """
    wc -l $sorted_bam > ${id}.markdup.bam
    """
    //    sambamba markdup --tmpdir /data/tmp -t ${task.cpus} $sorted_bam ${id}.markdup.bam

}

// split output channels, one for qc one for variantcalling
bam_markdup.into {
    bam_marked1
    bam_marked2
}

// post alignement qc
process post_align_qc {
    cpus 6

    input:
    set group, id, analysis_dir, file(markdup_bam) from bam_marked1
    output:
    set id, analysis_dir, file("${id}.bwa.qc") into bam_qc


    script:
    """
    echo "${analysis_dir}" > ${id}.bwa.qc
    """
  //  /data/bnf/scripts/postaln_qc.pl $markdup_bam $regions_bed ${id} ${task.cpus} $regions_bed $genome_file > ${id}.bwa.qc


}

// upload qc to CDM if upload flag is in use
process upload {
    publishDir "${OUTDIR}/postmap/exome", mode: 'copy', overwrite: 'true'
    input:
    set id, analysis_dir, file(qc) from bam_qc

    output:
    set id, file("${id}.qc.upload")
    when:
    params.upload =~ /true/

    script:
    """
    echo "${analysis_dir} ${qc}" > ${id}.qc.upload
    """
    
    // /data/bnf/scripts/register_sample.pl --run-folder ${analysis_dir} --sample-id ${id} --assay exome --qc ${qc}
}


process gvcf_combine {
    cpus 1

    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, id, analysis_dir, file(bam) from bam_marked2.groupTuple()
    
    output:
    file("${group}.concat.count") into results

    script:
    gvcfs = bam.collect{"$it"}.join(' --variant ')
    """
    echo "${gvcfs}" > ${group}.concat.count
    """



}

// skapa en pedfil, ändra input istället för sök ersätt?
process create_ped {
    input:
    set group, id, sex, mother, father, phenotype from ped

    output:
    file("${group}.ped") into ped_ch

    script:
    if ( sex =~ /F/) {
        sex = "2"
    }
    else {
        sex = "1"
    }
    if ( phenotype =~ /affected/ ) {
        phenotype = "2"
    }
    else {
        phenotype = "1"
    }
    if ( father == "" ) {
        father = "0"
    }
    if ( mother == "" ) {
        mother = "0"
    }
    """
    echo "${group}\t${id}\t${father}\t${mother}\t${sex}\t${phenotype}" > ${group}.ped
    """
}

ped_ch
    .collectFile(storeDir: "${OUTDIR}/ped/exome")
    .set{ lines }
    

// madeleine ped om familj
process madeleine_ped {
    publishDir "${OUTDIR}/ped/exome", mode: 'copy' , overwrite: 'true'
    input:
    file(ped) from lines

    output:
    file("hej")

    script:
    // funktion som räknar antalet individer i ped
    def bla = file("${OUTDIR}" + "/ped/exome/" + ped)
    no_lines =  bla.countLines()
    if (no_lines >= 2) {

    """
    echo "FAMILJ" > hej
    """
    }
    else {
    """
    echo "singel" > hej
    """
    }

}