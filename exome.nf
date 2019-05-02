#!/usr/bin/env nextflow

params.fasta = "/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
params.bed = "/data/bnf/ref/agilent_clinical_research_exome_v2/S30409818_Regions.nochr.bed"

OUTDIR = file(params.outdir)

// might move these to config
genome_file = file(params.fasta)
regions_bed = file(params.bed)




// SOFTWARE-BIN
VCFMULTI = "/data/bnf/sw/vcflib/bin/vcfbreakmulti"
BCFTOOLS = "/data/bnf/sw/bcftools-1.3.1/bcftools"

// VEP 
VEP = "/data/bnf/sw/ensembl-vep-95/vep"
    FIX_VEP = "/data/bnf/scripts/fix_empty_vep_vcf.pl"
    CADD = "/data/bnf/sw/.vep/PluginData/whole_genome_SNVs_1.4.tsv.gz"
    VEP_FASTA = "/data/bnf/sw/.vep/homo_sapiens/87_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
    MAXENTSCAN = "/data/bnf/sw/ensembl-vep-95/.vep/Plugins/MaxEntScan_scripts"
    VEP_CACHE = "/data/bnf/sw/ensembl-vep-95/.vep"
    GNOMAD = "/data/bnf/ref/b37/gnomad.exomes.r2.0.1.sites.vcf___.gz,gnomADg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH"
    GERP = "/data/bnf/ref/annotations_dbs/VEP_conservation/All_hg19_RS.bw,GERP,bigwig"
    PHYLOP =  "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phyloP100way.bw,phyloP100way,bigwig"
    PHASTCONS = "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phastCons.bw,phastCons,bigwig"
SNPSIFT = "java -jar /data/bnf/sw/snpEff/4.3/SnpSift.jar"
    CLINVAR = "/data/bnf/ref/annotations_dbs/clinvar_20190225.vcf.gz"










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

process DNAscope {
    cpus 6
    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, id, analysis_dir, file(bam) from bam_marked2
    output:
    set group, id, file("${id}.${group}.gvcf") into gvcf
    script:
    """
    echo "${id} hej hej hej" > ${id}.${group}.gvcf
    """
// --emit_mode GVCF
}



process gvcf_combine {
    cpus 6

    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, id, file(bam) from gvcf.groupTuple()

    output:
    set group, file("${group}.concat.count") into g_gvcf

    script:
    // Om fler än en vcf, GVCF combine annars döp om och skickade vidare
    // lägg till lista för alla outcomes, ensam vcf/multiple vcf
    bam = bam + [1]
    // antal element
    bam_l =  bam.size()
    // om 3 eller fler == multivcf
    if (bam_l >= 3) {
        bam = bam - [1]
        ggvcfs = bam.join(' --variant ')
    """
    echo "${ggvcfs}" > ${group}.concat.count
    """
    }
    // annars ensam vcf, skicka vidare
    else {
    bam = bam - [1]
    ggvcf = bam.join('')
    """
    mv ${ggvcf} ${group}.concat.count
    """
    }

}

process gvcf_genotype {
    cpus 6
    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'

    input:
    set group, file(vcf) from g_gvcf
    output:
    set group, file("${group}.vcf") into split
    """
    echo "FINAL" > ${group}.vcf
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
// collects each individual's ped-line and creates on ped-file
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


// # Splitting & normalizing variants: 6427-13

process split_normalize {
    
    input:
    set group, file(vcf) from split
    output:
    set group, file("${group}.norm.DPAF.vcf") into vep


    """
    echo "$VCFMULTI ${vcf}" > ${group}.multibreak.vcf
    echo "$BCFTOOLS norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf" > ${group}.norm.vcf
    echo "/data/bnf/scripts/exome_DPAF_filter.pl ${group}.norm.vcf" > ${group}.norm.DPAF.vcf
    """
}

// # Annotating variants with VEP: 6427-13

process annotate_vep {
    publishDir "${OUTDIR}/tmp/exome", mode: 'copy' , overwrite: 'true'
    input:
    set group, file(vcf) from vep

    output:
    set group, file("${group}.vep.vcf") into snpsift
    
    """
    echo \\
    "$VEP \\
    -i ${vcf} \\
    -o ${group}.vep.vcf \\
    --offline \\
    --merged \\
    --plugin CADD $CADD \\
    --plugin LoFtool \\
    --plugin MaxEntScan,$MAXENTSCAN,SWA,NCSS \\
    --fasta $VEP_FASTA \\
    --dir_cache $VEP_CACHE \\
    --dir_plugins $VEP_CACHE/Plugins \\
    --distance 200 \\
    -cache \\
    -custom $GNOMAD \\
    -custom $GERP \\
    -custom $PHYLOP \\
    -custom $PHASTCONS \\
    " > ${group}.vep.vcf
    echo "$FIX_VEP $vcf ${group}.vep.vcf" > ${group}.vep.vcf
    """

}
// # Annotating variants with SnpSift 2: 6427-13

process snp_sift {

    input:
    set group, file(vcf) from snpsift
    output:
    set group, file("${group}.clinvar.vep") into clinmod
    """
    echo "$SNPSIFT annotate $CLINVAR -info CLNSIG,CLNACC,CLNREVSTAT $vcf" > ${group}.clinvar.vep
    """

}

// # Modifying CLNSIG field to allow it to be used by genmod score properly: 6427-13

// # Adding SweGen allele frequencies: 6427-13

// # Annotating variants with Genmod: 6427-13

// # Marking splice INDELs: 6427-13

// # Annotating delins with cadd: 6427-13

// # Annotating variant inheritance models: 6427-13

// # Extracting most severe consequence: 6427-13

// # Modifying annotations by VEP-plugins, and adding to info-field: 6427-13

// # Adding loqusdb allele frequency to info-field: 6427-13

// # Scoring variants: 6427-13

// # Adjusting compound scores: 6427-13

// # Sorting VCF according to score: 6427-13

// # Bgzipping and indexing VCF: 6427-13

// # Renaming VCF: 6427-13

// # Indexing VCF: 6427-13

// # Creating madeline pedigree: 6427-13

// # Running PEDDY: 6427-13

// # Running gSNP: 6427-13

// # Uploading case to scout: 6427-13

// # Registering contamination data: 6427-13

// # Uploading variants to loqusdb: 6427-13

