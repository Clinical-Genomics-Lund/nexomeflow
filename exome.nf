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
SENT = "/data/bnf/sw/sentieon/sentieon-genomics-201808.01/bin/sentieon"
GATK = "/data/bnf/sw/gatk/4.1.0.0/gatk --java-options '-Xmx12G'"
//PICARD = "java -jar /data/bnf/sw/picard/2.8.1/picard.jar"
PICARD = "java -Xmx12g -jar /data/bnf/sw/picard/2.18.15/picard/build/libs/picard.jar"
SAMBAMBA = "/data/bnf/sw/sambamba/0.6.7/sambamba"
// PERL
POSTQC = "/trannel/proj/cmd-pipe/nexomeflow/bin/postaln_qc_nexomeflow.pl"
// VEP 
VEP = "/data/bnf/sw/ensembl-vep-95/vep"
    CADD = "/data/bnf/sw/.vep/PluginData/whole_genome_SNVs_1.4.tsv.gz"
    VEP_FASTA = "/data/bnf/sw/.vep/homo_sapiens/87_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
    MAXENTSCAN = "/data/bnf/sw/ensembl-vep-95/.vep/Plugins/MaxEntScan_scripts"
    VEP_CACHE = "/data/bnf/sw/ensembl-vep-95/.vep"
    GNOMAD = "/data/bnf/ref/b37/gnomad.exomes.r2.0.1.sites.vcf___.gz,gnomADg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH"
    GERP = "/data/bnf/ref/annotations_dbs/VEP_conservation/All_hg19_RS.bw,GERP,bigwig"
    PHYLOP =  "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phyloP100way.bw,phyloP100way,bigwig"
    PHASTCONS = "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phastCons.bw,phastCons,bigwig"
// SNPSIFT
SNPSIFT = "java -jar /data/bnf/sw/snpEff/4.3/SnpSift.jar"
    CLINVAR = "/data/bnf/ref/annotations_dbs/clinvar_20190225.vcf.gz"
    CLINMOD = "/data/bnf/scripts/modify_CLNSIG_NEW.pl"
    SWEGEN = "/data/bnf/ref/annotations_dbs/swegen_20170823/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz"
// GENMOD
GENMOD = "/home/bjorn/miniconda2/bin/genmod"
    SPIDEX = "/data/bnf/ref/annotations_dbs/hg19_spidex.tsv.gz"
    MARKSPLICE = "/data/bnf/scripts//mark_spliceindels.pl"
    ADDCADD = "/data/bnf/scripts/add_missing_CADDs_1.4.sh"
    MOSTSEV = "/data/bnf/scripts/most_severe_consequence.pl"



csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
println(mode)

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
    cpus 8
    publishDir "$OUTDIR/tmp/exome/", mode: 'copy', overwrite: 'true'

    input: 
    set clarity_sample_id, id, assay, sex, diagnosis, phenotype, group, father, mother, clarity_pool_id, platform, read1, read2, analysis_dir from fastq

    output:
    set group, id, analysis_dir, file("${id}_bwa.sort.bam"), file("${id}_bwa.sort.bam.bai") into bwa_bam

    script:
	"""
    $SENT bwa mem -M \\
    -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
    -t ${task.cpus} \\
    $genome_file \\
    ${read1} ${read2} | \\
    $SENT util sort \\
    -r $genome_file \\
    -o ${id}_bwa.sort.bam \\
    -t ${task.cpus} \\
    --sam2bam -i -
    """
    

}

// mark duplicates
process markdup {
    cpus 10
    input:
	set group, id, analysis_dir, file(sorted_bam), file(bai) from bwa_bam

    output:
	set group, val(id), analysis_dir, file("${id}.markdup.bam"), file("${id}.markdup.bam.bai") into bam_markdup
    
    script:
    """
    $SENT driver -t ${task.cpus} -i $sorted_bam \\
    --algo LocusCollector --fun score_info SCORE.gz
    $SENT driver -t ${task.cpus} -i $sorted_bam \\
    --algo Dedup --rmdup --score_info SCORE.gz  \\
    --metrics DEDUP_METRIC_TXT ${id}.markdup.bam   
    """
     //sambamba markdup --tmpdir /data/tmp -t ${task.cpus} $sorted_bam ${id}.markdup.bam
}

// split output channels, 4 for qc 1 for variantcalling
bam_markdup.into {
    bam_marked1
    bam_marked2
    bam_marked3
    bam_marked4
    bam_marked5
}

// post alignement qc BLOCK
process hsmetrics {
    cpus 6
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked1
    output:
    file("${id}.markdup.bam.hsmetrics") into qc_hsmetrics
    script:
    bait_i = regions_bed + ".interval_list"
    bed_i = regions_bed + ".interval_list"
    """
    $PICARD CollectHsMetrics I=$bam O=${id}.markdup.bam.hsmetrics R=$genome_file BAIT_INTERVALS=$bait_i TARGET_INTERVALS=$bed_i
    """
  //  /data/bnf/scripts/postaln_qc.pl $markdup_bam $regions_bed ${id} ${task.cpus} $regions_bed $genome_file > ${id}.bwa.qc


}

process reads {
    cpus 6
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked3
    output:
    file("${id}.markdup.bam.reads") into qc_reads
    script:
    """
    $SAMBAMBA flagstat -t ${task.cpus} $bam > ${id}.markdup.bam.reads
    """
}

process insertSize {
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked4
    output:
    file("${id}.markdup.bam.inssize") into qc_inssize
    script:
    """
    $PICARD CollectInsertSizeMetrics I=${bam} O=${id}.markdup.bam.inssize H=${id}.markdup.bam.ins.pdf STOP_AFTER=1000000
    """
}

process depthstats {
    cpus 6
    publishDir "${OUTDIR}/postmap/exome", mode: 'copy', overwrite: 'true'
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked5
    file(hs) from qc_hsmetrics
    file(reads) from qc_reads
    file(ins) from qc_inssize    
    output:
    set id, analysis_dir, file("${id}.bwa.qc") into qc_done
    """
    $SAMBAMBA depth base -c 0 -t ${task.cpus} -L $regions_bed $bam > ${id}.basecov.bed
    $POSTQC $hs $reads $ins ${id}.basecov.bed > ${id}.bwa.qc
    """
}


// upload qc to CDM if upload flag is in use
process upload {
    publishDir "${OUTDIR}/postmap/exome", mode: 'copy', overwrite: 'true'
    input:
    set id, analysis_dir, file(qc) from qc_done

    output:
    set id, file("${id}.qc.upload")
    when:
    params.upload =~ /true/

    script:
    """
    cat $qc > ${id}.qc.upload
    """
    
    // /data/bnf/scripts/register_sample.pl --run-folder ${analysis_dir} --sample-id ${id} --assay exome --qc ${qc}
}

process DNAscope {
    cpus 6
   // publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, id, analysis_dir, file(bam),file(bai) from bam_marked2
    output:
    set group, id, file("${id}.${group}.gvcf"), file("${id}.${group}.gvcf.idx") into gvcf
    script:
    
    """
    $SENT driver -t ${task.cpus} -r $genome_file -i $bam \\
    --algo DNAscope --emit_mode GVCF ${id}.${group}.gvcf
    """
    
// --emit_mode GVCF
//echo "${id} hej hej hej" > ${id}.${group}.gvcf
}



process gvcf_combine {
    cpus 6

    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, id, file(vcf), file(idx) from gvcf.groupTuple()

    output:
    set group, file("${group}.combined.gvcf"), file("${group}.combined.gvcf.idx") into g_gvcf

    script:
    // Om fler än en vcf, GVCF combine annars döp om och skickade vidare
    if (mode =~ /family/ ) {
    ggvcfs = vcf.join(' --variant ')
    """
    echo "${ggvcfs}" > ${group}.concat.count
    """
    }
    // annars ensam vcf, skicka vidare
    else {
    ggvcf = vcf.join('')
    gidx = idx.join('')
    """
    mv ${ggvcf} ${group}.combined.gvcf
    mv ${gidx} ${group}.combined.gvcf.idx
    """
    }

}

process gvcf_genotype {
    cpus 6
    publishDir "${OUTDIR}/tmp/exome/", mode: 'copy', overwrite: 'true'

    input:
    set group, file(vcf), file(idx) from g_gvcf
    output:
    set group, file("${group}.gvcf"), file("${group}.gvcf.idx") into split
    """
    $GATK GenotypeGVCFs \\
    -R $genome_file \\
    -L $regions_bed \\
    --variant $vcf \\
    -O ${group}.gvcf
    """
//echo "FINAL" > ${group}.vcf
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
    .into{ ped_mad; peddy; ped_inher }
    

// madeleine ped om familj
process madeleine {
    publishDir "${OUTDIR}/ped/exome", mode: 'copy' , overwrite: 'true'
    input:
    file(ped) from ped_mad

    output:
    file("hej")

    script:
    if (mode =~ /family/ ) {

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


// Splitting & normalizing variants:

process split_normalize {
    
    input:
    set group, file(vcf), file(idx) from split
    output:
    set group, file("${group}.norm.DPAF.vcf") into vep


    """
    $VCFMULTI ${vcf} > ${group}.multibreak.vcf
    $BCFTOOLS norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
    /data/bnf/scripts/exome_DPAF_filter.pl ${group}.norm.vcf > ${group}.norm.DPAF.vcf
    """
}

// Annotating variants with VEP: 

process annotate_vep {
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
    -custom $PHASTCONS" > ${group}.vep.vcf
    """
}

// Annotating variants with SnpSift 2

process snp_sift {

    input:
    set group, file(vcf) from snpsift
    output:
    set group, file("${group}.clinvar.vcf") into clinmod
    """
    echo "$SNPSIFT annotate $CLINVAR -info CLNSIG,CLNACC,CLNREVSTAT $vcf" > ${group}.clinvar.vcf
    """

}

// # Modifying CLNSIG field to allow it to be used by genmod score properly:
process clinsigmod {
    input:
    set group, file(vcf) from clinmod
    output:
    set group, file("${group}.clinmod.vcf") into sweall
    """
    echo "$CLINMOD $vcf" > ${group}.clinmod.vcf
    """
}
// Adding SweGen allele frequencies
process swegen_all {
    input:
    set group, file(vcf) from sweall
    output:
    set group, file("${group}.swegen.vcf") into splice
    """
    echo "$SNPSIFT annotate $SWEGEN -name swegen -info AF $vcf" > ${group}.swegen.vcf
    """
}
// Annotating variants with Genmod
// Marking splice INDELs: 
// Annotating delins with cadd: 
process annotate_genmod {
    input:
    set group, file(vcf) from splice
    output:
    set group, file("${group}.genmod.marksplice.addcadd.vcf") into inheritance
    """
    echo "$GENMOD annotate --spidex $SPIDEX --annotate_regions $vcf -o ${group}.genmod.vcf" 
    echo "$MARKSPLICE ${group}.genmod.vcf" > ${group}.genmod.marksplice.vcf
    echo "$ADDCADD -i ${group}.genmod.marksplice.vcf -o ${group}.genmod.marksplice.addcadd.vcf -t ${OUTDIR}/tmp/exome/${group}.addcadd" > ${group}.genmod.marksplice.addcadd.vcf
    """
}

// # Annotating variant inheritance models:
// # Extracting most severe consequence: 
process inher_models {
    //publishDir "${OUTDIR}/tmp/exome", mode: 'copy' , overwrite: 'true'

    input:
    set group, file(vcf) from inheritance
    file(ped) from ped_inher
    output:
    set group, file("${group}.models.mostsev.vcf") into vepmod
    """
    echo "$GENMOD models $vcf -f $ped" > ${group}.models.vcf
    echo "$MOSTSEV < ${group}.models.vcf > ${group}.models.mostsev.vcf" > ${group}.models.mostsev.vcf
    """
}


// # Modifying annotations by VEP-plugins, and adding to info-field: 

// # Adding loqusdb allele frequency to info-field: 

// # Scoring variants: 

// # Adjusting compound scores: 

// # Sorting VCF according to score: 

// # Bgzipping and indexing VCF: 

// # Running PEDDY: 

// # Running gSNP:

// # Uploading case to scout:
// process scout {

//     input:


// }
// # Registering contamination data:

// # Uploading variants to loqusdb:

