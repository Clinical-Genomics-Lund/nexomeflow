#!/usr/bin/env nextflow

params.fasta = "/data/bnf/ref/b37/human_g1k_v37_decoy.fasta"
params.bed = "/data/bnf/ref/agilent_clinical_research_exome_v2/S30409818_Regions.nochr.bed"
capture_kit = "Agilent_SureSelectCRE.V2"
OUTDIR = file(params.outdir)

// might move these to config
genome_file = file(params.fasta)
regions_bed = file(params.bed)
rank_model_s = "/rank_models/rank_model_cmd_v3_single_withoutmodels.ini"
rank_model = "/rank_models/rank_model_cmd_v3.ini"

REGCDM = "/data/bnf/scripts/register_sample.pl"
// VEP 
CADD = "/data/bnf/sw/.vep/PluginData/whole_genome_SNVs_1.4.tsv.gz"
VEP_FASTA = "/data/bnf/sw/.vep/homo_sapiens/87_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
MAXENTSCAN = "/data/bnf/sw/ensembl-vep-95/.vep/Plugins/MaxEntScan_scripts"
VEP_CACHE = "/data/bnf/sw/ensembl-vep-95/.vep"
GNOMAD = "/data/bnf/ref/b37/gnomad.exomes.r2.0.1.sites.vcf___.gz,gnomADg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH"
GERP = "/data/bnf/ref/annotations_dbs/VEP_conservation/All_hg19_RS.bw,GERP,bigwig"
PHYLOP =  "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phyloP100way.bw,phyloP100way,bigwig"
PHASTCONS = "/data/bnf/ref/annotations_dbs/VEP_conservation/hg19.100way.phastCons.bw,phastCons,bigwig"

SNPSIFT = "java -jar /opt/conda/envs/exome_general/share/snpsift-4.3.1t-1/SnpSift.jar"
CLINVAR = "/data/bnf/ref/annotations_dbs/clinvar_20190225.vcf.gz"
SWEGEN = "/data/bnf/ref/annotations_dbs/swegen_20170823/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz"
SPIDEX = "/data/bnf/ref/annotations_dbs/hg19_spidex.tsv.gz"




csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .set { fastq }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis) }
    .into { ped; all_ids; yml_diag }



// alignment using sentieon-bwa
process bwa_align {
    cpus 12
    input: 
    set clarity_sample_id, id, assay, sex, diagnosis, phenotype, group, father, mother, clarity_pool_id, platform, read1, read2, analysis_dir from fastq
    output:
    set group, id, analysis_dir, file("${id}_bwa.sort.bam"), file("${id}_bwa.sort.bam.bai") into bwa_bam
    script:
	"""
    sentieon bwa mem -M \\
    -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
    -t ${task.cpus} \\
    $genome_file \\
    ${read1} ${read2} | \\
    sentieon util sort \\
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
    sentieon driver -t ${task.cpus} -i $sorted_bam \\
    --algo LocusCollector --fun score_info SCORE.gz
    sentieon driver -t ${task.cpus} -i $sorted_bam \\
    --algo Dedup --rmdup --score_info SCORE.gz  \\
    --metrics DEDUP_METRIC_TXT ${id}.markdup.bam   
    """
}

// split output channels, 4 for qc 1 for variantcalling
bam_markdup.into {
    bam_marked1
    bam_marked2
    bam_marked3
    bam_marked4
    bam_marked5
    bam_marked6
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
    picard -Xmx12g CollectHsMetrics I=$bam O=${id}.markdup.bam.hsmetrics R=$genome_file BAIT_INTERVALS=$bait_i TARGET_INTERVALS=$bed_i
    """
  //  /data/bnf/scripts/postaln_qc.pl $markdup_bam $regions_bed ${id} ${task.cpus} $regions_bed $genome_file > ${id}.bwa.qc


}

process reads {
    cpus 1
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked3
    output:
    file("${id}.markdup.bam.reads") into qc_reads
    script:
    """
    sambamba flagstat -t ${task.cpus} $bam > ${id}.markdup.bam.reads
    """
}

process insertSize {
    errorStrategy 'ignore'
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked4
    output:
    file("${id}.markdup.bam.inssize") into qc_inssize
    script:
    """
    picard -Xmx12g CollectInsertSizeMetrics I=${bam} O=${id}.markdup.bam.inssize H=${id}.markdup.bam.ins.pdf STOP_AFTER=1000000
    """
}

process depthstats {
    cpus 10
    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked5
    output:
    set id, analysis_dir, file("${id}.basecov.bed") into qc_depth
    """
    sambamba depth base -c 0 -t ${task.cpus} -L $regions_bed $bam > ${id}.basecov.bed
    """
}
process combine_qc {
    input:
    set id, analysis_dir, file(depth) from qc_depth
    file(hs) from qc_hsmetrics
    file(reads) from qc_reads
    file(ins) from qc_inssize    
    output:
    set id, analysis_dir, file("${id}.bwa.qc") into qc_done
    """
    postaln_qc_nexomeflow.pl $hs $reads $ins $depth $id > ${id}.bwa.qc
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
    
    // $REGCDM --run-folder ${analysis_dir} --sample-id ${id} --assay exome --qc ${qc}
}

process DNAscope {
    cpus 12
    input:
    set group, id, analysis_dir, file(bam),file(bai) from bam_marked2
    output:
    set group, id, file("${id}.${group}.gvcf"), file("${id}.${group}.gvcf.idx") into gvcf
    script:
    type = mode == "family" ? "--emit_mode GVCF" : ""
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam \\
    --interval $regions_bed --algo DNAscope $type ${id}.${group}.gvcf
    """
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
    if (mode == "family" ) {
    ggvcfs = vcf.join(' -v ')
    """
    sentieon driver -t ${task.cpus} -r $genome_file --algo GVCFtyper \\
    -v $ggvcfs ${group}.combined.gvcf
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

// skapa en pedfil, ändra input istället för sök ersätt?
process create_ped {
    input:
    set group, id, sex, mother, father, phenotype, diagnosis from ped
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
    .collectFile(sort: true, storeDir: "${OUTDIR}/ped/exome")
    .into{ ped_mad; ped_peddy; ped_inher; ped_scout }
    


// madeline ped om familj
process madeline {
    //conda '/data/bnf/sw/miniconda3/envs/genmod'
    publishDir "${OUTDIR}/ped/exome", mode: 'copy' , overwrite: 'true'
    input:
    file(ped) from ped_mad
    output:
    set file("${ped}.madeline"), file("${ped}.madeline.xml") into madeline_ped
    when:
    mode == "family"
    script:
    """
    ped_parser -t ped $ped --to_madeline -o ${ped}.madeline
    madeline2 -L "IndividualId" ${ped}.madeline -o ${ped}.madeline -x xml
    """
}


// Splitting & normalizing variants:

process split_normalize {
    
    input:
    set group, file(vcf), file(idx) from g_gvcf
    output:
    set group, file("${group}.norm.DPAF.vcf") into split


    """
    vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
    bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
    /data/bnf/scripts/exome_DPAF_filter.pl ${group}.norm.vcf > ${group}.norm.DPAF.vcf
    """
}

// Annotating variants with VEP: 

process annotate_vep {
    container = 'container_VEP.sif'
    cpus 6
    input:
    set group, file(vcf) from split
    output:
    set group, file("${group}.vep.vcf") into vep
    """
    vep \\
    -i ${vcf} \\
    -o ${group}.vep.vcf \\
    --offline \\
    --merged \\
    --everything \\
    --vcf \\
    --no_stats \\
    --fork ${task.cpus} \\
    --force_overwrite \\
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
    -custom $PHASTCONS
    """
}

// Annotating variants with SnpSift 2

process snp_sift {
    //conda '/opt/conda/envs/exome_general'
    input:
    set group, file(vcf) from vep
    output:
    set group, file("${group}.clinvar.vcf") into snpsift
    """
    $SNPSIFT annotate $CLINVAR -info CLNSIG,CLNACC,CLNREVSTAT $vcf > ${group}.clinvar.vcf
    """

}

// Adding SweGen allele frequencies
process swegen_all {
    input:
    set group, file(vcf) from snpsift
    output:
    set group, file("${group}.swegen.vcf") into sweall
    """
    $SNPSIFT annotate $SWEGEN -name swegen -info AF $vcf > ${group}.swegen.vcf
    """
}
// Annotating variants with Genmod
process annotate_genmod {
    //conda '/data/bnf/sw/miniconda3/envs/genmod'
    input:
    set group, file(vcf) from sweall
    output:
    set group, file("${group}.genmod.vcf") into genmod
    """
    genmod annotate --spidex $SPIDEX --cadd_file $CADD --annotate_regions $vcf -o ${group}.genmod.vcf
    """
}

// # Annotating variant inheritance models:
process inher_models {
    //conda '/data/bnf/sw/miniconda3/envs/genmod'
    input:
    set group, file(vcf) from genmod
    file(ped) from ped_inher
    output:
    set group, file("${group}.models.vcf") into inhermod
    """
    genmod models $vcf -f $ped > ${group}.models.vcf
    """
}


// Extracting most severe consequence: 
// Modifying annotations by VEP-plugins, and adding to info-field: 
// Modifying CLNSIG field to allow it to be used by genmod score properly:
process modify_vcf {
    input:
    set group, file(vcf) from inhermod
    output:
    set group, file("${group}.mod.vcf") into mod_vcf
    """
    /opt/bin/modify_vcf_nexomeflow.pl $vcf > ${group}.mod.vcf
    """
} 


// Adding loqusdb allele frequency to info-field: 
// ssh needs to work from anywhere, filesystems mounted on cmdscout
process loqdb {
    //container = 'container_mongodb.sif'
    input:
    set group, file(vcf) from mod_vcf
    output:
    set group, file("${group}.loqdb.vcf") into loqdb_vcf
    """
    /opt/bin/loqus_db_filter.pl $vcf > ${group}.loqdb.vcf
    """
    // ssh cmdscout1.lund.skane.se 
}
// Marking splice INDELs: 
// Annotating delins with cadd: 
process splicecadd {
    input:
    set group, file(vcf) from loqdb_vcf
    output:
    set group, file("${group}.marksplice.cadd.vcf") into splice_cadd
    """
    /opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
    add_missing_CADDs_1.4.sh -i ${group}.marksplice.vcf -o ${group}.marksplice.cadd.vcf -t ${OUTDIR}/tmp/exome/${group}.addcadd
    """
}
// Scoring variants: 
// Adjusting compound scores: 
// Sorting VCF according to score: 

process genmodscore {
    //conda '/data/bnf/sw/miniconda3/envs/genmod'
    input:
    set group, file(vcf) from splice_cadd
    output:
    set group, file("${group}.scored.vcf") into scored_vcf
    script:
    if (mode == "family") {
        """
        genmod score -i $group -c $rank_model -r $vcf -o ${group}.score1.vcf
        genmod compound ${group}.score1.vcf > ${group}.score2.vcf
        genmod sort -p -f $group ${group}.score2.vcf -o ${group}.scored.vcf
        """
    }
    else {
        """
        genmod score -i $group -c $rank_model_s -r $vcf -o ${group}.score1.vcf
        genmod sort -p -f $group ${group}.score1.vcf -o ${group}.scored.vcf
        """
    }

}

// Bgzipping and indexing VCF: 
process vcf_completion {
    cpus 5
    publishDir "${OUTDIR}/vcf/exome/", mode: 'copy', overwrite: 'true'
    input:
    set group, file(vcf) from scored_vcf
    output:
    set group, file("${group}.scored.vcf.gz"), file("${group}.scored.vcf.gz.tbi") into vcf_done
    """
    bgzip -@ ${task.cpus} $vcf -f
    tabix ${vcf}.gz -f
    """
}

vcf_done.into {
    vcf_done1
    vcf_done2
    vcf_done3
}

// Running PEDDY: 
process peddy {
    //conda '/opt/conda/envs/peddy'
    publishDir "${OUTDIR}/ped/exome", mode: 'copy' , overwrite: 'true'
    cpus 6
    input:
    file(ped) from ped_peddy
    set group, file(vcf), file(idx) from vcf_done1
    output:
    set file("${group}.ped_check.csv"),file("${group}.background_pca.json"),file("${group}.peddy.ped"),file("${group}.html"), file("${group}.het_check.csv"), file("${group}.sex_check.csv"), file("${group}.vs.html") into peddy_files
    """
    source activate peddy
    python -m peddy -p ${task.cpus} $vcf $ped --prefix $group
    """
}


// Running gSNP:
process gnsp {
    container = 'container_mongodb.sif'
    publishDir "${OUTDIR}/tmp/exome/gSNP", mode: 'copy' , overwrite: 'true'
    input:
    set group, file(vcf), file(idx) from vcf_done3
    set group, id, sex, mother, father, phenotype, diagnosis from all_ids.groupTuple()
    output:
    set file("${group}_gSNP.tsv"), file("${group}_gSNP.png")
    when:
    mode == "family"
    script:
    ids = id.join(' ')
    """
    source activate exome_general
    perl /opt/bin/gSNP.pl $vcf ${group}_gSNP $ids
    """
}

// Uploading case to scout:
process create_yaml {

    input:
    set group, id, analysis_dir, file(bam), file(bai) from bam_marked6.groupTuple(sort: true)
    set group, file(vcf), file(idx) from vcf_done2
    set file(ped_check),file(json),file(peddy_ped),file(html), file(hetcheck_csv), file(sexcheck), file(vs_html) from peddy_files
    file(ped) from ped_scout
    set file(madeline), file(xml) from madeline_ped.ifEmpty()
    set group, id, sex, mother, father, phenotype, diagnosis from yml_diag
    output:
    set group, file("${group}.yml") into yaml
    script:
    bams = bam.join(',')
    madde = xml.name != '' ? "$xml" : "single"
    """
    /trannel/proj/cmd-pipe/nexomeflow/bin/create_yml.pl $bams $ped $group $vcf $madde $peddy_ped $ped_check $sexcheck $OUTDIR $diagnosis > ${group}.yml
    """

}
// # Registering contamination data:

// # Uploading variants to loqusdb:


//email results
