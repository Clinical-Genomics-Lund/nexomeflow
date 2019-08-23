#!/bin/bash
# A script to run cadd v.1.3 on delins found in a genmod annotate processed file.
# Created by PP 20180522.
#
# 20180823, path for score.sh changed from /data/cadd/1.3/bin/score.sh to /trannel/cadd/1.3/bin/score.sh
#
# check for arguments
if [[ $# != 6 ]] ; then
    echo 'How to use this script:'
    echo 'add_missing_CADDs.sh -i <genmod input> -o <output> -t <tmpdir>'
    exit 0
fi
#
# set python enviroment
#source /data/bnf/sw/miniconda3/envs/cadd/bin/activate cadd
source activate cadd-env
#
# create outdir, cp file and unzip file
tmpdir=$6 # dir temporary files
export tmpdir # export your shell variable so it is visible to subprocesses
>&2 echo "create tmpdir, cp file..."
if [ ! -d "${tmpdir}" ]; then
    mkdir ${tmpdir}
fi
cp $2 ${tmpdir}
gm_vcf=$(basename $2) # genmod file (*.genmod)
export gm_vcf # export your shell variable so it is visible to subprocesses
gm_out=$4 # genmod file out
>&2 echo "Done!"
# pull out indels from the file
>&2 echo "pull out all indels from vcf and collect them in separate file..."
perl -naF'\t' -e ' next if(/^#/); if( (length($F[3]) > 1) || (length($F[4]) > 1) ){ print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\n"; } ' ${tmpdir}/${gm_vcf} | gzip -c > ${tmpdir}/${gm_vcf}.indels.vcf.gz
>&2 echo "Done!"
# run score.sh (cadd scores without anno)
>&2 echo "start score.sh on all ins/dels from ${gm_vcf}.indels..."
CADD_2.sh -g GRCh37 -o ${tmpdir}/${gm_vcf}.indel_CADDs.gz ${tmpdir}/${gm_vcf}.indels.vcf.gz
>&2 echo "Done!"
# add CADD-scores to original file (*.genmod)
gunzip -f -k ${tmpdir}/${gm_vcf}.indel_CADDs.gz
>&2 echo "Complement original vcf with indels CADD from ${gm_vcf}.indel_CADDs..."
perl -naF'\t' -e '$i++; if($i == 1){ @cadd_file = `cat $ENV{tmpdir}/$ENV{gm_vcf}.indel_CADDs`; for $row (@cadd_file){ chomp($row); @col = split(/\t/, $row); $cadd{"$col[0]\t$col[1]\t$col[2]\t$col[3]"} = $col[5]; } } if( ($cadd{"$F[0]\t$F[1]\t$F[3]\t$F[4]"}) && (!/CADD=/) ){ s/\tGT:AD:/;CADD=$cadd{"$F[0]\t$F[1]\t$F[3]\t$F[4]"}\tGT:AD:/; print; }else{ print; }' ${tmpdir}/${gm_vcf} > ${gm_out}
>&2 echo "Done! Original vcf with inserted CADD scores for indels can be found in ${gm_out}"
