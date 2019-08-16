/data/bnf/sw/nextflow/nextflow run exome.nf \
			    --csv /trannel/proj/cmd-pipe/test_real_trio.csv \
				--outdir /data/bnf/dev/viktor/NF_output \
				-with-singularity container_bwamem2.sif \
			    -with-dag test.dag.png \
				-with-timeline \
				-resume



