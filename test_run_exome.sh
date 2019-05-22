/data/bnf/sw/nextflow/nextflow run exome.nf \
			        --csv /trannel/proj/cmd-pipe/test_real_trio.csv \
				--outdir /trannel/proj/cmd-pipe/test-files/run_singularity \
				-with-singularity container_2019-05-14.sif \
			      	-with-dag test.dag.png \
				-with-timeline \
				-resume



