Bootstrap:docker
From:nfcore/base

%labels
	MAINTAINER Viktor henmyr <viktor.henmyr@skane.se>
	DESCRIPTION Singularity container for CMD exome pipeline
	VERSION 0.0.1

%environment
	PATH=/opt/conda/envs/exome_general/bin:/opt/sentieon-genomics-201808.01/bin/:$PATH
	PICARD_HOME=/opt/conda/envs/CMD-twist/share/picard-2.18.26-0/


%files
        exome_conda.yml /
        bin/ /opt
	/data/bnf/sw/sentieon/sentieon-genomics-201808.01 /opt

%post
	/opt/conda/bin/conda env create -f /exome_conda.yml
	apt -y install libz-dev
	apt -y install build-essential
	/opt/conda/envs/exome_general/bin/pip install genmod
	/opt/conda/envs/exome_general/bin/pip install ped_parser
	#/opt/conda/bin/conda clean -a