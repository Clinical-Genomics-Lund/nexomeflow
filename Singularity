Bootstrap:docker
From:nfcore/base

%labels
	MAINTAINER Viktor henmyr <viktor.henmyr@skane.se>
	DESCRIPTION Singularity container for CMD exome pipeline
	VERSION 0.0.1

%environment
	PATH=/opt/conda/envs/exome_general/bin:/opt/sentieon-genomics-201808.05/bin/:/opt/bin:/opt/conda/envs/peddy/bin:$PATH
	PICARD_HOME=/opt/conda/envs/exome_conda/share/picard-2.18.26-0/
	export PERL5LIB=$PERL5LIB:/opt/conda/envs/exome_general/lib/site_perl/5.26.2/
	export PERL5LIB=$PERL5LIB:/opt/conda/envs/exome_general/lib/site_perl/5.26.2/x86_64-linux-thread-multi/
	export PERL5LIB=$PERL5LIB:/opt/bin/

%files
    exome_conda.yml /
	exome_lower_python.yml /
	cadd_environment.yml /
    bin/ /opt
	/data/bnf/sw/sentieon/sentieon-genomics-201808.05 /opt
	rank_models/ /
%post
	/opt/conda/bin/conda env create -f /cadd_environment.yml
	/opt/conda/bin/conda env create -f /exome_conda.yml
	apt -y install libz-dev build-essential gettext cmake libxml2-dev libcurl4-openssl-dev libssl-dev make
	/opt/conda/envs/exome_general/bin/cpanm Path::Tiny --force
	/opt/conda/envs/exome_general/bin/cpanm MongoDB::Collection
	git clone https://github.com/piratical/Madeline_2.0_PDE.git
	cd Madeline_2.0_PDE
	./configure --with-include-gettext
	make
	make install
	/opt/conda/envs/exome_general/bin/pip install genmod
	/opt/conda/bin/conda env create -f /exome_lower_python.yml
	/opt/conda/bin/conda clean -a