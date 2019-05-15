Bootstrap:docker
From:nfcore/base

%labels
	MAINTAINER Viktor henmyr <viktor.henmyr@skane.se>
	DESCRIPTION Singularity container for CMD exome pipeline
	VERSION 0.0.1

%environment
	PATH=/opt/conda/envs/exome_general/bin:/opt/sentieon-genomics-201808.01/bin/:$PATH
	PICARD_HOME=/opt/conda/envs/exome_conda/share/picard-2.18.26-0/


%files
    exome_conda.yml /
	exome_lower_python.yml /
    bin/ /opt
	/data/bnf/sw/sentieon/sentieon-genomics-201808.01 /opt
	MongoDB.pm /usr/share/perl5

%post
	/opt/conda/bin/conda env create -f /exome_conda.yml
	#/opt/conda/envs/exome_general/bin/cpan MongoDB
	/opt/conda/bin/conda env create -f /exome_lower_python.yml
	apt -y install libz-dev build-essential gettext \
	cmake libxml2-dev libcurl4-openssl-dev libssl-dev ant \
	#gradle openjdk-8-jdk
	#git clone https://github.com/broadinstitute/picard.git
	#cd picard
	#./gradlew shadowJar
	#cd ..
	git clone https://github.com/piratical/Madeline_2.0_PDE.git
	cd Madeline_2.0_PDE
	./configure --with-include-gettext
	make
	make install
	/opt/conda/envs/exome_general/bin/pip install genmod
	/opt/conda/bin/conda clean -a