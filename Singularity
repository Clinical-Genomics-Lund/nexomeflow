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

%post
	#/opt/conda/bin/conda env create -f /exome_lower_python.yml
	apt -y install libz-dev
	apt -y install build-essential
	apt -y install gettext
	apt -y install cmake
	apt -y install libxml2-dev
	apt -y install libcurl4-openssl-dev
	apt -y install libssl-dev
	apt -y install r-base ant gradle openjdk-8-jdk
	git clone https://github.com/broadinstitute/picard.git
	cd picard
	./gradlew shadowJar
	cd ..
	/opt/conda/bin/conda env create -f /exome_conda.yml
	git clone https://github.com/piratical/Madeline_2.0_PDE.git
	cd Madeline_2.0_PDE
	./configure --with-include-gettext
	make
	make install
	/opt/conda/envs/exome_general/bin/pip install genmod
	#/opt/conda/bin/conda clean -a