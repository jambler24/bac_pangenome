FROM ubuntu:xenial

LABEL authors="Jon Ambler" \
      description="Docker image containing all requirements for bacterial genome assembly"


# for bandage :/ otherwise it complains about missing libGL.so.1
RUN apt-get update --fix-missing -qq

RUN apt-get upgrade -y

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    aufs-tools \
    automake \
    build-essential \
    bowtie2 \
    cmake \
    curl \
    g++ \
    openjdk-8-jre \
    openjdk-8-jdk \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    libx11-dev \
    pkg-config \
    zlib1g-dev \
    bzip2 \
    r-base \
    default-jre \
    git-core \
    bc \
    unzip \
    wget \
    xutils-dev \
    rsync \
    tar \
    software-properties-common \
    && apt-get clean \
    && apt-get purge

RUN apt-get install -y -q bedtools mcl cd-hit mafft prank fasttree parallel cpanminus
RUN cpanm -f Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict Devel::OverloadInfo Digest::MD5::File --installdeps
RUN cpanm -f List::Util --installdeps
RUN cpanm -f Bio::Roary --installdeps
RUN cpanm -f Bio::Roary::CommandLine::Roary --installdeps

RUN ls
ENV DST=/usr/bin
ENV URL=https://github.com/sanger-pathogens/Roary/tarball/master
RUN wget $URL -O $DST/sanger-pathogens-Roary.tar.gz && tar xvzf $DST/sanger-pathogens-Roary.tar.gz -C $DST

RUN apt-get install -y roary

RUN apt-get install -y -q libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl

# Roary works from here

RUN apt-get update -y

# Why did I have to install blast like this? Specific version requirements? Or specific path requirements?
RUN mkdir /tools/
RUN ls
RUN wget --no-verbose -O /tools/ncbi-blast-2.8.1+-x64-linux.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz
RUN tar -zxvf /tools/ncbi-blast-2.8.1+-x64-linux.tar.gz -C /tools/
ENV PATH=/tools/ncbi-blast-2.8.1+/bin:$PATH
RUN /tools/ncbi-blast-2.8.1+/bin/makeblastdb -h
RUN rm /usr/bin/makeblastdb
RUN which makeblastdb

RUN apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
RUN cpan Bio::Perl

RUN git clone https://github.com/rrwick/Unicycler.git $HOME/Unicycler && cd $HOME/Unicycler && python3 setup.py install && cd ~


RUN wget -O /tools/SPAdes-3.14.1-Linux.tar.gz http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz && tar -xzf /tools/SPAdes-3.14.1-Linux.tar.gz -C /tools/
ENV PATH=/tools/SPAdes-3.14.1-Linux/bin/:$PATH

RUN apt-get install -y fastqc

RUN wget -O /tools/Trimmomatic-0.38.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && cd /tools/ && unzip Trimmomatic-0.38.zip && cd ~
ENV PATH=/tools/Trimmomatic-0.38:$PATH

RUN curl -fksSL https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar > /usr/local/bin/picard.jar && \
    chmod +x /usr/local/bin/picard.jar

RUN R -e 'install.packages( c("reshape2","optparse", "BiocManager"), repos="http://cloud.r-project.org/");' && \
    apt-get update && apt-get install r-cran-ggplot2 -y -q
RUN R -e 'BiocManager::install("dupRadar");'

# install STAR
RUN curl -fksSL https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz | tar xz && \
    cp STAR-2.5.2b/bin/Linux_x86_64/* /usr/local/bin

# install bwe-mem
RUN git clone https://github.com/lh3/bwa.git && cd bwa; make; cp bwa /usr/local/bin

RUN apt-get install liblzma-dev -y

# SamTools
RUN wget -O /tools/htslib-1.9.tar.bz2 https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -vxjf /tools/htslib-1.9.tar.bz2 -C /tools/ && cd /tools/htslib-1.9/ && make

RUN wget -O /tools/samtools-1.10.tar.bz2 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && tar -vxjf /tools/samtools-1.10.tar.bz2 -C /tools/ && cd /tools/samtools-1.10 && make

RUN curl -O -L https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2 && tar xvfj samtools-0.1.18.tar.bz2 && cd samtools-0.1.18 && make

RUN wget -O /tools/bcftools-1.9.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar -vxjf /tools/bcftools-1.9.tar.bz2 -C /tools/ && cd /tools/bcftools-1.9 && make

ENV PATH="$PATH:/tools/bcftools-1.9"
ENV PATH="$PATH:/tools/samtools-1.10"
ENV PATH="$PATH:/tools/htslib-1.9"


RUN apt-get install -y python3-pip
RUN python3 -m pip install --upgrade cutadapt

RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o /tools/trim_galore.tar.gz && tar xvzf /tools/trim_galore.tar.gz -C /tools/
ENV PATH="$PATH:/tools/TrimGalore-0.6.5"

#gffread
RUN git clone https://github.com/gpertea/gffread /tools/gffread && cd /tools/gffread/ && make release && cd ~
ENV PATH="$PATH:/tools/gffread"


RUN python3 -m pip install quast

RUN apt-get install qtbase5-dev libqt5svg5-dev -y
ENV QT_SELECT=5
RUN git clone https://github.com/rrwick/Bandage.git /tools/Bandage && cd /tools/Bandage/ && qmake && make install


RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install python3.7-dev -y
RUN apt-get install python3.7 -y
RUN python3.7 -m pip install biopython
RUN git clone https://github.com/nigyta/dfast_core.git /tools/dfast
RUN ln -s /tools/dfast/dfast /usr/local/bin/
RUN ln -s /tools/dfast/scripts/dfast_file_downloader.py /usr/local/bin/

#RUN conda install -c bioconda multiqc

RUN apt-get install bowtie2 -y

RUN echo 'java -jar /usr/local/bin/picard.jar "$@"' > /usr/bin/picard && \
    chmod +x /usr/bin/picard

ENV PATH="/tools/samtools-1.10:$PATH"

# Install pilon

RUN mkdir /tools/pilon/
RUN curl -fkSL https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar > /usr/local/bin/pilon-1.23.jar && chmod +x /usr/local/bin/pilon-1.23.jar

RUN echo 'java -jar /usr/local/bin/pilon-1.23.jar "$@"' > /usr/bin/pilon && \
    chmod +x /usr/bin/pilon

RUN apt-get install python-setuptools pkg-config libfreetype6-dev libpng-dev python-matplotlib -y
RUN apt-get install python-pip -y
RUN pip install --upgrade pip==20.3.3
RUN pip install joblib
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz && tar -xzf quast-5.0.2.tar.gz && cd quast-5.0.2 && ./setup.py install
RUN cd ~

RUN git clone https://github.com/tseemann/prokka.git /tools/prokka
RUN /tools/prokka/bin/prokka --setupdb
ENV PATH=/tools/prokka/bin/:$PATH

RUN pip install multiqc
RUN apt-get install libgsl-dev -y

RUN pip install scoary
# One day...
#RUN python3 -m pip install scoary


# Preseq install
RUN Rscript -e 'install.packages("preseqR", repos="https://cloud.r-project.org")'
RUN curl -fsSL https://github.com/smithlabcode/preseq/releases/download/v2.0.3/preseq_v2.0.3.tar.bz2 -o /tools/preseq_v2.0.3.tar.bz2 && tar xvjf /tools/preseq_v2.0.3.tar.bz2 -C /tools/
RUN cd /tools/preseq/ && make && make install
ENV PATH=/tools/preseq/:$PATH

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 1

RUN update-alternatives --config python3

RUN python3 -m pip install --upgrade pip

# Install SPAdes (Must be version 3.13.0 for unicycler
RUN wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
RUN tar -xzf SPAdes-3.13.0-Linux.tar.gz
RUN cd SPAdes-3.13.0-Linux/bin/
ENV PATH=/SPAdes-3.13.0-Linux/bin/:$PATH

RUN python3 -m pip install --upgrade cutadapt

RUN git clone https://github.com/rrwick/Unicycler.git && cd Unicycler && python3 setup.py install && cd ~

RUN wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz
RUN gunzip linux.tbl2asn.gz
RUN mv linux.tbl2asn /tools/tbl2asn
RUN chmod a+x /tools/tbl2asn
ENV PATH=/tools/:$PATH

RUN spades.py --test

RUN unicycler --version
RUN cutadapt --version
RUN roary -w
RUN prokka -v
RUN tbl2asn -

#docker build -t bacterial_pangenome -f Dockerfile .

#docker run -v /var/run/docker.sock:/var/run/docker.sock -v /Volumes/External/bacterial_genome:/output --privileged -t --rm singularityware/docker2singularity:v2.6 -m "/shared_fs /custom_mountpoint2" bacterial_pangenome:latest

# docker run -it bacterial_pangenome:latest /bin/bash

