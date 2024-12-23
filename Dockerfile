FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    gzip \
    bwa \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    fastp \
    bc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.9.tar.bz2 samtools-1.9

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
    && tar -xjf bcftools-1.9.tar.bz2 \
    && cd bcftools-1.9 \
    && ./configure \
    && make \
    && make install \
    && ln -s /usr/local/bin/bcftools /usr/bin/bcftools \
    && cd .. \
    && rm -rf bcftools-1.9.tar.bz2 bcftools-1.9

RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -xjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && ./configure \
    && make \
    && make install \
    && ln -s /usr/local/bin/bgzip /usr/bin/bgzip \
    && cd .. \
    && rm -rf htslib-1.9.tar.bz2 htslib-1.9

RUN mkdir -p /usr/local/files

COPY AMOSTRA_A.bam /usr/local/files/AMOSTRA_A.bam

COPY VariGene.sh /usr/local/VariGene.sh

COPY report.sh /usr/local/report.sh

#ADD INTERMEDIATE FILES
#COPY intermediate_files/* /usr/local/

RUN chmod +x /usr/local/VariGene.sh \
&& chmod +x /usr/local/report.sh

WORKDIR /usr/local

ENTRYPOINT ["/usr/local/VariGene.sh"]