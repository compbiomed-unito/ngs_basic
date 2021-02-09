# BUILD STAGE: build from source
FROM python:3.8-slim-buster AS build
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev gcc g++ make autoconf git cmake wget

RUN mkdir -p /opt/bin

# bwa
RUN git clone --branch v0.7.17 https://github.com/lh3/bwa.git && cd bwa && make && mv bwa /opt/bin/bwa

# htslib
RUN git clone --branch 1.11 https://github.com/samtools/htslib.git && cd htslib && autoheader && autoconf && ./configure --prefix=/opt && make -j4 install

# samtools 
RUN git clone --branch 1.11 https://github.com/samtools/samtools.git && cd samtools && autoheader && autoconf && ./configure --prefix=/opt --with-htslib=/opt && make -j4 install

# bcftools
RUN git clone --branch 1.11 https://github.com/samtools/bcftools.git && cd bcftools && autoheader && autoconf && ./configure --prefix=/opt --with-htslib=/opt && make -j4 install

# bedtools
RUN git clone --branch v2.29.2 https://github.com/arq5x/bedtools2.git && cd bedtools2 && make -j4 && mv bin/* /opt/bin/

# bamtools (needed for covtobed)
RUN git clone --branch v2.5.1 https://github.com/pezmaster31/bamtools.git && mkdir bamtools/build && cd bamtools/build && cmake -DCMAKE_INSTALL_PREFIX=/opt .. && make -j4 install

# covtobed
RUN git clone --branch v1.2.0 https://github.com/telatin/covtobed.git && cd covtobed && c++ -std=c++11 -O3 -o covtobed *.cpp -lbamtools -I/opt/include/bamtools -L/opt/lib/ -lz && mv covtobed /opt/bin/

# bcl2fastq, from source, it is very long since it wants to build boost and cmake. The source archive from illumina is also very big (200Mb)
#RUN apt-get install -y cmake libboost-dev
#ADD bcl2fastq2-v2.20.0.422-Source.tar.gz /
#RUN mkdir /bcl2fastq/build && cd /bcl2fastq/build && export BOOST_INCLUDEDIR=/usr/include/boost/ && ../src/configure --prefix=/opt --with-cmake=/usr/bin/cmake && make install -j4

# bcl2fastq from rpm, seems the simplest way
COPY Binaries/bcl2fastq2-v2-20-0-linux-x86-64.zip /
RUN apt-get install -y --no-install-recommends rpm unzip && unzip bcl2fastq2-v2-20-0-linux-x86-64.zip && rpm -ivh --nodeps --prefix=/opt bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm

# fastqc from binary
COPY Binaries/fastqc_v0.11.9.zip /
RUN cd /opt && unzip /fastqc_v0.11.9.zip && chmod +x /opt/FastQC/fastqc && ln -s /opt/FastQC/fastqc /opt/bin/

# module psutil (needed by snakemake) must be compiled in the build stage (needs gcc) and packaged as binary wheel for the deployment stage
RUN pip install --upgrade --no-cache-dir pip
RUN mkdir -p /opt/wheels && cd /opt/wheels && pip wheel --no-cache-dir psutil

# DEPLOYMENT STAGE
FROM python:3.8-slim-buster

# fastqc launcher is in perl, so we need perl (25mb more) or to rewrite the launcher. Maybe installing only perl-base (which appears to be included in the base image) and the needed modules would be enough, but without the cpan module I do not know how to install other modules :(
RUN apt-get update && apt-get install -y --no-install-recommends perl

RUN apt-get install -y --no-install-recommends nano pigz

# add java runtime
#ADD Binaries/jre-8u271-linux-x64.tar.gz /opt/
RUN wget 'https://javadl.oracle.com/webapps/download/AutoDL?BundleId=244058_89d678f2be164786b292527658ca1605' -O jre.tar.gz && tar -xf jre.tar.gz -C /opt/

# copy programs compiled from the build stage
COPY --from=build /opt /opt

COPY Pipelines /opt/pipelines

# snakemake and cutadapt
RUN pip install --upgrade --no-cache-dir pip \
		&& pip install --no-cache-dir /opt/wheels/* snakemake cutadapt multiqc pandas openpyxl

# add paths to global profile
RUN echo 'export PATH="/opt/bin:/opt/jre1.8.0_271/bin:$PATH"; export LIBRARY_PATH="/opt/lib"; export LD_LIBRARY_PATH="/opt/lib"' >> /etc/profile

# check installed software, did not find a way to check bwa and java
RUN . /etc/profile && cutadapt --version && snakemake --version && bedtools --version && samtools --version && bcftools --version && fastqc --version && multiqc --version
