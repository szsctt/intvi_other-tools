#https://hub.docker.com/_/python

#  $ docker build . -t szsctt/polyidus:latest -t szsctt/polyidus:3
#  $ docker run --rm -it szsctt/polyidus:latest /bin/bash
#  $ docker push szsctt/polyidus:latest
#  $ docker push szsctt/polyidus:3

# cd ~/Documents/Projects/viInt/experiments/expt6_simulations/polyidus
# docker run -v$(pwd):/home --rm szsctt/polyidus:3 python /usr/src/app/src/polyidus.py /home/data/hg38/hg38_bwt2_index /home/data/hpv16/hpv16_bowt_ind --fastq /home/data/fastqfiles/SiHa_R1.fastq.gz /home/data/fastqfiles/SiHa_R2.fastq.gz --outdir /home/polyidusOutput

FROM continuumio/miniconda3:4.8.2

WORKDIR /usr/src/app

COPY polyidus.yml .

RUN /opt/conda/bin/conda install \
		-c bioconda -c conda-forge \
		bowtie2=2.2 \
		bedtools=2.26 \
		samtools=1.9 \
		pysam=0.15 \
		python=3 \
		psutil=5.4 \
		pandas=0.22 &&\
	/opt/conda/bin/conda clean --all -y 

COPY polyidus .

CMD [ "python3", "src/polyidus.py" ]


