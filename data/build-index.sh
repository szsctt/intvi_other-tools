#!/bin/bash
set -e

module load  bowtie/2.3.4
module load blast+/2.2.26

mkdir -p bowtie-index

bowtie2-build references/test_AAV.fa bowtie-index/test_AAV

mkdir -p blast-index

makeblastdb -in references/test_human.fa -dbtype nucl -out blast-index/test_human
makeblastdb -in references/test_AAV.fa -dbtype nucl -out blast-index/test_AAV
