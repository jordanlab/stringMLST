#!/bin/bash

testdir='tests/fastqs'

echo -e "Downloading fastq files from EBI..."

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_1.fastq.gz -P $testdir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_2.fastq.gz -P $testdir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR027/ERR027250/ERR027250_1.fastq.gz -P $testdir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR027/ERR027250/ERR027250_2.fastq.gz -P $testdir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036104/ERR036104_1.fastq.gz -P $testdir
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036104/ERR036104_2.fastq.gz -P $testdir

echo -e "Done downloading.\n\n"

for file in `ls $testdir/*` ; do
    echo -e "Attempting to unzip $file..."
    gunzip $file;
    echo -e "Done unzipping $file.\n"
done

echo -e "Done downloading and extrating test read files."
