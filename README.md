# stringMLST

Fast k-mer based tool for multi locus sequence typing (MLST)
stringMLST is a tool for detecting the MLST of an isolate directly from the genome sequencing reads. stringMLST predicts the ST of an isolate in a completely assembly and alignment free manner. The tool is designed in a light-weight, platform-independent fashion with minimum dependencies.

Some portions of the allele selection algorithm in stringMLST are patent
pending.  Please refer to the PATENTS file for additional inforamation
regarding licencing and use.


Reference  
*http://jordan.biology.gatech.edu/page/software/stringmlst/*

Abstract  
*http://bioinformatics.oxfordjournals.org/content/early/2016/09/06/bioinformatics.btw586.short?rss=1*

Application Note  
*http://bioinformatics.oxfordjournals.org/content/early/2016/09/06/bioinformatics.btw586.full.pdf+html*

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/stringmlst/README.html)  [![PyPI version](https://badge.fury.io/py/stringMLST.svg)](https://badge.fury.io/py/stringMLST)  ![downloads](https://img.shields.io/conda/dn/bioconda/stringmlst.svg?style=flat) [![container ready](https://quay.io/repository/biocontainers/stringmlst/status)](https://quay.io/repository/biocontainers/stringmlst)



**stringMLST is a *tool* not a *database*, always use the most up-to-date database files as possible.** To facilitate
keeping your databases updated, stringMLST can download and build databases from pubMLST using the most recent allele
and profile definitions. Please see the "Included databases and automated retrieval of databases from pubMLST" section
below for instructions. *The databases bundled here are for convenience only, do not rely on them being up-to-date*.

stringMLST is licensed and distributed under [CC Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0)  
and is free for academic users and requires permission before any commercial use for any version of this code/algorithm.  
If you are a commercial user, please contact king.jordan@biology.gatech.edu for permissions

## Recommended installation method

```
pip install stringMLST

```

#### Installation via git (Not recommended for most users)

```
git clone https://github.com/jordanlab/stringMLST
# Optional, download prebuilt databases  
# We don't recommend this method, instead build the databases locally
cd stringMLST
git submodule init
git submodule update
```

## Quickstart guide

```bash  
pip install stringMLST  
mkdir -p stringMLST_analysis; cd stringMLST_analysis  
stringMLST.py --getMLST -P neisseria/nmb --species neisseria  
# Download all available databases with:  
# stringMLST.py --getMLST -P mlst_dbs --species all    
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026529/ERR026529_2.fastq.gz  
stringMLST.py --predict -P neisseria/nmb -1 ERR026529_1.fastq.gz -2 ERR026529_2.fastq.gz
Sample  abcZ    adk     aroE    fumC    gdh     pdhC    pgm     ST
ERR026529       231     180     306     612     269     277     260     10174

```

## Python dependencies and external programs

stringMLST does not require any python dependencies for basic usage (Building databases and predicting STs). 

For advanced used (genome coverage), stringMLST depends on the `pyfaidx` python module and `bamtools`, `bwa`, and `samtools`.  
See the coverage section for more information

stringMLST has been tested with:
```
pyfaidx: 0.4.8.1
samtools: 1.3 (Using htslib 1.3.1)  [Requires the 1.x branch of samtools]
bedtools: v2.24.0
bwa: 0.7.13-r1126
```

### To install the dependencies

```bash
# pyfaidx
pip install --user pyfaidx
# samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o samtools-1.3.1.tar.bz2
tar xf samtools-1.3.1.tar.bz2
cd samtools-1.3.1.tar
make
make prefix=$HOME install
# bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2; make
cp ./bin/* ~/bin
# bwa
git clone https://github.com/lh3/bwa.git
cd bwa; make
cp bwa ~/bin/bwa
export PATH=$PATH:$HOME/bin
```


## Usage for Example Read Files (Neisseria meningitidis)

* Download stringMLST.py, example read files (ERR026529, ERR027250, ERR036104) and the dataset for Neisseria meningitidis (Neisseria_spp.zip).
### Build database:

```
# Add dir to path
export PATH=$PATH:$PWD
# Will connect to EBI's SRA servers
download_example_reads.sh
````

* Extract the MLST loci dataset.

```
unzip datasets/Neisseria_spp.zip -d datasets
```

* Create or use a config file specifying the location of all the locus and profile files.
Example config file (Neisseria_spp/config.txt):

```
[loci]
abcZ  datasets/Neisseria_spp/abcZ.fa
adk datasets/Neisseria_spp/adk.fa
aroE  datasets/Neisseria_spp/aroE.fa
fumC  datasets/Neisseria_spp/fumC.fa
gdh datasets/Neisseria_spp/gdh.fa
pdhC  datasets/Neisseria_spp/pdhC.fa
pgm datasets/Neisseria_spp/pgm.fa
[profile]
profile datasets/Neisseria_spp/neisseria.txt
```

* Run stringMLST.py --buildDB to create DB. Choose a k value and prefix (optional).

```
stringMLST.py --buildDB -c databases/Neisseria_spp/config.txt -k 35 -P NM 
```

### Predict:

#### Single sample :
```
stringMLST.py --predict -1 tests/fastqs/ERR026529_1.fastq -2 tests/fastqs/ERR026529_2.fastq -k 35 -P NM   
```
#### Batch mode (all the samples together):
```
stringMLST.py --predict -d ./tests/fastqs/ -k 35 -P NM  
```
#### List mode:
Create a list file (list_paired.txt) as :
```
tests/fastqs/ERR026529_1.fastq  tests/fastqs/ERR026529_2.fastq
tests/fastqs/ERR027250_1.fastq  tests/fastqs/ERR027250_2.fastq
tests/fastqs/ERR036104_1.fastq  tests/fastqs/ERR036104_2.fastq
```
Run the tool as:
```
stringMLST.py --predict -l list_paired.txt -k 35 -P NM
```
#### Working with gziped files
```
stringMLST.py --predict -1 tests/fastqs/ERR026529_1.fq.gz -2 tests/fastqs/ERR026529_2.fq.gz -p -P NM -k 35 -o ST_NM.txt
```
## Usage Documentation

stringMLST's workflow is divided into two routines:
* Database building and
* ST discovery

*Database building:* Builds the stringMLST database which is used for assigning STs to input sample files. This step is required once for each organism. Please note that stringMLST is capable of working on a custom user defined typing scheme but its efficiency has not been tested on other typing scheme.

*ST discovery:* This routine takes the database created in the last step and predicts the ST of the input sample(s). Please note that the database building is required prior to this routine. stringMLST is capable of processing single-end and paired-end files. It can run in three modes:
* Single sample mode - for running stringMLST on a single sample
* Batch mode - for running stringMLST on all the FASTQ files present in a directory
* List mode - for running stringMLST on all the FASTQ files provided in a list file


```
Readme for stringMLST
=============================================================================================
Usage
./stringMLST.py 
[--buildDB]
[--predict]
[-1 filename_fastq1][--fastq1 filename_fastq1]
[-2 filename_fastq2][--fastq2 filename_fastq2]
[-d directory][--dir directory][--directory directory]
[-l list_file][--list list_file]
[-p][--paired]
[-s][--single]
[-c][--config]
[-P][--prefix]
[-z][--fuzzy]
[-a]
[-C][--coverage]
[-k]
[-o output_filename][--output output_filename]
[-x][--overwrite]
[-t]
[-r]
[-v]
[-h][--help]
==============================================================================================

There are two steps to predicting ST using stringMLST.
1. Create DB : stringMLST.py --buildDB
2. Predict : stringMLST --predict

1. stringMLST.py --buildDB

Synopsis:
stringMLST.py --buildDB -c <config file> -k <kmer length(optional)> -P <DB prefix(optional)>
  config file : is a tab delimited file which has the information for typing scheme ie loci, its multifasta file and profile definition file.
    Format : 
      [loci]
      locus1    locusFile1
      locus2    locusFile2
      [profile]
      profile   profileFile
  kmer length : is the kmer length for the db. Note, while processing this should be smaller than the read length.
    We suggest kmer lengths of 35, 66 depending on the read length.
  DB prefix(optional) : holds the information for DB files to be created and their location. This module creates 3 files with this prefix.
    You can use a folder structure with prefix to store your db at particular location.

Required arguments
--buildDB
  Identifier for build db module
-c,--config = <configuration file>
  Config file in the format described above. 
  All the files follow the structure followed by pubmlst. Refer extended document for details. 

Optional arguments  
-k = <kmer length>
  Kmer size for which the db has to be formed(Default k = 35). Note the tool works best with kmer length in between 35 and 66
  for read lengths of 55 to 150 bp. Kmer size can be increased accordingly. It is advised to keep lower kmer sizes 
  if the quality of reads is not very good.
-P,--prefix = <prefix>
  Prefix for db and log files to be created(Default = kmer). Also you can specify folder where you want the dbb to be created.
-a
        File location to write build log
-h,--help
  Prints the help manual for this application

 --------------------------------------------------------------------------------------------
 
2. stringMLST.py --predict
  
stringMLST --predict : can run in three modes
  1) single sample (default mode)
  2) batch mode : run stringMLST for all the samples in a folder (for a particular specie)
  3) list mode : run stringMLST on samples specified in a file
stringMLST can process both single and paired end files. By default program expects paired end files.

Synopsis
stringMLST.py --predict -1 <fastq file> -2 <fastq file> -d <directory location> -l <list file> -p -s -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x

Required arguments
--predict
  Identifier for predict miodule
  
Optional arguments
-1,--fastq1 = <fastq1_filename>
  Path to first fastq file for paired end sample and path to the fastq file for single end file.
  Should have extension fastq or fq.
-2,--fastq2 = <fastq2_filename>
  Path to second fastq file for paired end sample.
  Should have extension fastq or fq.
-d,--dir,--directory = <directory>
  BATCH MODE : Location of all the samples for batch mode.
-C,--coverage
  Calculate seqence coverage for each allele. Turns on read generation (-r) and turns off fuzzy (-z 1)
  Requires bwa, bamtools and samtools be in your path
-k = <kmer_length>
  Kmer length for which the db was created(Default k = 35). Could be verified by looking at the name of the db file. 
  Could be used if the reads are of very bad quality or have a lot of N's.
-l,--list = <list_file>
  LIST MODE : Location of list file and flag for list mode.
  list file should have full file paths for all the samples/files.
  Each sample takes one line. For paired end samples the 2 files should be tab separated on single line.
-o,--output = <output_filename>
  Prints the output to a file instead of stdio.
-p,--paired
  Flag for specifying paired end files. Default option so would work the same if you do not specify for all modes.
  For batch mode the paired end samples should be differentiated by 1/2.fastq or 1/2.fq
-P,--prefix = <prefix>
  Prefix using which the db was created(Defaults = kmer). The location of the db could also be provided.
-r
  A seperate reads file is created which has all the reads covering all the locus.
-s,--single
  Flag for specifying single end files.
-t
  Time for each analysis will also be reported.
-v
  Prints the version of the software.
-x,--overwrite
  By default stringMLST appends the results to the output_filename if same name is used.
  This argument overwrites the previously specified output file.
-z,--fuzzy = <fuzzy threshold int>
  Threshold for reporting a fuzzy match (Default=300). For higher coverage reads this threshold should be set higher to avoid
  indicating fuzzy match when exact match was more likely. For lower coverage reads, threshold of <100 is recommended
-h,--help
  Prints the help manual for this application

 --------------------------------------------------------------------------------------------

3. stringMLST.py --getMLST

Synopsis:  
stringMLST.py --getMLST --species= <species> [-k kmer length] [-P DB prefix]

Required arguments
--getMLST
    Identifier for getMLST module
--species= <species name>
    Species name from the pubMLST schemes (use --schemes to get list of available schemes)
    "all" will download and build all 

Optional arguments
-k = <kmer length>
    Kmer size for which the db has to be formed(Default k = 35). Note the tool works best with kmer length in between 35 and 66
    for read lengths of 55 to 150 bp. Kmer size can be increased accordingly. It is advised to keep lower kmer sizes
    if the quality of reads is not very good.
-P,--prefix = <prefix>
    Prefix for db and log files to be created(Default = kmer). Also you can specify folder where you want the db to be created.
    We recommend that prefix and config point to the same folder for cleanliness but this is not required
--schemes
    Display the list of available schemes
-h,--help
  Prints the help manual for this application

```


**stringMLST expects paired end reads to be in [Illumina naming convention](http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm), minimally ending with _1.fq and _2.fq to delineate read1 and read2:**

*Periods (.) are disallowed delimiters except for file extensions*

```
Illumina FASTQ files use the following naming scheme:

<sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz

For example, the following is a valid FASTQ file name:

NA10831_ATCACG_L002_R1_001.fastq.gz
```

## Running stringMLST

#### Included databases and automated retrieval of databases from pubMLST

stringMLST includes all the pubMLST databases as of **February 15, 2017**, built with the default kmer (*35*). They can be found in the `datasets/` folder.  
Simply unzip the databases you need and begin using stringMSLT as described below.

All the databases from pubMLST can be downloaded and prepared with your kmer choice

*Getting all pubMLST schemes*  
```
stringMLST.py --getMLST -P datasets/ --species all
```


Individual databases from pubMLST can also be downloaded as needed, using the scheme identifiers

*Downloading a scheme*  
```
# List available schemes
stringMLST.py --getMLST --schemes

# Download the Neisseria spp. scheme

stringMLST.py --getMLST -P datasets/nmb --species neisseria

```



#### Database Preparation
In order to create the database, files can be downloaded from the database page.

If the organism of interest is not present in the provided link, the required files can be downloaded from PubMLST as follows:
* On your browser, navigate to http://pubmlst.org/
* Navigate to "Download MLST definitions" link or go to http://pubmlst.org/data/
* Scroll to the species of interest. For each species, user may find the file for typing definitions and multi-FASTA files for each locus. Download these files.

E.g.:

Species of interest: Neisseria spp.
Corresponding definition file: http://pubmlst.org/data/profiles/neisseria.txt
Corresponding multi fasta locus files: 
http://pubmlst.org/data/alleles/neisseria/abcZ.tfa
http://pubmlst.org/data/alleles/neisseria/adk.tfa
http://pubmlst.org/data/alleles/neisseria/aroE.tfa
http://pubmlst.org/data/alleles/neisseria/fumC.tfa
http://pubmlst.org/data/alleles/neisseria/gdh.tfa
http://pubmlst.org/data/alleles/neisseria/pdhC.tfa
http://pubmlst.org/data/alleles/neisseria/pgm.tfa

Download these files at a desired location.


Custom user files can also be used for building database. The database building routine requires the profile definition file and allele sequence file. The profile definition file is a tab separated file that contains the ST and the allele profile corresponding to the ST. An example of the profile definition file is shown below:
```
ST  abcZ  adk aroE  fumC  gdh pdhC  pgm clonal_complex
1 1 3 1 1 1 1 3 ST-1 complex/subgroup I/II
2 1 3 4 7 1 1 3 ST-1 complex/subgroup I/II
3 1 3 1 1 1 23  13  ST-1 complex/subgroup I/II
4 1 3 3 1 4 2 3 ST-4 complex/subgroup IV
```
The allele sequence file is a standard multi-FASTA with the description being the loci name with the allele number. An example abcZ allele sequence is shown below:
```
>abcZ_1
TTTGATACTGTTGCCGA...
>abcZ_2
TTTGATACCGTTGCCGA...
>abcZ_3
TTTGATACCGTTGCGAA...
>abcZ_4
TTTGATACCGTTGCCAA...
```

These files can be obtained from PubMLST/BIGSdb or can be create by the user themselves. 

In either case, an accompanying configuration file is also required to describe the profile definition and allele sequence files. An example configuration file is shown below:
```
[loci]
abcZ  /data/home/stringMLST/pubmlst/Neisseria_sp/abcZ.fa
adk /data/home/stringMLST/pubmlst/Neisseria_sp/adk.fa
aroE  /data/home/stringMLST/pubmlst/Neisseria_sp/aroE.fa
fumC  /data/home/stringMLST/pubmlst/Neisseria_sp/fumC.fa
gdh /data/home/stringMLST/pubmlst/Neisseria_sp/gdh.fa
pdhC  /data/home/stringMLST/pubmlst/Neisseria_sp/pdhC.fa
pgm /data/home/stringMLST/pubmlst/Neisseria_sp/pgm.fa

[profile]
profile /data/home/stringMLST/pubmlst/Neisseria_sp/neisseria.txt
```

This file is pre-packed on stringMLSTs website and can easily be created by the user for custom database.

#### Database Building 
The next step is for database building is running the buildDB module to create the database files. buildDB module requires the user to specify the config file. The default k-mer size is 35 but can be changed using the -k option. Specifying the prefix for the created database files is optional but is recommended.

The choice of k-mer depends on the size of the sequencing read. In general, the value of k can never be greater than the read length. The application has been tested on a number of read lengths ranging from 55 to 150 bps using k-mer sizes of 21 to 66. In our testing, the k-mer size does not affect the accuracy of the read length. A smaller k-mer size will increase the runtime and a larger k-mer size will increase the file size. The user should ideally pick a k-mer with a length around half of the average read length. For lower quality data, it also advised to choose smaller k-mer values to reduce false hits.
```
stringMLST.py --buildDB --config <config file> -k  <k-mer length> -P <prefix>
```
Example:
```
stringMLST.py --buildDB --config config.txt -k 35 -P NM
```
This command will produce 3 database files and a log file. The log file is used for debugging purposes in the event an error is encountered. The 3 database files created are:
* <prefix>_<k-mer>.txt : The main database file for the application. This is a tab delimited file describing k-mer to locus relationship.
* <prefix>_weight.txt : Contains the weight factors for alleles which differ in lengths by more than 5%. Will be empty otherwise.
* <prefix>_profile.txt : Profile definition file used for finding the ST from the predicted allelic profile.

For the example above, the following files will be created:
NM_35.txt, NM_weight.txt and NM_profile.txt

Please note that in the prediction routine the database is identified with the prefix. 

ST discovery routine
As discussed earlier, StringMLST has 3 running modes 
* Single sample mode - for running stringMLST on a single sample
* Batch mode - for running stringMLST on all the FASTQ files present in a directory
* List mode - for running stringMLST on all the FASTQ files provided in a list file

####  Single sample mode: 
This is the default mode for stringMLST and takes in one sample at a time. The sample can be single-end or paired-end. The sample has to be in FASTQ format. In order to run, the user should know the prefix of the database created and the k-mer size. 

By default, the tool expects paired-end samples.
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
*For single-end samples:*
```
stringMLST.py --predict -1 <single-end file> -s --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
####  Batch Mode: 
This mode can be used for processing multiple files with one command. All the samples will be queried against the same database. Also all samples should be in the same directory. All the samples will be treated either as single-end or paired-end. The paired-end samples should be differentiated with the character _1 and _2 at the end (E.g.: sampleX_1.fastq and sampleX_2.fastq).

*Paired-end samples:*
```
stringMLST.py --predict -d <directory for samples> -p --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```

*Single-end samples:*
```
stringMLST.py --predict -d <directory for samples> -s --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
#### List Mode: 
This mode could be used if user has samples at different locations or if the paired-end samples are not stored in traditional way. All the samples will be queried against the same database. All the samples will be treated either as single-end or paired-end. This mode requires the user to provide a list file which has the list of all samples along with the location. Each line in the list file represents a new sample.
A sample list file for single-end sample looks like the following.
```
<full path of sample 1 fastq file>
<full path of sample 2 fastq file>
<full path of sample 3 fastq file>
.
.
<full path of sample n fastq file>
```
A sample list file for paired-end sample looks like the following.

```
<full path of sample 1 fastq file 1>  <full path of sample 1 fastq file 2>
<full path of sample 2 fastq file 1>  <full path of sample 2 fastq file 2>
<full path of sample 3 fastq file 1>  <full path of sample 3 fastq file 2>
.
.
<full path of sample n fastq file 1>  <full path of sample n fastq file 2>
```

Once the user has the list file, he can directly use the tool.

*Paired-end samples:*
```
stringMLST.py --predict -l <full path to list file> -p --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
*Single-end samples:*
```
stringMLST.py --predict -l <full path to list file > -s --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```

#### Gene coverage and match confidence

stringMLST provides two, complimentary methods for determining confidence in an inferred ST. There's the `-C|--coverage` flag and `-z|--fuzzy` threshold option.

stringMLST determines an allele based on its kmer support; the more kmers seen for allele 1, the more likely that allele 1 is the allele present in the genome. Unlike SRST2 and other mapping/BLAST based tools, stringMLST always infers an ST, using the maximimally supported allele (allele with most kmer hits). The difference between the maximum support (the reported allele) and the second support (next closest allele) can be informative for low coverage reads. The `-z|--fuzzy` threshold (Default = 300), assigns significance to the difference between supports. Much like SRST2 and Torsten Seemann's popular [pubMLST script](https://github.com/tseemann/mlst), stringMLST reports potentially new or closely supported alleles in allele* syntax. For high coverage reads, we suggest a fuzzy threshold >500. For low coverage reads, a fuzzy threshold of <50.

Coverage mode requires `bedtools`, `bwa`, and `samtools` in your PATH and an additional python module, `pyfaidx` (See the dependencies section for installion information).  Coverage mode by default disables display of fuzzy alleles in favor of sequence coverage information made by mapping potential reads to the putative allele sequence. In our testing, coverage mode slightly increases prediction time (<1 sec increase per sample). 

**Please note:** stringMLST *always* infers the ST from the reads, fuzzy matches and/or <100% coverage do not necessarily mean a new allele has been found.

*Getting gene coverage from reads*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -r -o <output file name>- -c <path to config> -C
```
*Changing the fuzziness of the search for low coverage reads*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -r -o <output file name>- -f 50
```

#### Other Examples : 

*Reporting time along with the output.*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -t -o <output file name>
```
*Getting reads file relevant to typing scheme.*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -r -o <output file name>
```
