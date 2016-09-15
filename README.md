# stringMLST
Fast k-mer based tool for multi locus sequence typing (MLST)
stringMLST is a tool for detecting the MLST of an isolate directly from the genome sequencing reads. stringMLST predicts the ST of an isolate in a completely assembly and alignment free manner. The tool is designed in a light-weight, platform-independent fashion with minimum dependencies.

Reference
*http://jordan.biology.gatech.edu/page/software/stringmlst/*

Abstract
*http://bioinformatics.oxfordjournals.org/content/early/2016/09/06/bioinformatics.btw586.short?rss=1*

Application Note
*http://bioinformatics.oxfordjournals.org/content/early/2016/09/06/bioinformatics.btw586.full.pdf+html*

## Usage for Example Read Files (Neisseria meningitidis)

* Download stringMLST.py, example read files (ERR026529, ERR027250, ERR036104) and the dataset for Neisseria meningitidis (Neisseria_spp.zip).
### Build database:

* Extract the dataset such that the locus and profile file are in Neisseria_spp folder.
* Create or use a config file specifying the location of all the locus and profile files.
Example config file (Neisseria_spp/config.txt):
```
[loci]
abcZ	Neisseria_spp/abcZ.fa
adk	Neisseria_spp/adk.fa
aroE	Neisseria_spp/aroE.fa
fumC	Neisseria_spp/fumC.fa
gdh	Neisseria_spp/gdh.fa
pdhC	Neisseria_spp/pdhC.fa
pgm	Neisseria_spp/pgm.fa
[profile]
profile	Neisseria_spp/neisseria.txt
```

* Run stringMLST.py --buildDB to create DB. Choose a k value and prefix (optional).
```
./stringMLST.py --buildDB -c Neisseria_spp/config.txt -k 35 -P NM	
```

### Predict:

* Extract all the read files to a folder (say ./example/) such that fastq files are in the folder as :
```
example/ERR026529_1.fastq
example/ERR026529._2fastq
example/ERR027250_1.fastq
example/ERR027250_2.fastq
example/ERR036104_1.fastq
example/ERR036104_2.fastq
```
* Predict:

#### Single sample :
```
./stringMLST.py --predict -1 example/ERR026529_1.fastq -2 example/ERR026529_2.fastq -k 35 -P NM 	
```
#### Batch mode (all the samples together):
```
./stringMLST.py --predict -d ./example/ -k 35 -P NM	
```
#### List mode:
Create a list file (list_paired.txt) as :
```
example/ERR026529_1.fastq	example/ERR026529_2.fastq
example/ERR027250_1.fastq	example/ERR027250_2.fastq
example/ERR036104_1.fastq	example/ERR036104_2.fastq
```
Run the tool as:
```
./stringMLST.py --predict -l list_paired.txt -k 35 -P NM
```
#### Working with gziped files
```
stringMLST.py --predict -1 example/ERR026529_1.fq.gz -2 example/ERR026529_2.fq.gz -p --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
## Usage Documentation

stringMLST's workflow is divided into two routines:
* Database building and
* ST discovery

*Database building:* Builds the stringMLST database which is used for assigning STs to input sample files. This step is required once for each organism. Please note that stringMLST is capable of working on a custom user defined typing scheme but its efficiency has not been tested on other typing scheme.

*ST discovery:* This routine takes the database created in the last step and predicts the ST of the input sample(s). Please note that the database building is required prior to this routine. stringMLST is capable of processing single-end and paired-end files. It can run in three modes:
*	Single sample mode - for running stringMLST on a single sample
*	Batch mode - for running stringMLST on all the FASTQ files present in a directory
*	List mode - for running stringMLST on all the FASTQ files provided in a list file

### Installation

stringMLST only requires Python 2.6+ to run. No additional dependencies are required. Simply place the script in one of your PATH directories and provide it executable permissions. The script can then be run by executing stringMLST.py command.


### Running stringMLST

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
ST	abcZ	adk	aroE	fumC	gdh	pdhC	pgm	clonal_complex
1	1	3	1	1	1	1	3	ST-1 complex/subgroup I/II
2	1	3	4	7	1	1	3	ST-1 complex/subgroup I/II
3	1	3	1	1	1	23	13	ST-1 complex/subgroup I/II
4	1	3	3	1	4	2	3	ST-4 complex/subgroup IV
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
abcZ	/data/home/stringMLST/pubmlst/Neisseria_sp/abcZ.fa
adk	/data/home/stringMLST/pubmlst/Neisseria_sp/adk.fa
aroE	/data/home/stringMLST/pubmlst/Neisseria_sp/aroE.fa
fumC	/data/home/stringMLST/pubmlst/Neisseria_sp/fumC.fa
gdh	/data/home/stringMLST/pubmlst/Neisseria_sp/gdh.fa
pdhC	/data/home/stringMLST/pubmlst/Neisseria_sp/pdhC.fa
pgm	/data/home/stringMLST/pubmlst/Neisseria_sp/pgm.fa

[profile]
profile	/data/home/stringMLST/pubmlst/Neisseria_sp/neisseria.txt
```

This file is pre-packed on stringMLSTs website and can easily be created by the user for custom database.

### Database Building 
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
*	<prefix>_<k-mer>.txt : The main database file for the application. This is a tab delimited file describing k-mer to locus relationship.
*	<prefix>_weight.txt : Contains the weight factors for alleles which differ in lengths by more than 5%. Will be empty otherwise.
*	<prefix>_profile.txt : Profile definition file used for finding the ST from the predicted allelic profile.

For the example above, the following files will be created:
NM_35.txt, NM_weight.txt and NM_profile.txt

Please note that in the prediction routine the database is identified with the prefix. 

ST discovery routine
As discussed earlier, StringMLST has 3 running modes 
*	Single sample mode - for running stringMLST on a single sample
*	Batch mode - for running stringMLST on all the FASTQ files present in a directory
*	List mode - for running stringMLST on all the FASTQ files provided in a list file

####	Single sample mode: 
This is the default mode for stringMLST and takes in one sample at a time. The sample can be single-end or paired-end. The sample has to be in FASTQ format. In order to run, the user should know the prefix of the database created and the k-mer size. 

By default, the tool expects paired-end samples.
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
*For single-end samples:*
```
stringMLST.py --predict -1 <single-end file> -s --prefix <prefix for the database> -k <k-mer size> -o <output file name>
```
####	Batch Mode: 
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
<full path of sample 1 fastq file 1>	<full path of sample 1 fastq file 2>
<full path of sample 2 fastq file 1>	<full path of sample 2 fastq file 2>
<full path of sample 3 fastq file 1>	<full path of sample 3 fastq file 2>
.
.
<full path of sample n fastq file 1>	<full path of sample n fastq file 2>
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
#### Other Examples : 

*Reporting time along with the output.*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -t -o <output file name>
```
*Getting reads file relevant to typing scheme.*
```
stringMLST.py --predict -1 <paired-end file 1> -2 <paired-end file 2> -p --prefix <prefix for the database> -k <k-mer size> -r -o <output file name>
```
