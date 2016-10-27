#!/usr/bin/env python
import getopt
import sys
import logging
import os
import time
import ast
import gzip
import re
import tempfile
import shutil
v = """ stringMLST v0.3 (updated : October 26, 2016) """
"""
LICENSE TERMS FOR stringMLST
1. INTENT/PURPOSE:
stringMLST is free for academic users and requires permission before any commercial use for
any version of this code/algorithm.
If you are a commercial user, please contact king.jordan@biology.gatech.edu for permissions.
2. LICENSEE:
Any person or organization who receives the software with a copy of this license.
3. INTELLECTUAL PROPERTY LICENSED:
The rights to use, copy, modify, and compile the software,
subject to the restrictions described in the present document.
5. SCOPE OF THE LICENSE
* NON-COMMERCIAL license for research and evaluation purposes ONLY.
* NO right to commercialize the software, or any derivative work,
 without separate agreement with the copyright owners.
6. MODIFICATION
* License permits licensee to modify the software for their
 research and evaluation purposes.
7. REDISTRIBUTION
* License permits licensee to redistribute verbatim copies
 of the software, accompanied with a copy of this license.
* License does not permit licensee to redistribute
 modified versions of the software.
* License does not permit licensee to commercialize the software
 or any derivative work of the software.
8. FEE/ROYALTY
* Licensee pays no royalty for non-commercial license
* Licensee and any third parties must enter a new agreement
 for any use beyond scope of license.
9. NO WARRANTY
The software is provided "as is" without warranty of any kind,
either expressed or implied, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose.
The entire risk as to the quality and performance of
the program is with the Licensee.
10. NO LIABILITY
In no event unless required by applicable law or agreed to in writing
will any copyright owner be liable to Licensee for damages, including any
general, special, incidental or consequential damages arising out of the
use or inability to use the program (including but not limited to loss of
data or data being rendered inaccurate or losses sustained by Licensee
or third parties or a failure of the program to operate with other programs),
even if such holder has been advised of the possibility of such damages.


The program has 3 basic modes :
    mainTool: for single sample (both single and paired end)
    batchTool: for multiple samples stored at a common location (both single and paired end samples)
    listTool: for multiple samples with location information stored in a list (both single and paired end samples)
predict part starts here
"""
#############################################################
# Function   : batchTool
# Input      : Directory name, paired or single, k value
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes all FASTQ files present in the input
#              directory
#############################################################
def batchTool(fdir, paired, k):
    fileList = []
    if not dir.endswith('/'):
        fdir += '/'
    for inputFile in os.listdir(fdir):
        if paired is True:
            if inputFile.endswith('1.fastq') or inputFile.endswith('1.fq') or inputFile.endswith('1.fq.gz') or inputFile.endswith('1.fastq.gz'):
                fastq1 = fdir+inputFile
                fastq2 = fdir+inputFile.replace('1.', '2.')
                fileList.append((fastq1, fastq2))
        else:
            if inputFile.endswith('.fastq') or inputFile.endswith('.fq') or inputFile.endswith('.fq.gz') or inputFile.endswith('.fastq.gz'):
                fastq1 = fdir + inputFile
                fileList.append(fastq1)
    results = multiSampleTool(fileList, paired, k)
    return results
#############################################################
# Function   : listTool
# Input      : List file, paired or single, k value
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes all FASTQ files present in the input
#              list file
#############################################################
def listTool(fList, paired, k):
    fileList = []
    listf = open(fList, 'r')
    samples = listf.readlines()
    for sample in samples:
        if paired is True:
            s = sample.strip().split()
            fastq1 = s[0]
            try:
                fastq2 = s[1]
            except IndexError:
                print("Error: Paired end files should be whitespace/tab seperated")
                exit(0)
            fileList.append((fastq1, fastq2))
        else:
            fastq1 = sample.rstrip()
            fileList.append(fastq1)
    results = multiSampleTool(fileList, paired, k)
    return results
#############################################################
# Function   : multiSampleTool
# Input      : List of files to process, paired or single, k value
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes all FASTQ files present in the input list
#############################################################
def multiSampleTool(fileList, paired, k):
    results = {}
    for sample in fileList:
        if paired is True:
            fastq1 = sample[0]
            fastq2 = sample[1]
        else:
            fastq1 = sample
            fastq2 = None
        results = singleSampleTool(fastq1, fastq2, paired, k, results)
    return results
#############################################################
# Function   : singleSampleTool
# Input      : fastq file 1 and 2, paired or single, k value, output dictionary
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes both FASTQ files passed to the function
#############################################################
def singleSampleTool(fastq1, fastq2, paired, k, results):
    if paired is True:
        fileName = fastq1.split('/')[-1].split('.')[0][:-1]
    else:
        fileName = fastq1.split('/')[-1].split('.')[0]
    if reads is True:
        readFileName = fileName + '_reads.fq'
        global readFile
        readFile = open(readFileName, 'w+')
    if paired is True:
        msg = "singleSampleTool : " + fastq1 + ' and ' + fastq2
    else:
        msg = "singleSampleTool : " + fastq1
    logging.debug(msg)
    global alleleCount
    alleleCount = {}
    t1 = time.time()
    if paired is True:
        logging.debug("singleSampleTool : paired True")
        logging.debug("singleSampleTool : fastq1 start")
        singleFileTool(fastq1, k)
        logging.debug("singleSampleTool : fastq1 done")
        logging.debug("singleSampleTool : fastq2 start")
        singleFileTool(fastq2, k)
        logging.debug("singleSampleTool : fastq2 done")
        if alleleCount == {}:
            string = "No k-mer matches were found for the sample " + fastq1 + " and "+ fastq2 + ".  Probable cause of the error:  low quality data/too many N's in the data"
            logging.error("singleSampleTool : " + string)
            print(string)
#           exit(0)
        profileCount = alleleCount
    else:
        logging.debug("singleSampleTool : paired False")
        logging.debug("singleSampleTool : fastq start")
        singleFileTool(fastq1, k)
        profileCount = alleleCount
        logging.debug("singleSampleTool : fastq done")
        if alleleCount == 0:
            string = "No k-mer matches were found for the sample " + fastq1 + ".  Probable cause of the error:  low quality data/too many N's in the data"
            logging.error("singleSampleTool : " + string)
            print(string)
    logging.debug("singleSampleTool : weightedProfile start")
    weightedProfile = weightedProf(profileCount, weightDict)
    logging.debug("singleSampleTool : weightedProfile finished")
    logging.debug("singleSampleTool : getMaxCount start")
    finalProfile = getMaxCount(weightedProfile, fileName)
    logging.debug("singleSampleTool : getMaxCount end")
    st = 0
    if profileFile != '':
        logging.debug("singleSampleTool : findST start")
        st = findST(finalProfile, stProfile)
        logging.debug("singleSampleTool : findST end")
    if reads is True:
        readFile.close()
    t3 = time.time()
    finalProfile['ST'] = st
    finalProfile['t'] = t3-t1
    results[fileName] = finalProfile
    return results
#############################################################
# Function   : singleFileTool
# Input      : fastq file, k value
# Output     : Edits a global dictionary - results
# Description: Processes the single fastq file
#############################################################
def singleFileTool(fastq, k):
    msg = "singleFileTool :" + fastq
    logging.debug(msg)
    if os.path.isfile(fastq):
        logging.debug("singleFileTool : fastq")
        non_overlapping_window = 1
        finalProfile = {}
        t1 = time.time()
        fileExplorer(fastq, k, non_overlapping_window)
        t3 = time.time()
    else:
        msg = "File does not exist: " + fastq
        logging.error("singleFileTool : msg")
        print(msg)
def fileExplorer(file, k, non_overlapping_window):
    if file.endswith('.gz'):
        f = gzip.open(file, 'rb')
    else:
        f = open(file)
    msg = "fileExplorer :" + file
    logging.debug(msg)
    lines = f.readlines()
    i = 1
    n_reads = 0
    try:
        if len(lines[1]) < k:
            m1 = "Read length " + len(lines[1])+" for file " + file + " smaller than " + k
            print(m1)
            print("Skipping to next file.")
            logging.debug(m1)
            return 0
    except Exception:
        m2 = "Check fastq file " + file
        print(m2)
        logging.debug(m2)
        return 0
    start = int((len(lines[1])-k)//2)
    end = int((len(lines[1])-k)//2)
    yesRead = False
    for line in lines:
        if i % 4 == 0 and yesRead:
            readFile.write(line)
        if i % 4 != 3:
            yesRead = False
        if i%4 == 1:
            head = line
        if i%4 == 2:
            s1 = str(line[start:k+start])
            sn_1 = str(line[-k-end:-end]).rstrip()
            if s1 in kmerDict[k]:
                n_reads += 1
                goodReads(line, k, non_overlapping_window)
                if reads is True:
                    readFile.write(head)
                    readFile.write(line)
                    readFile.write('+\n')
                    yesRead = True
        i += 1
#############################################################
# Function   : goodReads
# Input      : sequence read, k, step size
# Output     : Edits the count of global variable alleleCount
# Description: Increment the count for each k-mer match
#############################################################
def goodReads(read, k, non_overlapping_window):
    n = 0
    line = read.rstrip()
    while n+k <= len(line):
        s = str(line[n:n+k])
        if s in kmerDict[k]:
            for probLoc in kmerDict[k][s]:
                if probLoc not in alleleCount:
                    alleleCount[probLoc] = {}
                a = kmerDict[k][s][probLoc]
                for allele in a:
                    allele = allele.rstrip()
                    if allele in alleleCount[probLoc]:
                        alleleCount[probLoc][allele] += 1
                    else:
                        alleleCount[probLoc][allele] = 1
        n += non_overlapping_window
#############################################################
# Function   : weightedProf
# Input      : allele count global var, weight factors
# Output/Desc: Normalizes alleleCount by weight factor
#############################################################
def weightedProf(alleleCount, weightDict):
    logging.debug("weightedProf")
    weightedDict = {}
    for loc in alleleCount:
        weightedDict[loc] = {}
        for allele in alleleCount[loc]:
            if loc in weightDict:
                if allele in weightDict[loc]:
                    weightedDict[loc][allele] = (alleleCount[loc][allele] / weightDict[loc][allele])
                else:
                    weightedDict[loc][allele] = alleleCount[loc][allele]
            else:
                weightedDict[loc][allele] = alleleCount[loc][allele]
    return weightedDict
#############################################################
# Function   : getMaxCount
# Input      : allele counts
# Output     : allelic profile and ST
# Description: Finds the alleles with maximum counts and
#              generates the allelic profile and ST
#############################################################
def getMaxCount(alleleCount, fileName):
    logging.debug("getMaxCount")
    max_n = {}
    secondMax = {}
    maxSupport = {}
    secondSupport = {}
    finalProfileCount = {}
    num = ''
    for loc in alleleCount:
        n = 0
        m = 0
        for num in alleleCount[loc]:
            if alleleCount[loc][num] >= n:
                m = n
                n = alleleCount[loc][num]
        if n-m < fuzzy:
            alleleCount[loc][num] = str(alleleCount[loc][num])+'*'
            max_n[loc] = str(n)+'*'
        else:
            max_n[loc] = n
        secondMax[loc] = m
    for loc in alleleCount:
        maxSupport[loc] = {}
        secondSupport[loc] = {}
        num_max = []
        num_max2 = []
        compare = int(re.sub("\*$", "", str(max_n[loc])))
        for num in alleleCount[loc]:
            if  alleleCount[loc][num] == compare:
                if "\*" in str(max_n[loc]):
                    insert = num + '*'
                    num_max.append(insert)
                else:
                    num_max.append(num)
                maxSupport[loc][num] = max_n[loc]
            if  alleleCount[loc][num] == secondMax[loc]:
                num_max2.append(num)
                secondSupport[loc][num] = secondMax[loc]
        try:
            finalProfileCount[loc] = num_max[0]
        except LookupError:
            finalProfileCount[loc] = '0'
    msgs = "Max Support :" + fileName + " : " + str(maxSupport)
    logging.debug(msgs)
    msgs = "Second Max Support :" + fileName + " : " + str(secondSupport)
    logging.debug(msgs)
    return finalProfileCount
#############################################################
# Function   : findST
# Input      : allelic profile for one sample and profiles for all STs
# Output     : ST number, or 0 if no ST match was found
# Description: Finds the ST number which best matches the given sample profile.
#############################################################
def findST(finalProfile, stProfile):
    if not stProfile:
        return 0
    oneProfile = next(iter(stProfile.values()))
    # The gene names in finalProfile may not exactly match those in stProfile. To deal with this,
    # each finalProfile gene is associated with the best matching gene in the ST profiles.
    finalGeneToSTGene = {}
    profileGenes = list(oneProfile.keys())
    for finalGene in list(finalProfile.keys()):
        if finalGene in profileGenes:  # exact match is preferable
            finalGeneToSTGene[finalGene] = finalGene
        else:  # failing an exact match, look for a case-sensitive containment
            for profileGene in profileGenes:
                if finalGene in profileGene:
                    finalGeneToSTGene[finalGene] = profileGene
                    break
        if finalGene not in finalGeneToSTGene:  # if there's still no match, try a case-insensitive containment
            for profileGene in profileGenes:
                if finalGene.lower() in profileGene.lower():
                    finalGeneToSTGene[finalGene] = profileGene
                    break
        if finalGene not in finalGeneToSTGene:
            print("ERROR: gene names in config file do not match gene names in profile file")
            exit(0)
    transformedFinalProfile = {}
    for gene, allele in finalProfile.items():
        allele = re.sub("\*", "", allele)
        transformedFinalProfile[finalGeneToSTGene[gene]] = allele
        # Check to see if the dictionary is empty, if so then means no allele were found at all
        if bool(transformedFinalProfile) is False:
            return 0
    # Find the best matching ST, considering only the genes in the sample's profile. This is to
    # allow for superfluous columns in the ST profile.
    logging.debug("findST")
    for stNum, profile in stProfile.items():
        if all(x in list(profile.items()) for x in list(transformedFinalProfile.items())):
            return stNum
    return 0
#############################################################
# Function   : loadModule
# Input      : k value and prefix of the DB file
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#              by calling other functions
#############################################################
def loadModule(k, dbPrefix):
    global dbFile
    dbFile = dbPrefix+'_'+str(k)+'.txt'
    global weightFile
    weightFile = dbPrefix+'_weight.txt'
    global profileFile
    profileFile = dbPrefix+'_profile.txt'
    global kmerDict
    kmerDict = {}
    kmerDict[k] = loadKmerDict(dbFile)
    global weightDict
    weightDict = loadWeightDict(weightFile)
    global stProfile
    stProfile = loadSTfromFile(profileFile)
#############################################################
# Function   : loadSTfromFile
# Input      : profile definition file
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadSTfromFile(profileF):
    with open(profileF, 'r') as definitionFile:
        st = {}
        index = {}
        lines = definitionFile.readlines()
        heads = lines[0].rstrip().split('\t')
        for locus in heads:
            index[locus] = heads.index(locus)
        for line in lines:
            pro = line.rstrip().split('\t')
            l = {}
            for locus in heads[1:]:
                try:
                    l[locus] = pro[index[locus]]
                except LookupError:
                    logging.debug("ERROR while loading ST")
                    pass
            st[pro[0]] = l
    return st
#############################################################
# Function   : loadKmerDict
# Input      : DB prefix
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadKmerDict(dbFile):
    kmerTableDict = {}
    with open(dbFile, 'r') as kmerTableFile:
        lines = kmerTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            kmerTableDict[array[0]] = {}
            kmerTableDict[array[0]][array[1]] = array[2][1:-1].rsplit(',')
    return kmerTableDict
#############################################################
# Function   : loadWeightDict
# Input      : Weight file prefix
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadWeightDict(weightFile):
    weightDict = {}
    with open(weightFile, 'r') as weightTableFile:
        lines = weightTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            loc = array[0].rsplit('_')[0]
            allele = array[0].rsplit('_')[1]
            if loc not in weightDict:
                weightDict[loc] = {}
            weightDict[loc][allele] = float(array[1])
    return weightDict
#############################################################
# Function   : loadConfig
# Input      : config file path from getopts
# Output     : Updates configDict
# Description: Used to find allele fasta files for getCoverage
#############################################################
def loadConfig(config):
    global configDict
    configDict = {}
    with open(config) as configFile:
        lines = configFile.readlines()
        head = ''
        for line in lines:
            if line.rstrip() == '':
                continue
            if line.rstrip() == '[loci]':
                head = 'loci'
                configDict[head] = {}
            elif line.rstrip() == '[profile]':
                head = 'profile'
                configDict[head] = {}
            else:
                arr = line.strip().split()
                configDict[head][arr[0]] = arr[1]
    for head in configDict:
        for element in configDict[head]:
            if not os.path.isfile(configDict[head][element]):
                print("ERROR: %s file does not exist at %s" % (element, configDict[head][element]))
                exit(0)
    return configDict
#############################################################
# Function   : getCoverage
# Input      : results dictionary
# Output     : Updates results to include coverage info
#############################################################
def getCoverage(results):
    tmpdir = tempfile.mkdtemp()
    for sample in results:
        file = tmpdir +'/'+ sample + '.fasta'
        bed = tmpdir +'/'+ sample + '.bed'
        sortedFile = tmpdir +'/'+ sample + '.sorted'
        covOut = tmpdir +'/'+ sample + '.out'
        with open(file, 'w') as tmpFasta:
            with open(bed, 'w') as bedFile:
                for gene in configDict['loci']:
                    genes = Fasta(configDict['loci'][gene])
                    allele = gene+'_'+re.sub("*", "", str(results[sample][gene]))
                    tmpFasta.write('>'+gene+'\n')
                    bedFile.write(gene+'\t0\t'+str(len(genes[allele]))+'\n')
                    for line in genes[allele]:
                        tmpFasta.write(str(line)+'\n')
        cmdIndex = "bwa index %s 2>/dev/null"%(file)
        os.system(cmdIndex)
        readBWA = sample+'_reads.fq'
        cmdBwaMem = "bwa mem %s %s 2>/dev/null| samtools view -uS - | samtools sort - -o %s"%(file, readBWA, sortedFile)
        os.system(cmdBwaMem)
        cmdCov = "bedtools coverage -a %s -b %s > %s"%(bed, sortedFile, covOut)
        os.system(cmdCov)
        with open(covOut, 'r') as cov:
            for line in cov.readlines():
                records = line.rstrip().rsplit('\t')
                gene = records[0]
                geneCov = float(records[6]) * 100
                results[sample][gene] = results[sample][gene] + " (" + str("%.2f" % geneCov) + ")"
    shutil.rmtree(tmpdir)
"""Prints the results in the format asked by the user."""
#############################################################
# Function   : printResults
# Input      : results, output file, overwrite?
# Output     : Prints on the screen or in a file
# Description: Prints the results in the format asked by the user
#############################################################
def printResults(results, output_filename, overwrite, timeDisp):
    if output_filename != None:
        if overwrite is False:
            outfile = open(output_filename, "a")
        else:
            outfile = open(output_filename, "w")
    heading = "Sample"
    for head in sorted(results[list(results.keys())[0]]):
        if head == 'ST' or head == 't':
            continue
        heading += '\t' + head
    heading += '\tST'
    if timeDisp is True:
        heading += '\tTime'
    if output_filename != None:
        outfile.write(heading)
        outfile.write('\n')
    else:
        print(heading)
    for s in results:
        sample = s.split("_")[0]
        for l in sorted(results[s]):
            if l == 'ST' or l == 't':
                continue
            sample += '\t'+results[s][l]
        if timeDisp is True:
            sample += '\t' + str(results[s]['ST']) + '\t%.2f ' %results[s]['t']
        else:
            sample += '\t' + str(results[s]['ST'])
        if output_filename != None:
            outfile.write(sample)
            outfile.write('\n')
        else:
            print(sample)
"""Predict part ends here"""
"""Build DB part starts"""
"""Returns the reverse complement of the sequence"""
def reverseComplement(seq):
    seqU = seq.upper()
    seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'N':'N'}
    try:
        return "".join([seq_dict[base] for base in reversed(seqU)])
    except Exception:
        strn = "Reverse Complement Error:" + seqU
        logging.debug(strn)
        pass
#############################################################
# Function   : getFastaDict
# Input      : locus file name
# Output     : dictionary with all the allele sequences
# Description: Stores each allele sequence in a dictionary
#############################################################
def getFastaDict(fullLocusFile):
    logging.debug("Create Fasta Dict")
    logging.debug(fullLocusFile)
    fastaFile = open(fullLocusFile, 'r').read()
    entries = [x for x in fastaFile.split('>') if len(x) != 0]
    fastaDict = {}
    for entry in entries:
        key = [x for x in entry.split('\n')[0].split() if len(x) != 0][0]
        sequence = ''.join(entry.split('\n')[1:]).rstrip()
        fastaDict[key] = {'sequence':sequence}
    return fastaDict
#############################################################
# Function   : formKmerDB
# Input      : configuration file, k value, output prefix
# Output     : stringMLST DB
# Description: Constructs the k-mer DB in both strand orientation
#############################################################
def formKmerDB(configDict, k, output_filename):
    dbFileName = output_filename+'_'+str(k)+'.txt'
    weightFileName = output_filename+'_weight.txt'
    kmerDict = {}
    mean = {}
    for locus in configDict['loci']:
        msgs = "formKmerDB :" +locus
        logging.debug(msgs)
        fastaDict = getFastaDict(configDict['loci'][locus])
        sum = 0
        n = 0
        for allele in list(fastaDict.keys()):
            seq = fastaDict[allele]['sequence'].strip()
            l = len(seq)
            sum += l
            n += 1
            try:
                (loc, num) = allele.split('_')
            except ValueError:
                print("Error : Allele name in locus file should be seperated by _")
                exit(0)
            splitId = allele.split('_')
            i = 0
            while i+k <= l:
                kmer = seq[i:i+k]
                revCompKmer = reverseComplement(kmer)
                if kmer not in kmerDict:
                    kmerDict[kmer] = {}
                    kmerDict[kmer][splitId[0]] = []
                    kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                else:
                    if splitId[0] not in kmerDict[kmer]:
                        kmerDict[kmer][splitId[0]] = []
                        kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                    else:
                        kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                if revCompKmer not in kmerDict:
                    kmerDict[revCompKmer] = {}
                    kmerDict[revCompKmer][splitId[0]] = []
                    kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                else:
                    if splitId[0] not in kmerDict[revCompKmer]:
                        kmerDict[revCompKmer][splitId[0]] = []
                        kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                    else:
                        kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                i += 1
        mean[locus] = sum/n*1.0
    with open(dbFileName, 'w') as kfile:
        for key in kmerDict:
            for key1 in kmerDict[key]:
                string = key+'\t'+key1+'\t'+str(kmerDict[key][key1]).replace(" ", "")+'\n'
                kfile.write(string)
    with open(weightFileName, 'w') as wfile:
        for locus in configDict['loci']:
            fastaDict = getFastaDict(configDict['loci'][locus])
            for allele in list(fastaDict.keys()):
                splitId = allele.split('_')
                seq = fastaDict[allele]['sequence']
                l = len(seq)
                fac = (l/mean[locus])
                s = allele  + '\t' + str(fac) + '\n'
                if fac > 1.05 or fac < 0.95:
                    wfile.write(s)
"""Copies the profile definition file as a new file"""
def copyProfileFile(profileDict, output_filename):
    profileFileName = output_filename+'_profile.txt'
    with open(profileDict['profile']) as f:
        lines = f.readlines()
        with open(profileFileName, "w") as f1:
            f1.writelines(lines)
#############################################################
# Function   : makeCustomDB
# Input      : configuration file, k value, output prefix
# Output     : None
# Description: Processes the config file and calls the relevant
#              function
#############################################################
def makeCustomDB(config, k, output_filename):
    configDict = {}
    if output_filename == None:
        output_filename = 'kmerDB'
    with open(config, 'r') as configFile:
        lines = configFile.readlines()
        head = ''
        for line in lines:
            if line.rstrip() == '':
                continue
            if line.rstrip() == '[loci]':
                head = 'loci'
                configDict[head] = {}
            elif line.rstrip() == '[profile]':
                head = 'profile'
                configDict[head] = {}
            else:
                arr = line.strip().split()
                configDict[head][arr[0]] = arr[1]
    for head in configDict:
        for element in configDict[head]:
            if not os.path.isfile(configDict[head][element]):
                print("ERROR: %s file does not exist at %s" % (element, configDict[head][element]))
                exit(0)
    formKmerDB(configDict, k, output_filename)
    copyProfileFile(configDict['profile'], output_filename)
"""Build DB part ends"""
"""Check Parameters"""
def checkParams(buildDB, predict, config, k, listMode, list, batch, dir, fastq1, fastq2, paired, dbPrefix):
    if predict is True and buildDB is True:
        print(helpTextSmall)
        print("Select either predict or buildDB module")
        exit(0)
    if predict is False and buildDB is False:
        print(helpTextSmall)
        print("Select either predict or buildDB module")
        exit(0)
    if predict is True:
        if config != None and coverage is False:
            print(helpTextSmall)
            print("Config parameter is not required for predict mode.")
            exit(0)
        elif config is None and coverage is True:
            print(helpTextSmall)
            print("Config parameter is required to for coverage prediction")
            exit(0)
        if not os.path.isfile(dbPrefix+'_'+str(k)+'.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_', str(k), '.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_weight.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_weight.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_profile.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_profile.txt or change DB prefix.')
            exit(0)
        if listMode is True:
            if not os.path.isfile(fList):
                print(helpTextSmall)
                print("Error: List file ("+fList+") does not exist!")
                exit(0)
        elif batch is True:
            if not os.path.isdir(dir):
                print(helpTextSmall)
                print("Error: Directory ("+dir+") does not exist!")
                exit(0)
        elif paired is True:
            if not os.path.isfile(fastq1):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
            if not os.path.isfile(fastq2):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq2+") does not exist!")
                exit(0)
        elif paired is False:
            if not os.path.isfile(fastq1):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
    if buildDB is True:
        try:
            if not os.path.isfile(config):
                print(helpTextSmall)
                print("Error: Configuration file ("+config+") does not exist!")
                exit(0)
        except Exception:
            print(helpTextSmall)
            print("Error: Specify Configuration file")
            exit(0)
helpText = """
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
            locus1      locusFile1
            locus2      locusFile2
            [profile]
            profile     profileFile
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
  Should have extention fastq or fq.
-2,--fastq2 = <fastq2_filename>
  Path to second fastq file for paired end sample.
  Should have extention fastq or fq.
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
=============================================================================================
Example usage:
./stringMLST.py --buildDB
1) Build DB
 ./stringMLST.py --buildDB --config config.txt -k 35 -P NM
 --------------------------------------------------------------------------------------------
./stringMLST.py --predict
1) Single sample, paired end
 ./stringMLST.py --predict -1 data/Neisseria/ERR017001_1.fastq -2 data/Neisseria/ERR017001_2.fastq -p --prefix NM -k 35 -o output.txt
2) Single sample, single end, overwrite output
  ./stringMLST.py --predict -1 data/Neisseria/ERR017001_1.fastq -s --prefix NM -k 35 -o output.txt -x
3) Multiple sample batch mode, paired end
   ./stringMLST.py --predict -d data/Neisseria/ -p --prefix NM -k 35 -o output.txt -x
4) Multiple samples list mode, paired end
   ./stringMLST.py --predict -l data/listFile.txt -p --prefix NM -k 35 -o output.txt -x
5) Single, high coverage sample, paired end
 ./stringMLST.py --predict -1 data/Neisseria/ERR017001_1.fastq -2 data/Neisseria/ERR017001_2.fastq -p --prefix NM -k 35 -z 1000 -o output.txt
"""
helpTextSmall = """
Usage
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
-h,--help
  Prints the help manual for this application
==============================================================================================
2. stringMLST.py --predict
Synopsis
stringMLST.py --predict -1 <fastq file> -2 <fastq file> -d <directory location> -l <list file> -p -s -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x
Required arguments
--predict
    Identifier for predict miodule
Optional arguments
-1,--fastq1 = <fastq1_filename>
  Path to first fastq file for paired end sample and path to the fastq file for single end file.
  Should have extention fastq or fq.
-2,--fastq2 = <fastq2_filename>
  Path to second fastq file for paired end sample.
  Should have extention fastq or fq.
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
=============================================================================================
"""
"""The Program Starts Execution Here"""
"""Default Params"""
buildDB = False
predict = False
output_filename = None
batch = False
listMode = False
overwrite = False
paired = True
fastq1 = None
fastq2 = None
user_k = False
config = None
timeDisp = False
reads = False
dbPrefix = 'kmer'
log = ''
k = 35
fuzzy = 300
coverage = False
#print 'ARGV      :', sys.argv[1:]
#exit(0)
"""Input arguments"""
options, remainder = getopt.getopt(sys.argv[1:], 'o:x1:2:k:l:bd:pshP:c:trva:z:C', [
    'buildDB',
    'predict',
    'output=',
    'config=',
    'prefix=',
    'overwrite',
    'batch',
    'list',
    'fastq1=',
    'fastq2=',
    'dir=',
    'directory=',
    'paired',
    'single',
    'help',
    'fuzzy=',
    'coverage'])
for opt, arg in options:
    if opt in ('-o', '--output'):
        output_filename = arg
    elif opt in ('-x', '--overwrite'):
        overwrite = True
    elif opt in '--buildDB':
        buildDB = True
    elif opt in ('-P', '--prefix'):
        dbPrefix = arg
    elif opt in '--predict':
        predict = True
    elif opt in ('-c', '--config'):
        config = arg
    elif opt in '-k':
        user_k = True
        try:
            k = int(arg)
        except ValueError:
            print("Error: Enter a numerical k value.")
            exit(0)
        # Check to make sure the arg is an int.
    elif opt in ('-l', '--list'):
        listMode = True
        fList = arg
    elif opt in ('-1', '--fastq1'):
        fastq1 = arg
    elif opt in ('-2', '--fastq2'):
        fastq2 = arg
    elif opt in ('-d', '--dir', '--directory'):
        dir = arg
        batch = True
    elif opt in ('-p', '--paired'):
        paired = True
        single = False
    elif opt in ('-s', '--single'):
        single = True
        paired = False
    elif opt in '-t':
        timeDisp = True
    elif opt in '-a':
        log = arg
    elif opt in '-r':
        reads = True
    elif opt in '-v':
        print(v)
        exit(0)
    elif opt in ('-C', '--coverage'):
        coverage = True
        reads = True
        fuzzy = 1
    elif opt in ('-z', '--fuzzy'):
        try:
            fuzzy = int(arg)
        except ValueError:
            print("You provided '" + arg + "' for your fuzziness threshold, which is not an integer value")
            exit(0)
    elif opt in ('-h', '--help'):
        print(helpText)
        exit(0)
checkParams(buildDB, predict, config, k, listMode, list, batch, dir, fastq1, fastq2, paired, dbPrefix)
if buildDB is True:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    if os.path.isfile(config):
        print("Info: Making DB for k = ", k)
        print("Info: Making DB with prefix =", dbPrefix)
        print("Info: Log file written to ", log)
        makeCustomDB(config, k, dbPrefix)
    else:
        print("Error: The input config file "+config +" does not exist.")
elif predict is True:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    loadModule(k, dbPrefix)
    if batch is True:
        results = batchTool(dir, paired, k)
    elif listMode is True:
        results = listTool(fList, paired, k)
    else:
        results = {}
        results = singleSampleTool(fastq1, fastq2, paired, k, results)
    if coverage is True:
        try:
            from pyfaidx import Fasta
        except ImportError:
            print("pyfaidx is required for coverage calculation\npip install pyfaidx")
            exit(0)
        loadConfig(config)
        getCoverage(results)
    printResults(results, output_filename, overwrite, timeDisp)
else:
    print(helpTextSmall)
    print("Error: Please select the mode: buildDB (for database building) or predict (for ST discovery) module")
logging.debug('Command :' + str(sys.argv))
