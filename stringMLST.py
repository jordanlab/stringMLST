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
import xml.etree.ElementTree as ET
try:
    from urllib.request import urlopen, urlretrieve
except ImportError:
    from urllib import urlopen, urlretrieve
import argparse
version = """ stringMLST v0.6.2 (updated : March, 5 2019) """
"""

stringMLST free for academic users and requires permission before any commercial
use for any version of this code/algorithm. If you are a commercial user, please
contact king.jordan@biology.gatech.edu for permissions

LICENSE TERMS FOR stringMLST
Adopted from: https://creativecommons.org/licenses/by-nc-sa/4.0/

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public
License

By exercising the Licensed Rights (defined below), You accept and agree to be
bound by the terms and conditions of this Creative Commons Attribution-
NonCommercial-ShareAlike 4.0 International Public License ("Public License"). To
the extent this Public License may be interpreted as a contract, You are granted
the Licensed Rights in consideration of Your acceptance of these terms and
conditions, and the Licensor grants You such rights in consideration of benefits
the Licensor receives from making the Licensed Material available under these
terms and conditions.

Section 1 - Definitions.

Adapted Material means material subject to Copyright and Similar Rights that is
derived from or based upon the Licensed Material and in which the Licensed
Material is translated, altered, arranged, transformed, or otherwise modified in
a manner requiring permission under the Copyright and Similar Rights held by the
Licensor. For purposes of this Public License, where the Licensed Material is a
musical work, performance, or sound recording, Adapted Material is always
produced where the Licensed Material is synched in timed relation with a moving
image. Adapter's License means the license You apply to Your Copyright and
Similar Rights in Your contributions to Adapted Material in accordance with the
terms and conditions of this Public License. BY-NC-SA Compatible License means a
license listed at creativecommons.org/compatiblelicenses, approved by Creative
Commons as essentially the equivalent of this Public License. Copyright and
Similar Rights means copyright and/or similar rights closely related to
copyright including, without limitation, performance, broadcast, sound
recording, and Sui Generis Database Rights, without regard to how the rights are
labeled or categorized. For purposes of this Public License, the rights
specified in Section 2(b)(1)-(2) are not Copyright and Similar Rights. Effective
Technological Measures means those measures that, in the absence of proper
authority, may not be circumvented under laws fulfilling obligations under
Article 11 of the WIPO Copyright Treaty adopted on December 20, 1996, and/or
similar international agreements. Exceptions and Limitations means fair use,
fair dealing, and/or any other exception or limitation to Copyright and Similar
Rights that applies to Your use of the Licensed Material. License Elements means
the license attributes listed in the name of a Creative Commons Public License.
The License Elements of this Public License are Attribution, NonCommercial, and
ShareAlike. Licensed Material means the artistic or literary work, database, or
other material to which the Licensor applied this Public License. Licensed
Rights means the rights granted to You subject to the terms and conditions of
this Public License, which are limited to all Copyright and Similar Rights that
apply to Your use of the Licensed Material and that the Licensor has authority
to license. Licensor means the individual(s) or entity(ies) granting rights
under this Public License. NonCommercial means not primarily intended for or
directed towards commercial advantage or monetary compensation. For purposes of
this Public License, the exchange of the Licensed Material for other material
subject to Copyright and Similar Rights by digital file-sharing or similar means
is NonCommercial provided there is no payment of monetary compensation in
connection with the exchange. Share means to provide material to the public by
any means or process that requires permission under the Licensed Rights, such as
reproduction, public display, public performance, distribution, dissemination,
communication, or importation, and to make material available to the public
including in ways that members of the public may access the material from a
place and at a time individually chosen by them. Sui Generis Database Rights
means rights other than copyright resulting from Directive 96/9/EC of the
European Parliament and of the Council of 11 March 1996 on the legal protection
of databases, as amended and/or succeeded, as well as other essentially
equivalent rights anywhere in the world. You means the individual or entity
exercising the Licensed Rights under this Public License. Your has a
corresponding meaning. Section 2 - Scope.

License grant. Subject to the terms and conditions of this Public License, the
Licensor hereby grants You a worldwide, royalty-free, non-sublicensable, non-
exclusive, irrevocable license to exercise the Licensed Rights in the Licensed
Material to: reproduce and Share the Licensed Material, in whole or in part, for
NonCommercial purposes only; and produce, reproduce, and Share Adapted Material
for NonCommercial purposes only. Exceptions and Limitations. For the avoidance
of doubt, where Exceptions and Limitations apply to Your use, this Public
License does not apply, and You do not need to comply with its terms and
conditions. Term. The term of this Public License is specified in Section 6(a).
Media and formats; technical modifications allowed. The Licensor authorizes You
to exercise the Licensed Rights in all media and formats whether now known or
hereafter created, and to make technical modifications necessary to do so. The
Licensor waives and/or agrees not to assert any right or authority to forbid You
from making technical modifications necessary to exercise the Licensed Rights,
including technical modifications necessary to circumvent Effective
Technological Measures. For purposes of this Public License, simply making
modifications authorized by this Section 2(a)(4) never produces Adapted
Material. Downstream recipients. Offer from the Licensor - Licensed Material.
Every recipient of the Licensed Material automatically receives an offer from
the Licensor to exercise the Licensed Rights under the terms and conditions of
this Public License. Additional offer from the Licensor - Adapted Material.
Every recipient of Adapted Material from You automatically receives an offer
from the Licensor to exercise the Licensed Rights in the Adapted Material under
the conditions of the Adapter's License You apply. No downstream restrictions.
You may not offer or impose any additional or different terms or conditions on,
or apply any Effective Technological Measures to, the Licensed Material if doing
so restricts exercise of the Licensed Rights by any recipient of the Licensed
Material. No endorsement. Nothing in this Public License constitutes or may be
construed as permission to assert or imply that You are, or that Your use of the
Licensed Material is, connected with, or sponsored, endorsed, or granted
official status by, the Licensor or others designated to receive attribution as
provided in Section 3(a)(1)(A)(i). Other rights.

Moral rights, such as the right of integrity, are not licensed under this Public
License, nor are publicity, privacy, and/or other similar personality rights;
however, to the extent possible, the Licensor waives and/or agrees not to assert
any such rights held by the Licensor to the limited extent necessary to allow
You to exercise the Licensed Rights, but not otherwise. Patent and trademark
rights are not licensed under this Public License. To the extent possible, the
Licensor waives any right to collect royalties from You for the exercise of the
Licensed Rights, whether directly or through a collecting society under any
voluntary or waivable statutory or compulsory licensing scheme. In all other
cases the Licensor expressly reserves any right to collect such royalties,
including when the Licensed Material is used other than for NonCommercial
purposes. Section 3 - License Conditions.

Your exercise of the Licensed Rights is expressly made subject to the following
conditions.

Attribution.

If You Share the Licensed Material (including in modified form), You must:

retain the following if it is supplied by the Licensor with the Licensed
Material: identification of the creator(s) of the Licensed Material and any
others designated to receive attribution, in any reasonable manner requested by
the Licensor (including by pseudonym if designated); a copyright notice; a
notice that refers to this Public License; a notice that refers to the
disclaimer of warranties; a URI or hyperlink to the Licensed Material to the
extent reasonably practicable; indicate if You modified the Licensed Material
and retain an indication of any previous modifications; and indicate the
Licensed Material is licensed under this Public License, and include the text
of, or the URI or hyperlink to, this Public License. You may satisfy the
conditions in Section 3(a)(1) in any reasonable manner based on the medium,
means, and context in which You Share the Licensed Material. For example, it may
be reasonable to satisfy the conditions by providing a URI or hyperlink to a
resource that includes the required information. If requested by the Licensor,
You must remove any of the information required by Section 3(a)(1)(A) to the
extent reasonably practicable. ShareAlike. In addition to the conditions in
Section 3(a), if You Share Adapted Material You produce, the following
conditions also apply.

The Adapter's License You apply must be a Creative Commons license with the same
License Elements, this version or later, or a BY-NC-SA Compatible License. You
must include the text of, or the URI or hyperlink to, the Adapter's License You
apply. You may satisfy this condition in any reasonable manner based on the
medium, means, and context in which You Share Adapted Material. You may not
offer or impose any additional or different terms or conditions on, or apply any
Effective Technological Measures to, Adapted Material that restrict exercise of
the rights granted under the Adapter's License You apply. Section 4 - Sui
Generis Database Rights.

Where the Licensed Rights include Sui Generis Database Rights that apply to Your
use of the Licensed Material:

for the avoidance of doubt, Section 2(a)(1) grants You the right to extract,
reuse, reproduce, and Share all or a substantial portion of the contents of the
database for NonCommercial purposes only; if You include all or a substantial
portion of the database contents in a database in which You have Sui Generis
Database Rights, then the database in which You have Sui Generis Database Rights
(but not its individual contents) is Adapted Material, including for purposes of
Section 3(b); and You must comply with the conditions in Section 3(a) if You
Share all or a substantial portion of the contents of the database. For the
avoidance of doubt, this Section 4 supplements and does not replace Your
obligations under this Public License where the Licensed Rights include other
Copyright and Similar Rights. Section 5 - Disclaimer of Warranties and
Limitation of Liability.

Unless otherwise separately undertaken by the Licensor, to the extent possible,
the Licensor offers the Licensed Material as-is and as-available, and makes no
representations or warranties of any kind concerning the Licensed Material,
whether express, implied, statutory, or other. This includes, without
limitation, warranties of title, merchantability, fitness for a particular
purpose, non-infringement, absence of latent or other defects, accuracy, or the
presence or absence of errors, whether or not known or discoverable. Where
disclaimers of warranties are not allowed in full or in part, this disclaimer
may not apply to You. To the extent possible, in no event will the Licensor be
liable to You on any legal theory (including, without limitation, negligence) or
otherwise for any direct, special, indirect, incidental, consequential,
punitive, exemplary, or other losses, costs, expenses, or damages arising out of
this Public License or use of the Licensed Material, even if the Licensor has
been advised of the possibility of such losses, costs, expenses, or damages.
Where a limitation of liability is not allowed in full or in part, this
limitation may not apply to You. The disclaimer of warranties and limitation of
liability provided above shall be interpreted in a manner that, to the extent
possible, most closely approximates an absolute disclaimer and waiver of all
liability. Section 6 - Term and Termination.

This Public License applies for the term of the Copyright and Similar Rights
licensed here. However, if You fail to comply with this Public License, then
Your rights under this Public License terminate automatically. Where Your right
to use the Licensed Material has terminated under Section 6(a), it reinstates:

automatically as of the date the violation is cured, provided it is cured within
30 days of Your discovery of the violation; or upon express reinstatement by the
Licensor. For the avoidance of doubt, this Section 6(b) does not affect any
right the Licensor may have to seek remedies for Your violations of this Public
License. For the avoidance of doubt, the Licensor may also offer the Licensed
Material under separate terms or conditions or stop distributing the Licensed
Material at any time; however, doing so will not terminate this Public License.
Sections 1, 5, 6, 7, and 8 survive termination of this Public License. Section 7
- Other Terms and Conditions.

The Licensor shall not be bound by any additional or different terms or
conditions communicated by You unless expressly agreed. Any arrangements,
understandings, or agreements regarding the Licensed Material not stated herein
are separate from and independent of the terms and conditions of this Public
License. Section 8 - Interpretation.

For the avoidance of doubt, this Public License does not, and shall not be
interpreted to, reduce, limit, restrict, or impose conditions on any use of the
Licensed Material that could lawfully be made without permission under this
Public License. To the extent possible, if any provision of this Public License
is deemed unenforceable, it shall be automatically reformed to the minimum
extent necessary to make it enforceable. If the provision cannot be reformed, it
shall be severed from this Public License without affecting the enforceability
of the remaining terms and conditions. No term or condition of this Public
License will be waived and no failure to comply consented to unless expressly
agreed to by the Licensor. Nothing in this Public License constitutes or may be
interpreted as a limitation upon, or waiver of, any privileges and immunities
that apply to the Licensor or You, including from the legal processes of any
jurisdiction or authority.



The program has 3 basic modes :
    mainTool: for single sample (both single and paired end)
    batchTool: for multiple samples stored at a common location (both single and paired end samples)
    listTool: for multiple samples with location information stored in a list (both single and paired end samples)
predict part starts here
"""
#############################################################
# Function   : get_links
# Input      : speciesName and schemes dict
# Output     : Dict containing links to alleles and profile
# Description: Gets the URLs from pubMLST for the required
#              files (alleles, profile)
#############################################################
def get_links(xmlData, savePath, speciesName):
    lociList = {}
    profileURL = None
    for species in xmlData:
        if re.search(re.escape(speciesName), species.text, re.IGNORECASE, ):
            for mlst in species:
                for database in mlst:
                    for child in database:
                        if child.tag == "profiles":
                            profileURL = child[1].text
                        if child.tag == "loci":
                            for locus in child:
                                lociList[locus.text.rstrip()] = locus[0].text
    if profileURL is None:
        profileError = "Parsing failed: could not find profiles file"
        print(profileError)
        print("This usually means the provided species, '{}', does not exist on PubMLST".format(speciesName))
        print("Use `{} --getMLST --species list` to list available species".format(sys.argv[0]))
        print("Or visit PubMLST for more information:\nhttps://pubmlst.org/data/")
        logging.debug(profileError)
        sys.exit(1)
    elif lociList == {}:
        lociError = "Parsing failed: could not find allele sequences"
        logging.debug(lociError)
        print(lociError)
        sys.exit(1)
    else:
        return profileURL, lociList
#############################################################
# Function   : get_files
# Input      : URLs from get_links
# Output     : Downloads files and builds database
#############################################################
def get_files(filePrefix, loci, profileURL, speciesName):
    with open(config, "w") as configFile:
        configFile.write("[loci]\n")
        for file in loci:
            localFile = filePrefix + "_" + file + ".tfa"
            try:
                localFile, headers = urlretrieve(loci[file], localFile)
            except:
                print('\033[91m' + "There was an error downloading " + file + '\033[0m')
                pass
            configFile.write(file + "\t" + filePrefix + "_" + file + ".tfa\n")
        localFile = filePrefix + "_profile.txt"
        localFile, headers = urlretrieve(profileURL, localFile)
        configFile.write("[profile]\n")
        configFile.write("profile\t" + filePrefix + "_profile.txt\n")
        configFile.close()
        try:
            makeCustomDB(config, k, filePrefix)
        except:
            print('\033[91m' + "Failed to create database " + speciesName + '\033[0m')
            pass
        else:
            print("\t" + '\033[92m' + "Database ready for " + speciesName + '\033[0m')
            print("\t" + filePrefix)
############################################################
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
        if sys.version_info[0] == 3:
            f = gzip.open(file, 'rt')
        else:
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
    for locus in alleleNames:
        finalProfileCount[locus] = {}
    num = ''
    for loc in alleleCount:
        n = 0
        m = 0
        for num in alleleCount[loc]:
            if alleleCount[loc][num] >= n:
                m = n
                n = alleleCount[loc][num]
        if n-m < fuzzy:
            try:
                alleleCount[loc][num]
            except:
                pass
            else:
                alleleCount[loc][num] = str(alleleCount[loc][num])+'*'
                max_n[loc] = str(n)+'*'
        else:
            max_n[loc] = n
        secondMax[loc] = m
    for loc in alleleCount:
        try:
            max_n[loc]
        except:
            pass
        else:
            maxSupport[loc] = {}
            secondSupport[loc] = {}
            num_max = []
            num_max2 = []
            compare = float(re.sub("\*$", "", str(max_n[loc])))
            for num in alleleCount[loc]:
                if  float(re.sub("\*$", "", str(alleleCount[loc][num]))) == compare:
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
                finalProfileCount[loc] = 'NA'
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
        if allele:
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
        global alleleNames
        alleleNames = set()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            kmerTableDict[array[0]] = {}
            kmerTableDict[array[0]][array[1]] = array[2][1:-1].rsplit(',')
            alleleNames.add(array[1])
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
            try:
               (loc, allele) = array[0].replace('-', '_').rsplit('_', 1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0)
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
                    allele = gene+'_'+re.sub('\*', "", str(results[sample][gene]))
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
            if results[s][l]:
                sample += '\t'+results[s][l]
            else:
                sample += '\tNA'
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
                (loc, num) = allele.replace('-', '_').rsplit('_', 1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0)
            splitId = allele.replace('-', '_').rsplit('_', 1)
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
    if predict is False and buildDB is False and downloadDB is False:
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
    Identifier for predict module
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
    Calculate sequence coverage for each allele. Turns on read generation (-r) and turns off fuzzy (-z 1)
    Requires bwa, bamtools and samtools be in your path
-k = <kmer_length>
  Kmer length for which the db was created(Default k = 35). Could be verified by looking at the name of the db file.
  Could be used if the reads are of very bad quality or have a lot of N's.
-l,--list = <list_file>
  LIST MODE : Location of list file and flag for list mode.
  list file should have full file paths for all the samples/files.
  Each sample takes one line. For paired end samples the 2 files should be tab separated on single line.
-o,--output = <output_filename>
  Prints the output to a file instead of stdout.
-p,--paired
  Flag for specifying paired end files. Default option so would work the same if you do not specify for all modes.
  For batch mode the paired end samples should be differentiated by 1/2.fastq or 1/2.fq
-P,--prefix = <prefix>
    Prefix using which the db was created(Defaults = kmer). The location of the db could also be provided.
-r
  A separate reads file is created which has all the reads covering all the locus.
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
--------------------------------------------------------------------------------------------
./stringMLST.py --getMLST
1) List available schemes
 ./stringMLST.py --getMLST --schemes
2) Download the Neisseria spp. pubMLST scheme
  ./stringMLST.py --getMLST --species=neisseria -P datasets/nmb
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
    Prefix for db and log files to be created(Default = kmer). Also you can specify folder where you want the db to be created.
-h,--help
  Prints the help manual for this application
==============================================================================================
2. stringMLST.py --predict
Synopsis
stringMLST.py --predict -1 <fastq file> -2 <fastq file> -d <directory location> -l <list file> -p -s -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x
Required arguments
--predict
    Identifier for predict module
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
    Calculate sequence coverage for each allele. Turns on read generation (-r) and turns off fuzzy (-z 1)
    Requires bwa, bamtools and samtools be in your path
-k = <kmer_length>
  Kmer length for which the db was created(Default k = 35). Could be verified by looking at the name of the db file.
  Could be used if the reads are of very bad quality or have a lot of N's.
-l,--list = <list_file>
  LIST MODE : Location of list file and flag for list mode.
  list file should have full file paths for all the samples/files.
  Each sample takes one line. For paired end samples the 2 files should be tab separated on single line.
-o,--output = <output_filename>
  Prints the output to a file instead of stdout.
-p,--paired
  Flag for specifying paired end files. Default option so would work the same if you do not specify for all modes.
  For batch mode the paired end samples should be differentiated by 1/2.fastq or 1/2.fq
-P,--prefix = <prefix>
    Prefix using which the db was created(Defaults = kmer). The location of the db could also be provided.
-r
  A separate reads file is created which has all the reads covering all the locus.
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
3. stringMLST.py --getMLST
Synopsis:
stringMLST.py --getMLST --species= <species> [-k kmer length] [-P DB prefix]
Required arguments
--getMLST
    Identifier for getMLST module
--species= <species name>
    Species name from the pubMLST schemes
    Use "show" or "list" to list available schemes
    "all" will download and build all available schemes
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
=============================================================================================

"""

"""The Program Starts Execution Here"""
"""Default Params"""
downloadDB = False
species = None
printSchemes = False
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
#print'ARGV      :', sys.argv[1:]
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
    'coverage',
    'getMLST',
    'schemes',
    'species='])
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
        print(version)
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
    elif opt in '--schemes':
        print("The `--schemes` option has been depreciated.  Please use `--species list` to see available schemes")
        exit(0)
    elif opt in '--getMLST':
        downloadDB = True
    elif opt in '--species':
        species = arg
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
elif downloadDB is True:
    dbURL = "http://pubmlst.org/data/dbases.xml"
    databaseXML = urlopen(dbURL)
    dbTree = ET.parse(databaseXML)
    dbRoot = dbTree.getroot()
    if species is None:
        print("Please refer to --help to more information")
        print()
        print("Expected command format:")
        print("stringMLST.py --getMLST --species= <species> [-k kmer length] [-P DB prefix]")
        print()
        print("To printavailable MLST Schemes use:")
        print("stringMLST.py --getMLST --species show")
        exit(0)
    elif species == "show" or species == "list":
        for species in dbRoot:
            print(species.text.rstrip())
    elif species == "all":
        print("Using a kmer size of " + str(k) + " for all databases.")
        for species in dbRoot:
            speciesName = species.text.rstrip()
            print('\033[1m' + "Preparing: " + speciesName + '\033[0m')
            if re.search('[/#. ()]', speciesName):
                normSpeciesName = re.sub('[/# ]', "_", speciesName)
                normSpeciesName = re.sub('[.()]', "", normSpeciesName)
                print('\t\033[33m' + "INFO: normalizing name to: " + normSpeciesName + '\033[0m')
            else:
                normSpeciesName = speciesName
            filePrefix = str(dbPrefix.rsplit("/", 1)[0]) + "/" + normSpeciesName
            # Move the rest of this informational message into the download handler
            # + " ( " + filePrefix + "/" + key + "_" +str(k) + " )")
            try:
                os.makedirs(filePrefix)
            except OSError:
               pass
            filePrefix = filePrefix + "/" + normSpeciesName
            config = filePrefix + "_config.txt"
            profileURL, loci = get_links(dbRoot, filePrefix, speciesName)
            get_files(filePrefix, loci, profileURL, speciesName)
    else:
        print('\033[1m' + "Preparing: " + species + '\033[0m')
        if re.search('[/#. ()]', species):
            normSpeciesName = re.sub('[/# ]', "_", species)
            normSpeciesName = re.sub('[.()]', "", normSpeciesName)
            print('\t\033[33m' + "INFO: normalizing name to: " + normSpeciesName + '\033[0m')
        else:
            normSpeciesName = species
        try:
            os.makedirs(dbPrefix.rsplit("/", 1)[0])
        except OSError:
            pass
        if len(re.findall("/", dbPrefix)) == 0:
            filePrefix = dbPrefix + "/" + normSpeciesName
        elif len(re.findall("/", dbPrefix)) == 1 and len(dbPrefix.rsplit("/", 1)[1]) > 0:
            filePrefix = dbPrefix
        elif len(re.findall("/", dbPrefix)) == 1 and len(dbPrefix.rsplit("/", 1)[1]) == 0:
            filePrefix = dbPrefix + normSpeciesName
        elif len(re.findall("/", dbPrefix)) > 1:
            if dbPrefix.endswith('/'):
                filePrefix = dbPrefix + normSpeciesName
            else:
                filePrefix = dbPrefix
        config = filePrefix + "_config.txt"
        profileURL, loci = get_links(dbRoot, filePrefix, species)
        get_files(filePrefix, loci, profileURL, species)
else:
    print(helpTextSmall)
    print("Error: Please select the mode: buildDB (for database building) or predict (for ST discovery) module")
logging.debug('Command :' + str(sys.argv))

