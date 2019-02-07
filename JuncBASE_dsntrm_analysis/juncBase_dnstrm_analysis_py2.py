#!/usr/bin/env python3
# Names: Daniel Schmelter(dschmelt) and Carmelle Catmura(mcatmura)
# Group Members: Carmelle Catmura(mcatmura) and Daniel Schmelter(dschmelt)
# Last Updated: 6/12/17

"""
JuncBase_outputAnalyzer.py is an easy to use tool that finds and plots significant
alternative splicing events correlated with mutations in a particular gene

This program takes JuncBase's output file and user options from command line
It counts and returns significant events by pValue threshold
It makes parallel violin plots of mutant and wild type PSI distributions with your selected gene.
It makes a bar graph counting different types of splicing events in your JuncBase output file.
It also prints runtime of the program. ~90 seconds

Example usage:
python3 JuncBase_outputAnalyzer.py -i KRAS_JuncBase_Output.tsv -p 0.01 -g COPA -m KRAS_Patients.txt -wt NoKRAS_Patients.txt
"""
from collections import defaultdict
import csv, sys, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

class CommandLine() :
    '''
    Takes input file and options.
    '''
    def __init__(self, inOpts=None) :
            '''
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            '''
            self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does',
                                                 epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                                 add_help = True, #default is True 
                                                 prefix_chars = '-', 
                                                 usage = '%(prog)s [options] -option1[default] <input >output'
                                                 )
            self.parser.add_argument('-p', '--pValue', type=float, action = 'store', help='desired p value', default=0.05)
            self.parser.add_argument('-i', '--inputFile', action = 'store', help = 'input file')
            self.parser.add_argument('-m', '--MutationPatients', action = 'store', help = 'input file')
            self.parser.add_argument('-wt', '--WildtypePatients', action = 'store', help = 'input file')
            self.parser.add_argument('-g','--inputGene', action = 'store', help='Select a gene to show violin plot', default='COPA')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)


class tsvReader:
    def __init__(self, myCommandline):
        self.pVal = myCommandline.args.pValue
        self.inputFile = myCommandline.args.inputFile
        self.inputGene = myCommandline.args.inputGene
        self.allPSI_List = []
        self.mutPatients = myCommandline.args.MutationPatients
        self.wtPatients = myCommandline.args.WildtypePatients
        self.wtPSIList = []
        self.mutPSIList = []
        self.buildDict(self.inputFile)
        #self.splicingEventBar(self.inputFile)
        #self.plotBar(self.inputFile)
        self.plotViolin()

    def buildDict(self, tsv_file):
        self.allPValDict = defaultdict(list)
        self.pValDictionary = defaultdict(list)
        self.psiDictionary = defaultdict(list)
        try:
            with open(tsv_file) as file:
                reader = csv.reader(file, dialect = 'excel-tab')
                next(reader, None)  #skips column headers
                for row in reader:
                    try:
                        self.allPValDict[row[2]].append(row[510])
                        if float(row[510]) < self.pVal: #corrected pValue column
                            self.pValDictionary[row[2]].append(row[510])

                            for columnNumber in range(14,506): #all the PSI values
                                self.psiDictionary[row[2]].append(row[columnNumber])
                                if self.inputGene in row:
                                    try:
                                        tempPSI = float(row[columnNumber])
                                        self.allPSI_List.append(tempPSI)
                                    except ValueError:
                                        continue
                    except IndexError:
                        pass
            print("\nYou have {} significant alternatively spliced genes at {} p-value."
                  .format(len(self.pValDictionary),self.pVal))
            print("Here is a dictionary of genes below your {} threshold and p-values for each event:".format(self.pVal))
            #print(dict(self.pValDictionary.items()))
        except TypeError:
            print("Usage: python3 JuncBase_outputAnalyzer.py -i <inputFile.tsv> [-p <P value>]")
            sys.exit()
        self.getMutNonmutPSI()


    def getMutNonmutPSI(self):
        mutIDList = []
        wtIDList = []
        with open(self.mutPatients) as mutFile:
            reader = csv.reader(mutFile, dialect = 'excel-tab')
            for mutRow in mutFile:
                mutIDList.append(mutRow[: -1])
        with open(self.wtPatients) as wtFile:
            reader2 = csv.reader(wtFile,dialect='excel-tab')
            for wtRow in wtFile:
                wtIDList.append(wtRow[: -1])
        #print("here's the ID list" + str(mutIDList) + str(wtIDList))


        mutIndexList = []
        wtIndexList = []
        initialized = False
        rowList = []
        temp = True

        with open(self.inputFile) as infile:
            reader3 = csv.reader(infile, dialect='excel-tab')
            for row in infile:
                """
                Takes infile, returns Index list of mut and wt patients
                """
                if temp:
                    #print(row)
                    temp=False
                rowList = row.split("\t")
                #print(rowList)
                if initialized==False:
                    initialized = True
                    for patient in wtIDList:
                        mutIndexList.append(rowList.index(patient))
                    for patient in mutIDList:
                        wtIndexList.append(rowList.index(patient))
                    #print(rowList)
                    continue
                """
                Takes index list and returns PSI list
                """
                if self.inputGene in rowList[2]:
                    try:
                        for wtPatientIndex in wtIndexList:     #loops through patientIndexList
                                if float(row[510]) < self.pVal:
                                    self.wtPSIList.append(float(rowList[wtPatientIndex]))
                    except ValueError:
                        pass
                    try:
                        for mutPatientIndex in mutIndexList:
                                self.mutPSIList.append(float(rowList[mutPatientIndex]))
                    except ValueError: pass
                else: continue
            #print(wtIndexList)

    def plotViolin(self):
        """
        plotViolin creates parallel violin o
        """
        #print("wtPSIs:\n")
        #print(self.wtPSIList)print(self.wtPSIList)
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
        allPSIs = self.allPSI_List
        mutPSIs = self.mutPSIList
        wtPSIs = self.wtPSIList
        axes[0].violinplot(wtPSIs, showmeans = False, showmedians = True)
        axes[1].violinplot(mutPSIs, showmeans = False, showmedians = True)
        axes[0].set_title("PSI distribution for "+ str(self.inputGene)+" with wt KRAS")
        axes[1].set_title("PSI distribution for "+str(self.inputGene)+" with mutant KRAS")
        axes[0].set_ylabel('Percent Spliced In (PSI)')
        plt.suptitle("{} p-value is: {}".format(str(self.inputGene),
                        self.allPValDict[self.inputGene]), fontsize=16)
        #legendary = mpatches.Patch(label="{} p-value is: {}".format(str(self.inputGene), self.allPValDict[self.inputGene][-1]))
        #plt.legend(handles=[legendary])
        plt.show()

    def plotBar(self, tsv_file):
        """
        plotBar takes in the tsv file output by JuncBase
        It counts the splicingEvent type from that file
        It then makes a bar Graph of those type counts
        """

        splicingEvent = {}
        with open(tsv_file) as file:
            reader = csv.reader(file, dialect='excel-tab')
            next(reader, None)
            for row in reader:
                splicingEvent[row[1]] = 0
        with open(tsv_file) as file:
            reader = csv.reader(file, dialect='excel-tab')
            next(reader, None)
            for row in reader:
                splicingEvent[row[1]] += 1
        plt.title('KRAS associated Splicing Event Types', fontsize=16)
        x = splicingEvent.keys()
        y = splicingEvent.values()
        y_pos = np.arange(len(x))
        plt.bar(y_pos, y)
        plt.xticks(y_pos, x, rotation='vertical')
        plt.subplots_adjust(bottom=0.35)
        plt.ylabel("Occurrences")


    def splicingEventBar(self, inputFile):
        self.mtIndexList = []
        self.wtIndexList = []
        self.mtEventDict = {}
        self.wtEventDict = {}
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            fileHeader = next(tsvReader)
            #To access the values, we have to create separate index lists
            with open(self.mutPatients) as mtFile: #creates index list of mutated patients
                mtReader = csv.reader(mtFile, dialect = 'excel-tab')
                for patientID in mtReader:
                    for patient in patientID:
                        self.mtIndexList.append(fileHeader.index(patient))
            with open(self.wtPatients) as wtFile: #creates index list of wild type patients
                mtReader = csv.reader(wtFile, dialect = 'excel-tab')
                for patientID in mtReader:
                    for patient in patientID:
                        self.wtIndexList.append(fileHeader.index(patient))
        #first we have to create the keys in the dictionary. we also have to take into account if the count is 0.
            for row in tsvReader:
                self.mtEventDict[row[1]] = 0
                self.wtEventDict[row[1]] = 0
        #Now, we have the index list, and the dict keys we can create a dictionary of alternative splicing events for each
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)
            for row in tsvReader:
                try:
                    if float(row[510]) < self.pVal:
                        for patient in self.mtIndexList:
                            try:
                                if float(row[patient]) > 0:
                                    self.mtEventDict[row[1]] += 1
                            except ValueError:
                                continue
                        for patient in self.wtIndexList:
                            try:
                                if float(row[patient]) > 0:
                                    self.mtEventDict[row[1]] += 1
                            except ValueError:
                                continue
                except IndexError:
                    continue
        print (self.wtEventDict)
        print (self.mtEventDict)
        

def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.
    '''

    myCommandLine = CommandLine(myCommandLine)
    myFileReader = tsvReader(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()


#scatterplot w/ KRAS median vs non median
#stacked barplot of KRAS and non-KRAS
#filter violin plots by event instead of all
#filter by p value and then call violin plot
#volcano plot