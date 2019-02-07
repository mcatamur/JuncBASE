# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 22:54:28 2017

@author: Carmelle
"""
#!/usr/bin/env python3

import csv, time, argparse
start_time = time.time()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
import numpy as np
import math

class CommandLine() :
    '''
    Takes input file and options.
    '''
    def __init__(self, inOpts=None) :
            '''
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            '''
            self.parser = argparse.ArgumentParser(description = 'outputAnalyzer.py arguments',
                                                 epilog = 'For more help contact: ',
                                                 add_help = True, #default is True
                                                 prefix_chars = '-',
                                                 usage = '%(prog)s [options] -option1[default]'
                                                 )
            self.parser.add_argument('-p', '--pValue', type=float, action = 'store', help='desired p value', default=0.01)
            self.parser.add_argument('-i', '--inputFile', action = 'store', help = 'input file')
            self.parser.add_argument('-sE', '--specificEvent', action = 'store', nargs='?', const=True, default=False, help='Event specific scatter plot')
            self.parser.add_argument('-e','--eventType', type=str, action = 'store', help = 'Select an alternative splicing event type. \n Omit space; replace with underscore.', default = 'cassette')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.pVal = myCommandLine.args.pValue
        self.scatterChoice = myCommandLine.args.specificEvent
        self.eventChoice = myCommandLine.args.eventType
        self.buildData(self.inputFile)
        self.plotEventsBar()

        if self.scatterChoice == True:
            self.verifyScatter(self.eventChoice)
        #else:
        #    self.plotAllEventsScatter()


    def buildData(self, inputFile):
        """
        Build dictionaries in this method for stacked bar plot and volcano scatter plot
        """
        self.allEvents = {}
        self.all_pVal = []
        self.filtered_pVal = []

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)
        #first we have to create the keys in the dictionary. we also have to take into account if the count is 0.
            for row in tsvReader:
                self.allEvents[row[1]] = 0

        self.filteredEvents = {key:0 for key,value in self.allEvents.items()}
        self.oneEvent_all = {key:[] for key,value in self.allEvents.items()}
        self.oneEvent_filtered = {key:[] for key,value in self.allEvents.items()}

        #Now, we have the dict keys we can create a dictionary of alternative splicing events for each
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)
            try:
                for row in tsvReader:
                    #all events
                    self.allEvents[row[1]] += 1

                    #filter events with user input p-value
                    if float(row[24]) < self.pVal:
                        self.filteredEvents[row[1]] += 1
                        self.filtered_pVal.append([float(row[22]), float(row[24])])
                        self.oneEvent_filtered[row[1]].append([float(row[22]), float(row[24])])
                    else:
                        self.all_pVal.append([float(row[22]), float(row[24])])
                        self.oneEvent_all[row[1]].append([float(row[22]), float(row[24])])

            except IndexError:
                pass

            except ValueError:
                pass


            #print (self.all_pVal)
            #print (self.allEvents)
            #print (self.filtered_pVal)
            #print (self.oneEvent_filtered)
            #print (self.oneEvent_all)


        #Will be used for the labels later and the modification we will be doing next
        self.totalEvents = 0
        self.totalAdjusted = 0

        for event in self.allEvents.keys():
            self.totalEvents += self.allEvents[event]
            self.totalAdjusted += self.filteredEvents[event]

        #print (self.filteredEvents)
        #print (self.totalAdjusted)
        #print (self.totalEvents)

        self.xlabel = ['All events. \n %s total events' % (self.totalEvents),
                   'Corrected p-value \n threshold p < %s. \n %s total events' % (self.pVal, self.totalAdjusted)]

        self.colors = ['#FFD54F', '#FFF176', '#81C784', '#4DD0E1', '#64B5F6', '#9575CD', '#F06292', '#EF5350', '#FF8A65', '#A1887F','#90A4AE']


        self.legend = ['All other events. \n %s total events' % (self.totalEvents - self.totalAdjusted),
                   'Corrected p-value \n threshold p < %s. \n %s total events' % (self.pVal, self.totalAdjusted)]



    def plotEventsBar(self):
        """
        This method creates a bar graph comparison between all significant alternative splicing events
        and events with a corrected p-value greater than the user defined input.
        """
        #We will do a little modification with the values for our bar plot
        #We will display the values as a fraction of the whole

        for event in self.allEvents.keys():
            self.allEvents[event] = self.allEvents[event] / self.totalEvents
            self.filteredEvents[event] = self.filteredEvents[event] / self.totalAdjusted

        #Multiply by 100, so we can have a percentage bar
        for event in self.allEvents.keys():
            self.allEvents[event] = self.allEvents[event] *100
            self.filteredEvents[event] = self.filteredEvents[event] *100

        #plot the results in a stacked bar graph


        #Turn the dict to a tuple. That way it is ordered and is subscriptable.
        a_events = list(self.allEvents.items())
        f_events = list(self.filteredEvents.items())
        x_vals = [0.5, 2.5]

        #Plot the Top bar first
        plt.bar(x_vals[0], a_events[0][1], bottom = 0, color = self.colors[0], label = a_events[0][0])
        plt.bar(x_vals[1], f_events[0][1], bottom = 0, color = self.colors[0])
        a_height = [a_events[0][1]]
        f_height = [f_events[0][1]]

        #Plot the rest
        for x in range(1, 20):
            try:
                plt.bar(x_vals[0], height = a_events[x][1], bottom=sum(a_height), color=self.colors[x], label=a_events[x][0])
                a_height.append(a_events[x][1])
                plt.bar(x_vals[1], height = f_events[x][1], bottom=sum(f_height), color=self.colors[x])
                f_height.append(f_events[x][1])
            except IndexError:
                continue


        plt.xticks(x_vals, self.xlabel)
        plt.ylabel('Proportion')
        plt.legend()
        plt.axis([0, 5.5, 0, 101])
        plt.show()
        plt.savefig('bar.jpeg')


    def verifyScatter(self, event):
        event = event.lower()

        for key in self.oneEvent_all.keys():
            if key == event:
                thisEvent = key
                break
            else:
                thisEvent = 'cassette'

        if thisEvent == 'cassette':
            print ('Default event selected: \'cassette\'')

        print (thisEvent)

        self.plotOneEventScatter(thisEvent)

    def plotOneEventScatter(self, event):
        y_all = []
        x_all = []

        for valueSet in self.oneEvent_all[event]:
            x_all.append(valueSet[0])
            y_all.append(-(math.log(valueSet[1])))

        #print (self.oneEvent_filtered)

        y_filtered = []
        x_filtered = []

        for valueSet in self.oneEvent_filtered[event]:
            x_filtered.append(valueSet[0])
            y_filtered.append(-(math.log(valueSet[1])))

        legend = ['All other \'%s\' events.\n %s total %s events' % (event, len(x_all), event),
                  '\'%s\' events with \n corrected p-value above \n the threshold. p < %s. \n %s total events' % (event, self.pVal, len(x_filtered))]

        plt.scatter(x_filtered,y_filtered, label = legend[1], color = self.colors[5], marker = '.', s = 60, alpha = 0.6)
        plt.scatter(x_all,y_all, label = legend[0], color = self.colors[1], marker = '.', s = 60, alpha = 0.5)

        plt.autoscale(enable=True, axis='both', tight=None)
        plt.ylabel('-log(corrected p-value)')
        plt.xlabel('Delta PSI Value')
        plt.legend()



    def plotAllEventsScatter(self):

        y_all = []
        x_all = []

        for value in self.all_pVal:
            x_all.append(value[0])
            y_all.append(-(math.log(value[1]))) #Modify the corrected p-values to -log(corrected p-value)

        y_filtered = []
        x_filtered = []

        #print (self.filtered_pVal)

        for value in self.filtered_pVal:
            x_filtered.append(value[0])
            y_filtered.append(-(math.log(value[1])))


        plt.scatter(x_filtered,y_filtered, label = self.legend[1], color = self.colors[5], marker = '.', s = 60, alpha = 0.6)
        plt.scatter(x_all,y_all, label = self.legend[0], color = self.colors[1], marker = '.', s = 60, alpha = 0.5)

        plt.autoscale(enable=True, axis='both', tight=None)
        #plt.axis([-50, 50, 0, 20])
        plt.ylabel('-log(corrected p-value)')
        plt.xlabel('Delta PSI Value')
        plt.legend()

        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
