#!/usr/bin/env python
"""
---------------------------------------------
Friend Matching for MGB

Author: {author}
Contact: {email}
---------------------------------------------
"""
import os
import string
import re
import sys
from itertools import *
from collections import defaultdict, Counter
import ConfigParser
import logging
import argparse

__author__ = 'Hyukjung Kwon'
#__company__ = "EDGC"
__email__ = "hjkwon@edgc.com"
__version__ = '0.9'

logger = logging.getLogger(__name__)
# Specifies format of log
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)

target_files = []
target_files_RS = []
target_files_GeneData = []
#Final_NormalizedData = {}
FinalDataSet = []

# Let's use constant value
# Modify later to get real size in data file
chr_size = {1:249250621, 2:243199373, 3:198022430, 4:191154276, 5:180915260, 6:171115067, 7:159138663, 8:146364022 \
            , 9:141213431, 10:135534747, 11:135006516, 12:133851895, 13:115169878, 14:107349540, 15:102531392\
            , 16:90354753, 17:81195210, 18:78077248, 19:59128983, 20:63025520, 21:48129895, 22:51304566}

def file_check(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {0} does not exist!".format(arg))
    else:
        return str(arg)

def makeRS_Autosomal(filename):

    if os.path.exists(filename+"_RS"):
        logger.info("{0}_RS is already exist!".format(filename))
        target_files_RS.append(filename+ "_RS")
        return

    new_file = open(filename+ "_RS", "w")
    target_files_RS.append(filename+ "_RS")

    with open(filename, "r") as f1:
        for line in f1:
            if line.startswith("rs"):

                a = [x.strip() for x in line.split('\t')]

                if re.search(r'[a-zA-Z]', a[1]) or a[3] == "--" or a[3] == "II":
                        continue
                new_file.write(line)
    new_file.close()
    logger.info("Make rsXXXX and Autosomal only ---------- [OK]")

#-------------------------------------------------------------------
def makeEqualize():
    filelist = []
    genlist = []
    rs_list = []

    # Open all RS files & Make generator list
    for i in range(len(target_files_RS)):
        with open(target_files_RS[i], "r") as f1:

            # Make generator for each file handler
            it1 = (line.split('\t') for line in f1)
            it2 = (word[0] for word in it1)
            genlist.append(it2)

            # Make list which has only rsXXXXX to make common set
            temp_list= []
            while True:
                try :
                    temp_list.append(next(genlist[i]))
                except StopIteration:
                    break
            rs_list.append(temp_list)

    # Make common RS set
    commonSet = set(rs_list[0])
    for s in rs_list[1:]:
        commonSet.intersection_update(s)

    logger.info("Make Common Set between files---------- [OK]")

    # Make final normalized data with block
    block_index = 0

    for RSfile in target_files_RS:
        line_cnt = 0

        with open(RSfile, "r") as source:
            Final_NormalizedData = {}
            for line in source:
                a = [x.strip() for x in line.split('\t')]
                if a[0] in commonSet:
                    # all uppercase and make it sorted to compare
                    genomeData = ''.join(sorted(a[3].upper()))
                    line_cnt += 1

                    block_index = int(a[2]) // 10000000

                    try:
                        Final_NormalizedData[int(a[1])][block_index] += [genomeData]

                    except IndexError: #new block_index

                        block_gap = block_index - len(Final_NormalizedData[int(a[1])])
                        if block_gap >= 1: # block_index skiped
                            for i in range(block_gap+1):
                                Final_NormalizedData[int(a[1])] += [[]]
                        else:
                            Final_NormalizedData[int(a[1])] += [[]]

                        Final_NormalizedData[int(a[1])][block_index] = [genomeData]

                    except KeyError: #new data

                        if block_index >= 1:
                            Final_NormalizedData[int(a[1])] = [[]]
                            for i in range(block_index):
                                Final_NormalizedData[int(a[1])] += [[]]
                        else:
                            Final_NormalizedData[int(a[1])] = [[]]

                        Final_NormalizedData[int(a[1])][block_index] = [genomeData]

            logger.info("Make Final Data structure (%s)  -- [OK]" %(RSfile))
            FinalDataSet.append(Final_NormalizedData)

    return "OK"

def findCopy(izipped_bindata):

    result = Counter()
    for pair in izipped_bindata:

        setData = set(pair[0]) & set(pair[1])
        if pair[0] == pair[1]:
            copy = "both"
        elif len(setData) == 1:
            copy = "one"
        else:
            copy = "neither"

        result[copy] += 1
    return result["both"], result["one"], result["neither"]

def calculateMatchingPercentage(outFile):
    """
    Calcuate the matching percentage for each binblocks

    Args:
        threshold : the value for making a decision 
    """
    source1_dict = FinalDataSet[0]
    source2_dict = FinalDataSet[1]

    result_dict = {}
    simple_result = {}
    scoring_dict = {}
    decision_list = []
    result_decision = defaultdict(list)

    for key, value in source1_dict.iteritems():
        outFile.write("\nchromosome #"+str(key)+"\n")
        outFile.write("block    both    one    neither    decision\n")

        chro_sum = 0
        one_sum = 0
        neither_sum = 0
        for binblock in range(len(value)):
            bindata = izip(value[binblock], source2_dict[key][binblock])
            (bothcopy, onecopy, neithercopy) = findCopy(bindata)

            result_dict[binblock+1] = [bothcopy, onecopy, neithercopy]
            #print result_dict[binblock+1]
            #simple_result[binblock+1] = [onecopy/float(bothcopy+onecopy+neithercopy), neithercopy/float(bothcopy+onecopy+neithercopy)]

        dict_sum = sum([sum(x) for x in list(result_dict.values())])
        centromere_criteria = int((dict_sum / float(len(result_dict))) * 0.5)

        scoring_list = []

        for binblock, result_data in result_dict.iteritems():

            if sum(result_data) < centromere_criteria:
                scoring_list.append(-1)
                decision = "ignore"
                outFile.write("%d\t%d\t%d\t%d\t%s\n" %(binblock, result_data[2], result_data[1], result_data[0], decision))
                continue
            chro_sum += result_data[2] + result_data[1] + result_data[0]
            one_sum += result_data[1]
            neither_sum += result_data[0]
            if result_data[2] != 0 and (result_data[2] / float(sum(result_data)) * 100) > args.neither_decision_ratio:
                scoring_list.append(0)
                decision = "neither"
            else:
                if result_data[1] == 0:
                    scoring_list.append(2)
                    decision = "both"
                elif result_data[1] > 0 and (result_data[1] / float(sum(result_data)) * 100) < args.one_drop_ratio:
                    scoring_list.append(2)
                    decision = "both"
                else:
                    scoring_list.append(1)
                    decision = "one"

            outFile.write("%d\t%d\t%d\t%d\t%s\n" %(binblock, result_data[2], result_data[1], result_data[0], decision))

        one_result_dict = {}
        #logger.debug(one_sum/float(chro_sum))
        scoring_dict[key] = [scoring_list.count(2), scoring_list.count(1), scoring_list.count(0), one_sum/float(chro_sum), neither_sum/float(chro_sum)]
        #logger.debug(scoring_dict)
    #logger.debug(scoring_dict)

    Final_List = []
    for key, value in scoring_dict.iteritems():
        if value[0] > sum(value) * (args.both_ack_ratio/float(100)):
            Final_List.append(2)
        elif value[1] > sum(value) * (args.one_ack_ratio/float(100)):
            Final_List.append(1)
        else:
            Final_List.append(0)

    OneRatio = sum(value[3] for key, value in scoring_dict.iteritems())/float(len(Final_List))*100
    NeitherRatio = sum(value[4] for key, value in scoring_dict.iteritems())/float(len(Final_List))*100

    #print Final_List
    print "\n[Final Decision]"
    print "Both    %.2f" %(Final_List.count(2)/float(len(Final_List))*100.0)
    print "One    %.2f" %(Final_List.count(1)/float(len(Final_List))*100.0)
    print "Neither    %.2f" %(Final_List.count(0)/float(len(Final_List))*100.0)
    print "\n[Matching Percentage]"
    print "Matching(One)    %.2f" %(OneRatio)
    print "Matching(Neither)    %.2f" %(NeitherRatio)

    outFile.write("\n[Final Decision]\n")
    outFile.write("Both    %.2f\n" %(Final_List.count(2)/float(len(Final_List))*100.0))
    outFile.write("One    %.2f\n" %(Final_List.count(1)/float(len(Final_List))*100.0))
    outFile.write("Neither    %.2f\n" %(Final_List.count(0)/float(len(Final_List))*100.0))
    outFile.write("\n[Matching Percentage]")
    outFile.write("Matching(One)    %.2f" %(OneRatio/float(len(Final_List))))
    outFile.write("Matching(Neither)    %.2f\n" %(NeitherRatio/float(len(Final_List))))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.format(author=__author__, email=__email__), formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

    pgroup = parser.add_argument_group("Input")
    pgroup.add_argument('f1', metavar='in-file1', type=lambda x:file_check(parser, x), help='23andMe format')
    pgroup.add_argument('f2', metavar='in-file2', type=lambda x:file_check(parser, x), help='23andMe format')
    pgroup.add_argument('output', metavar='output-file', help='detailed')

    ogroup = parser.add_argument_group("Options for Family")
    ogroup.add_argument('-nd', dest='neither_decision_ratio', default=1.0, type=float)
    ogroup.add_argument('-od', dest='one_drop_ratio', default=1.0, type=float)
    ogroup.add_argument('-bar', dest='both_ack_ratio', default=10.0, type=float)
    ogroup.add_argument('-oar', dest='one_ack_ratio', default=20.0, type=float)
    ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
    ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')
    ogroup.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    makeRS_Autosomal(args.f1)
    makeRS_Autosomal(args.f2)

    result = makeEqualize()
    if result != "OK":
        logger.ERROR("Error Occurred!")
        sys.exit(0)

    with open(args.output, "w") as output:
        calculateMatchingPercentage(output)
