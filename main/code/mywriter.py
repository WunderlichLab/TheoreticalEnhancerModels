# This script writes CRN files to be run
# for the CERENA Matlab package.

# Important notes:
# - keep mRNA 'R' to be the last species 
# - don't put any species between the TFs and R

import string
import math
import itertools
import sys 
import time
from collections import defaultdict 
import os 
import time


# Inputs:
# dupBool: the bsT1 and bsT2 arguments become bsT1 and bsT2
# per enhancer.
# binding sites 
# bsT1: total number of binding sites for T1
# bsT2: total number of binding sites for T2
# numEnhancers: how many enhancers
# fileTag: model tag number
# adairBool: to generate an Adair version of the model
# multiBinding: makes TFs share all their binding sites

def FileGenerator(dupBool,bsT1,bsT2,numEnhancers,fileTag,adairBool, \
                    multiBinding, extraModel):
    print('')
    print('')
    print('')

    print('fileTag:')
    print(fileTag)
    print('binding sites T1:')
    print(bsT1)
    print('binding sites T2:')
    print(bsT2)
    print('enhancers:')
    print(numEnhancers)

    rateNames = [] 


    if bsT1 == 0:
        rateNames = rateNames + ['kon' + str(2)]
        rateNames = rateNames + ['koff' + str(2)]
        rateNames = rateNames + ['r' + str(2)]
        numTFs = 1

    elif bsT2 == 0:
        rateNames = rateNames + ['kon' + str(1)]
        rateNames = rateNames + ['koff' + str(1)]
        rateNames = rateNames + ['r' + str(1)]
        numTFs = 1
    else:
        for num in range(1,3):
            rateNames = rateNames + ['kon' + str(num)]
            rateNames = rateNames + ['koff' + str(num)]
            rateNames = rateNames + ['r' + str(num)]
        numTFs = 2


    # adding mRNA degradation
    rateNames = rateNames + ['myalpha']     


    if bsT1 == 0:
        rateNames = rateNames + ['beta' + str(2)]
        rateNames = rateNames + ['betam' + str(2)]
    elif bsT2 == 0:
        rateNames = rateNames + ['beta' + str(1)]
        rateNames = rateNames + ['betam' + str(1)]
    else:
        for num in range(1,numTFs + 1):
            rateNames = rateNames + ['beta' + str(num)]
            rateNames = rateNames + ['betam' + str(num)]

    otherNames = ['Omega', 'time']

    # string for file name about fluc TF

    if multiBinding:
        multiStr = 'Multi'
    else:
        multiStr = ''

    if adairBool: 
        adairStr = 'Adair'
    else:
        adairStr = ''

    if dupBool: 
        dupStr = 'Dup'
    else:
        dupStr = ''

    if not adairBool and not dupBool:
        dupStr = 'Normal'


    bsT1Str = str(bsT1)
    bsT2Str = str(bsT2)

    # enhancer names are just drawn from a sorted alphabet
    alphabet_string = string.ascii_uppercase

    enhancerLetters = list(alphabet_string) 

    # using number of enhancers (should be less than 26)
    enhancerLettersUsed = enhancerLetters[:numEnhancers]

    totalBs = bsT1 + bsT2
    numRuns = 1
    sortPepper = False
    if numEnhancers > 1 and numEnhancers < totalBs and \
        bsT1 > 0 and bsT2 > 0 and not dupBool:
        sortPepper = True
        # we generate the sorted and peppered files
        numRuns = 2

    for runIndex in range(numRuns):
        if runIndex == 0:
            fileTagStr = str(fileTag)
        else:
            fileTagStr = str(fileTag + 1)

        # defining string of file name 
        fileName = 'modelDef_' + fileTagStr + multiStr + dupStr + \
                    adairStr + '_' + bsT1Str + bsT2Str + str(numEnhancers)

        fileName += '.m'

        # creating all possible permutations to add to the enhancer letters later
        permsTF = []

        # dup case. bsTi is binding sites of Ti per enhancer
        if dupBool:
            T1Array = []
            T2Array = []
            for i in range(numEnhancers):
                T1Array.append(bsT1)
                T2Array.append(bsT2)

        # non-dup case. bsTi is total binding sites of Ti 
        else:
            # building T1 and T2 arrays
            if sortPepper:
                T1Array = []
                T2Array = []
                # sort
                if runIndex == 0:
                    # figuring out if there are enough
                    # to just assign all T1s or T2s to
                    # the first enhancer.
                    T1Array = [0] * numEnhancers
                    T2Array = [0] * numEnhancers

                    # assign all bsT2 to the first enhancer
                    # if there are enough bsT1 for 
                    # the remaining enhancers
                    if bsT1 >= bsT2:
                        if bsT1 >= numEnhancers - 1:
                            T2Array[0] = bsT2
                            # starting from the second enhancer
                            tracker = 1
                            for i in range(bsT1):
                                T1Array[tracker] += 1 
                                tracker += 1
                                if tracker == numEnhancers: 
                                    tracker = 1
                        else:
                            # how many enhancers need
                            # single T2 binding sites?
                            helpEnhancers = numEnhancers - bsT1
                            T2Array[0] = bsT2 - helpEnhancers
                            tracker = 1
                            # this already makes the assumption
                            # that there cannot be empty enhancers.
                            # i.e. that totalBs >= numEnhancers
                            for i in range(helpEnhancers):
                                T2Array[tracker] += 1  
                                tracker += 1

                            # starting from the enhancer that
                            # doesn't need sites. The first enhancer that got
                            # the majority of T2 sites is index 0
                            tracker = 1 + helpEnhancers
                            for i in range(bsT1):
                                T1Array[tracker] += 1 
                                tracker += 1
                                if tracker == numEnhancers: 
                                    tracker = 1 + helpEnhancers     


                    # assign all bsT1 to the first enhancer
                    # if there are enough bsT2 for 
                    # the remaining enhancers
                    # bsT1 < bsT2
                    else:
                        if bsT2 > numEnhancers - 1:
                            T1Array[0] = bsT1
                            # starting from the second enhancer
                            tracker = 1
                            for i in range(bsT2):
                                T2Array[tracker] += 1 
                                tracker += 1
                                if tracker == numEnhancers: 
                                    tracker = 1
                        else:
                            # how many enhancers need
                            # single T1 binding sites?
                            helpEnhancers = numEnhancers - bsT2
                            T1Array[0] = bsT1 - helpEnhancers
                            tracker = 1
                            # this already makes the assumption
                            # that there cannot be empty enhancers.
                            # i.e. that totalBs >= numEnhancers
                            for i in range(helpEnhancers):
                                T1Array[tracker] += 1  
                                tracker += 1

                            # The first enhancer that got
                            # the majority of T2 sites is index 0
                            tracker = 1 + helpEnhancers
                            for i in range(bsT2):
                                T2Array[tracker] += 1 
                                tracker += 1
                                if tracker == numEnhancers: 
                                    tracker = 1 + helpEnhancers 


                elif runIndex == 1:
                    T1Array = [0] * numEnhancers
                    T2Array = [0] * numEnhancers
                    # in both cases if there are more
                    # binding sites for a Ti than enhancers
                    # we just restart the tracker.

                    # giving all T1s to the first
                    # bsT1 enhancers starting from
                    # the first enhancer
                    tracker = 0
                    for i in range(bsT1):
                        T1Array[tracker] += 1 
                        tracker += 1
                        if tracker == numEnhancers: 
                            tracker = 0 

                    # giving all T2s to the first
                    # bsT2 enhancers starting from the 
                    # last enhancer
                    tracker = numEnhancers - 1
                    for i in range(bsT2):
                        T2Array[tracker] += 1 
                        tracker -= 1
                        if tracker == -1: 
                            tracker = numEnhancers - 1  


                else:
                    print('Something is wrong with numRuns')
                    return                      

            # normal case 
            else:
                if numEnhancers == 1:
                    T1Array = [bsT1]
                    T2Array = [bsT2]
                elif numEnhancers == totalBs:
                    T1Array = []
                    T2Array = []
                    for i in range(bsT1):
                        T1Array.append(1)
                        T2Array.append(0)

                    for i in range(bsT2):
                        T1Array.append(0)
                        T2Array.append(1)
                # distribute as evenly as possible since this case
                # implies that they are all identical
                elif numEnhancers > 1 and numEnhancers < totalBs:
                    # initializing arrays 
                    T1Array = [0] * numEnhancers
                    T2Array = [0] * numEnhancers
                    if bsT1 == 0:
                        tracker = 0
                        for i in range(bsT2):
                            T2Array[tracker] += 1 
                            tracker += 1
                            if tracker == numEnhancers: 
                                tracker = 0
                    elif bsT2 == 0:
                        tracker = 0
                        for i in range(bsT1):
                            T1Array[tracker] += 1 
                            tracker += 1
                            if tracker == numEnhancers: 
                                tracker = 0
                    else: 
                        print('Something is wrong with bsT1 and bsT2 \
                             in the not sAndP case')
                else:
                    print('Something is wrong with T1/T2Array')
                    return


        # We assign TFs for each enhancer.
        # T1Array specifies how many T1 bind to enhancer 1, then
        # enhancer 2, .... e.g. [T1inEnhancer1, T1inEnhancer2,...]
        # T2Array does the same for T2.
        if adairBool:
            permsTFtemp = []
            enhancerTracker = 0
            for i in range(numEnhancers):
                permsTFtemp = []
                # in Adair, the number of binding sites 
                # determines the upper bound of the subscript
                holder = []
                if T1Array[enhancerTracker] != 0:
                    for k in range(T1Array[enhancerTracker] + 1):
                        holder.append(k)
                    permsTFtemp.append(holder)

                holder = []
                if T2Array[enhancerTracker] != 0:
                    for k in range(T2Array[enhancerTracker] + 1):
                        holder.append(k)
                    permsTFtemp.append(holder)

                enhancerTracker += 1
                # print(permsTFtemp)
                # time.sleep(1)

                # creating the combinations of all the list of lists where
                # each element corresponds to the permutation of a TF 
                permsTFtemp = list(itertools.product(*permsTFtemp))
                permsTF.append(permsTFtemp)     

                # print(permsTF)    
                # time.sleep(1)

        # general model with all sites explicit
        else: 

            # holds the TF permutations for each enhancer
            # before being added to the main list which is 
            # permsTF
            permsTFtemp = []
            trackerTF   = 1
            # these are the TF numbers for each of the remaining 
            # enhancers 
            # iterating over enhancer
            # note that for individualBool we already know how many
            # enhancers are because TArrays specify num binding sites 
            # for T1 and T2 in each enhancer
            enhancerTracker = 0
            for i in range(numEnhancers):
                permsTFtemp = []

                holder = [0, 1]

                if T1Array[enhancerTracker] != 0:
                    permsTFtemp.append([p for p in \
                        itertools.product(holder, \
                            repeat = T1Array[enhancerTracker])])

                holder = [0, 2]

                if T2Array[enhancerTracker] != 0:
                    permsTFtemp.append([p for p in \
                        itertools.product(holder, \
                            repeat = T2Array[enhancerTracker])])

                enhancerTracker += 1 

                # creating the combinations of all the list of lists where
                # each element corresponds to the permutation of a TF 
                permsTFtemp = list(itertools.product(*permsTFtemp))
                permsTF.append(permsTFtemp)

            # print(permsTF)
            # time.sleep(10)

        # writing enhancer names with their possible binding 
        # configurations (numBindingSites is the number of binding
        # sites PER TF.
        enhancerNames = []
        # defining the name for the TFs
        tfNames = []
        # creating a dictionary for each enhancer where each TF configuration
        # of the enhancer in a string form such as 'X_{x1,x2,x3...}' is
        # a key for the rank of such enhancer (defined as the number of nonzero
        # elements in the enhancer subscript)
        enhancerDicts = []
        # creating one big dictionary with all enhancers names
        # (useful for later)
        enhancerBigDict = defaultdict(list)
        # creating a dictionary that just lists the enhancers alphabetical 
        # order in numbers so that we can associate rates in case 
        # diffentRates is given as an option 
        miniEnhancerDict = defaultdict(list)
        # creating an enhancer dictionary that stores 
        # a list of the indices of nonzero values in
        # the enhancer
        locIndexDicts = [] 
        # stores which Tfs belong to which enhancer. If splitTF = False,
        # then every enhancer has every TF. As such, length is numEnhancers
        enhancerTFAssocs = [] 
        # stores which TFs bind to which sites. As such, it is a list of lists of length numTFs if
        # splitTF is false OR if splitTF is true then it is a list (length enhancerNums) of lists
        # (length TFnums/numEnhancers (where first enhancer gets the surplus if odd)) 
        # of lists (length binding sites for this TF)
        # if splitTF = False then [[1,2,3],[4,5,6],....]
        # if splitTF = True then [[[1,2,3],[4,5,6]],[[1,2,3]]....]
        TFLocsAssocs = []

        # making a list with all TF names 
        if bsT1 == 0:
            tfNames.append('T2')
        elif bsT2 == 0:
            tfNames.append('T1')
        else:
            for i in range(1,numTFs + 1):
                tfNames.append('T' + str(i))

        # the ith dictionary corresponds to the ith enhancer
        for i in range(numEnhancers):
            enhancerDicts.append(defaultdict(list))
            locIndexDicts.append(defaultdict(list))

        # creating enhancerTFAssocs which is a list of lists describing which enhancer owns
        # which TFs. When splitTF = false then all enhancers own all TFs and therefore we
        # just need a list of all tfs for each enhancer.
        # tracks the TFs in tfNames
        for enhancer in range(numEnhancers):
            tempTFList = []

            if T1Array[enhancer] != 0:
                tempTFList.append('T1')             
    
            if T2Array[enhancer] != 0:
                tempTFList.append('T2')             

            enhancerTFAssocs.append(tempTFList)

        # print(enhancerTFAssocs)
        # time.sleep(5)

        if not adairBool:
            # creating TFLocsAssocs which is a list of lists of lists that describes 
            # the locations/indices of the binding sites
            # that belong to each TF for each enhancer, the first element is for E_1 and within are
            # the elements for T_1, then T_2 and so on until we exhaust the TFs for E_1 and the next
            # element would be the TFs for E_2 and the binding sites that they own but we need to reset
            # the tracker to count from 0 since these TFs are exclusive to the enhancer in splitTF = True
            for enhancer in range(numEnhancers):
                tracker = 0
                # holds the list of lists of TFs and their sites 
                # for this enhancer
                tempTFs = []

                if T1Array[enhancer] != 0:
                    tempSites = []
                    for site in range(T1Array[enhancer]):
                        tempSites.append(tracker)
                        tracker += 1
                    tempTFs.append(tempSites)               
        
                if T2Array[enhancer] != 0:
                    tempSites = []
                    for site in range(T2Array[enhancer]):
                        tempSites.append(tracker)
                        tracker += 1
                    tempTFs.append(tempSites)           

                TFLocsAssocs.append(tempTFs)

        # making the enhancer letter names.
        # going over permutations for these enhancer. permsTF in this case
        # should be a meta-perms TF in the case of splitTF = True
        tracker = 0
        # going over enhancer names
        for letter in enhancerLettersUsed:
            # stores the letter of the enhancer and returns a number in case
            # of different rates
            miniEnhancerDict[letter] = str(tracker) 
            permsTFtemp = permsTF[tracker] 
            # print(permsTFtemp)
            # time.sleep(1)
            if adairBool:
                for permType in permsTFtemp: 
                    # holds the str of the TF configuration
                    if type(permType) is tuple:
                        holder = permType
                    else:
                        # case where there is only 1 TF in the enhancer
                        holder = [permType]

                    numStr = ''

                    # list of nonzero subcript locations in enhancer config
                    # at first stores 1 or 0 if TF bound of not respectively
                    locTracker = []
                    # sums the nonzero elements in this TF configuration
                    # to form the rank
                    numTracker = 0
                    # in the Adair case the subscript gives the number
                    # of TFs of a given kind currently bound
                    for num in holder:
                        numStr = numStr + str(num)
                        if num != 0:
                            numTracker += num
                            locTracker.append(num)
                        else:
                            locTracker.append(num)

                    # getting the name of this configuration
                    enhancerNames.append(letter + numStr)
                    enhancerDicts[tracker][letter + numStr] = numTracker
                    enhancerBigDict[letter + numStr] = numTracker 
                    # list of all indices in this enhancer
                    locIndexDicts[tracker][letter + numStr] = locTracker

            else:
                for permType in permsTFtemp: 
                    holder = [item for sublist in permType for item in sublist] 
                    # holds the str of the TF configuration
                    numStr = ''
                    # sums the nonzero elements in this TF configuration
                    # to form the rank
                    numTracker = 0
                    # list of nonzero subcript locations in enhancer config
                    # at first stores 1 or 0 if TF bound of not respectively
                    locTracker = []
                    for num in holder:
                        numStr = numStr + str(num)
                        if num != 0:
                            numTracker += 1
                            locTracker.append(1)
                        else:
                            locTracker.append(0)

                    # getting the name of this configuration
                    enhancerNames.append(letter + numStr)
                    #  adding the enhancer's rank to the dictionaries
                    enhancerDicts[tracker][letter + numStr] = numTracker
                    enhancerBigDict[letter + numStr] = numTracker 
                    # getting list of indices that are nonzero (= 1) 
                    # for this enhancer config and storing them in a list
                    nonZeroIndeces = [k for k, e in enumerate(locTracker) if e != 0]
                    # adding the list of nonzero indices that 
                    locIndexDicts[tracker][letter + numStr] = nonZeroIndeces
            tracker += 1

        enhancerNames = list(dict.fromkeys(enhancerNames))
        # print(enhancerNames) 
        # time.sleep(10)

        # opening the file to be created
        file = open("models" + extraModel + "/" + fileName,'w+') 

        # start writing up file
        # writing the symbolic variables
        # holds all symbolic variables 
        symsStr = 'syms '
        # adding all enhancer names
        for name in enhancerNames:
            symsStr += name + ' '       
        # adding all tf names
        for name in tfNames:
            symsStr += name + ' '
        # adding reaction rate names
        for name in rateNames:
            symsStr += name + ' ' 
        # adding all remaining parameters
        for name in otherNames:
            symsStr += name + ' ' 

        file.write(symsStr + 'R' + '\n')
        file.write('\n')
        file.write('System.time = time;\n')
        file.write('System.compartments = {\'cell\'};\n')
        file.write('System.volumes = [Omega];\n')

        # holds all the species with semicolons in between
        speciesStr = ''
        # holds compartments
        compartmentsStr = ''
        # holds the types (stochastic or moment)
        typeStr = ''
        # holds minimum/maximum values for each species
        minStr = ''
        maxStr = ''
        # holds mean string
        mu0Str = ''
        # temp tracker 
        trackerTF = True
        for name in enhancerNames:
            speciesStr += name + '; '
            compartmentsStr += '\'cell\'; '
            typeStr += '\'stochastic\';'
            minStr += '0; '  
            maxStr += '1; '
            # if the rank is 0 the enhancer has no TFs bound.
            # we assume that all enhancers start with no TFs bound
            if enhancerBigDict[name] == 0:
                mu0Str += '1; '
            else:
                mu0Str += '0; '

        # Always have TFs at the end followed by R
        # adding elements of T
        for name in tfNames:
            speciesStr += str(name) + '; '
            compartmentsStr += '\'cell\'; '
            typeStr += '\'stochastic\';'
            minStr += '0; '  
            # this doesn't matter for MM methods
            maxStr += '10; '
            # all tfs just start at 0
            mu0Str += '0; '

        # adding the elements for R
        speciesStr += 'R'
        compartmentsStr += '\'cell\''
        typeStr += '\'moment\''
        minStr += '0' 
        maxStr += '150'
        mu0Str += '0'

        file.write('System.state.variable = [' + speciesStr + '];\n')
        file.write('System.state.compartment = {' + compartmentsStr + '};\n')
        file.write('System.state.type = {' + typeStr + '};\n')
        file.write('System.state.xmin = [' + minStr + '];\n')
        file.write('System.state.xmax = [' + maxStr + '];\n')
        file.write('System.state.mu0 = [' + mu0Str + '];\n')
        file.write('System.state.C0 = zeros(length(System.state.variable)' + \
                '*(length(System.state.variable)+1)/2,1);\n')


        # finding constraints (enhnacers can only be in one state bound at a time)
        # holds list of constraints (all the indexes in the System.state variable that
        # should add up to 1) (should be n elements for n enhancers)
        listIndexes = []
        counter = 1
        indexer = 0
        for indexDict in enhancerDicts:
            listIndexes.append([])
            # in MATLAB indexes start at 1 
            for index in indexDict:
                listIndexes[indexer].append(counter) 
                counter += 1
            indexer += 1

        constraintStr = ''
        # creating the list of strings for each constraint
        for indexes in listIndexes:
            constraintStr += '('
            for index in range(len(indexes)):
                if index == len(indexes) - 1:
                    constraintStr += 'x(' + str(indexes[index]) + ')) == 1 && '
                else:
                    constraintStr += 'x(' + str(indexes[index]) + ') + '

        # removing the last && and the two whitespaces
        constraintStr = constraintStr[:-4]
        file.write('System.state.constraint = @(x) ('  + constraintStr + ');')
        file.write('\n')

        # writing system parameter rate variables
        rateStr = ''
        for name in rateNames:
            rateStr += name + '; ' 

        # removing the ; and whitespaces at the end
        rateStr = rateStr[:-2]
        file.write('System.parameter.variable = [' + rateStr + '];\n')
        file.write('System.kappa.variable = [Omega];\n')
        file.write('System.scaleIndicator = \'microscopic\';\n')

        file.write('\n')

        # writing the reactions 

        # defining all reactions based on a few principles:
        # 1. An enhancer A_{x1 x2 x3 ... xn} is assigned a rank
        # which is sum_{i = 1}^{n} (if x_i = 0, then 0, else x_i = 1)
        # 2. Except for when rank equals 0 or n (edge cases, need to be handled
        # differently but still simple implementation), then there exists 
        # a reaction between all species with rank X and those with rank X - 1
        # or X + 1  (in this case you might to include a '+ T' if fluc TF; also
        # we have to multiply by the number of available binding sites
        # at each step)
        # 3. The reaction for mRNA degradation (R -> 0) can be added manually
        # and production-degradation of TFs (T_i <> 0) can be added iteratively
        # 4. For now, k_on and k_off will be the same and a single value of 
        # mRNA production 'r'' will be defined. Then, the rank will be used to 
        # determine the value of mRNA production for each enhancer which
        # will be defined as r/i for enhancer with rank = i.

        # keeps track of reaction numbers (which have arbitrary
        # order but a unique numerical identifier is needed for each reaction)
        # (remember MATLAB indices start at 1).
        reactionID = 1
        #  keeps track of the list of lists of tfs for each enhancer 
        # where each top level element corresponds to an enhancer.
        tfAssocsTracker = 0
        # going over ranks and locations of nonzero indices in enhancer
        # subscript (each dicts is for a different enhancer)
        for mydict, myDictLoc in zip(enhancerDicts, locIndexDicts):
            # print(mydict)
            # print(myDictLoc)
            # time.sleep(1)
            # selecting which TFs are in control of this enhancer
            currentTFs = enhancerTFAssocs[tfAssocsTracker]
            # TFLocAssocs is a list of lists of lists (top level each enhancer)
            # if splitTF is true. But if splitTF is false, then all enhancers
            # own all TFs. In particular, TFLocAssocs describes the binding sites
            # controlled by each TF.
            if not adairBool:
                currentSitesOfTFs = TFLocsAssocs[tfAssocsTracker]

            # print("current sites of TFs:")
            # print(currentSitesOfTFs)
            # time.sleep(1)
            tfAssocsTracker += 1
            # going over each rank in the dictionary.
            # finding the largest rank possible; for now
            # it would be numTFs * numBindingSites but 
            # we will keep this as general as possible
            # since we want to account for spliTF case 
            # where the first enhancer might have more TFs
            # than others and thus a larger rank
            largestRank = max(mydict.values())
            for rank in range(largestRank + 1):
                # finding all species at the current rank
                currentRank = [k for k,v in mydict.items() if v == rank]
                # print("current rank species:")
                for currentSpecies in currentRank:
                    # we decompose the string of the current species to find
                    # the relevant TF from the index given differences sites
                    # below. This should also work for the non-mutiBinding 
                    # case but we keep it separate for now while testing.
                    # list on a string decomposes it into an array of 
                    # individual letters
                    decompCurrentSpecies = list(currentSpecies)
                    # getting rid of the letter since we will use this
                    # only to get the relevant TF
                    currentLetter = decompCurrentSpecies.pop(0) 
                    # for differentRates we will just assign rate k_{i-1} to 
                    # enhancer i. All enhancers already exist in their own linkage
                    # classes
                    # rateNumber = miniEnhancerDict[currentLetter]
                    # print("current species:")
                    # print(currentSpecies)
                    # time.sleep(1)
                    # check which location is different between two lists python
                    indicesCurrent = myDictLoc[currentSpecies]
                    # print("indices current:")
                    # print(indicesCurrent)
                    # time.sleep(5)
                    # collecting all species that 1 rank below (rank - 1)
                    # into a list
                    belowRank = [k for k,v in mydict.items() if v == rank - 1]
                    # creating all reactions that link a species with 
                    # k_on to all species in the current rank
                    if currentRank != 0:
                        for belowSpecies in belowRank:
                            print("below species:")
                            print(belowSpecies)
                         
                            # decomposing the current species below the current one.
                            decompBelowSpecies = list(belowSpecies)
                            decompBelowSpecies.pop(0)
                            # keeping only those in belowRank that are at most one TF different
                            # in their configuration
                            indicesBelow = myDictLoc[belowSpecies]
                            print("indices below:")
                            print(indicesBelow)
                          
                            if adairBool:
                                differencesSites = []
                                for index, (first, second) in enumerate(zip(indicesCurrent, indicesBelow)):
                                    if first != second:
                                        differencesSites.append(index)
                            else:                       
                                # finding TF locations/indeces that appear in either 
                                # the currentSpecies and the belowSpecies   
                                elementsOR = (set(indicesCurrent) | set(indicesBelow))
                                # print("elementsOR:")
                                # print(elementsOR)
                                # time.sleep(1)
                                elementsAND = (set(indicesCurrent) & set(indicesBelow))
                                # print("elementsAND:")
                                # print(elementsAND)
                                # time.sleep(1)
                                differencesSites = list(elementsOR - elementsAND)
                            # this is used for the kon rates which are 
                            # multiplied by the number of available sites  
                            if adairBool:
                                # finding the current TF
                                # in the Adair case differencesSites stores the index of 
                                # the subscript that we are updating
                                # T1Array specifies how many T1 bind to enhancer 1, then
                                # enhancer 2, .... e.g. [T1inEnhancer1, T1inEnhancer2,...]
                                # T2Array does the same for T2..
                                currentEnhancerNum = int(miniEnhancerDict[currentLetter])
                                currentBsT1 = T1Array[currentEnhancerNum]
                                currentBsT2 = T2Array[currentEnhancerNum]

                                if currentBsT1 == 0:
                                    relevantLargestRank = currentBsT2
                                elif currentBsT2 == 0:
                                    relevantLargestRank = currentBsT1
                                else:
                                    # case where we have both TF
                                    print('differences sites:')
                                    tempTf = differencesSites[0]
                                    print(differencesSites[0])
                                   
                                    if tempTf == 0:
                                        relevantLargestRank = currentBsT1
                                    else:
                                        relevantLargestRank = currentBsT2

                              
                                print('enhancer number:')
                                print(miniEnhancerDict[currentLetter])
                     

                                # finding the current enhancer through miniEnhancerDict
                                # that returns a number based on the enhancer latter
                                
                 
                                # the number of available sites to the below species 
                                # would be the total (given by the largest rank) minus
                                # how many are already taken (given by the index of the
                                # below species).
                                # turning this into an int to do the operation below
                                decompBelowSpeciesTemp = int(decompBelowSpecies[differencesSites[0]])
                                # the relevant largest rank is that which is for the current TF
                                adairNumber = relevantLargestRank - decompBelowSpeciesTemp
                                # turning back into string for printing it into the CRN file
                                adairNumber = str(adairNumber);
                            # if we are doing multiBinding then we have to make sure
                            # that every binding site is the same except for the one 
                            # we are going to update
                            multiBool = True
                            if not adairBool:
                                if len(differencesSites) == 1:
                                    # placeholder that changes if any site besides 
                                    # the one to be updated is different  
                                    multiBool = True
                                    # comparing each site except the one that will be 
                                    # updated.
                                    for siteIndex in range(0,len(decompCurrentSpecies)):
                                        if siteIndex != differencesSites[0]:
                                            if decompCurrentSpecies[siteIndex] != \
                                                 decompBelowSpecies[siteIndex]:
                                                 multiBool = False
                            # print("differences Sites:")
                            # print(differencesSites)
                            # time.sleep(1)
                            # this implies that all TFs are in the same location as the currentSpecies
                            # except for a single difference (a single displacement of the arrow in the 
                            # stoichiometric 2D plane)
                            if len(differencesSites) == 1 and multiBool:
                                if multiBinding:
                                    relevantTF = "T" + decompCurrentSpecies[differencesSites[0]]
                                elif adairBool:
                                    # in the Adair case differencesSites stores the index of 
                                    # the subscript that we are updating
                                    # print(differencesSites[0])
                                    # time.sleep(1)
                                    relevantTF = currentTFs[differencesSites[0]]

                                else:
                                    # getting the TF that corresponds to this binding site index
                                    # (note that differencesSites gives the index of the binding site
                                    # that is different between the currentSpecies and the belowSpecies).
                                    # To find which TF owns this site, we need to find the index of this
                                    # site in the top level list of currentTFs.
                                    # getting tuple where the first index is 
                                    # index in top level list and second index is index 
                                    # in the low level list
                                    # (note the distinction between currentSite and currentSites)
                                    reactionTFIndex = [(i, currentSite.index(differencesSites[0])) \
                                                        for i, currentSite in enumerate(currentSitesOfTFs)  \
                                                        if differencesSites[0] in currentSite]
                                    # print("reactionTFIndex")
                                    # print(reactionTFIndex)
                                    # time.sleep(1)
                                    # getting top level index so we know which TF owns this site
                                    reactionTFIndex = reactionTFIndex[0][0]
                                    # print("reactionTFIndex top level:")
                                    # print(reactionTFIndex)
                                    # time.sleep(1)                             
                                    relevantTF = currentTFs[reactionTFIndex]
                                # splitting the TF into two so that we can do the different rates       
                                rateNumber = list(relevantTF)
                                rateNumber.pop(0)
                                rateNumber = rateNumber[0]
                                # rates start at 0 while TFs start at 1
                                rateNumber = str(int(rateNumber));
                                # print("relevant TF:")
                                # print(relevantTF)
                                # time.sleep(1) 
                                file.write('System.reaction(' + str(reactionID) + ').educt = [' \
                                     + belowSpecies + ',' + relevantTF + '];\n')                    
                                file.write('System.reaction(' + str(reactionID) + ').product = ' \
                                    + currentSpecies + ';\n') 
                                # in adair models the rate constants are multiplied by the number of 
                                # available sites (kon) or occupied sites (koff)
                                if adairBool:
                                    file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                                        + adairNumber + ' * ' + 'kon' + rateNumber + ' * ' + belowSpecies + \
                                        ' * ' + relevantTF +';\n')                                
                                else:                  
                                    file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                                        + 'kon' + rateNumber + ' * ' + belowSpecies + ' * ' + relevantTF +';\n')                                

                                file.write('\n')
                                reactionID += 1

                    # collecting all species that are 1 rank above (rank + 1)
                    # into a list
                    aboveRank = [k for k,v in mydict.items() if v == rank + 1]
                    if currentRank != largestRank:
                        for aboveSpecies in aboveRank:
                            decompAboveSpecies = list(aboveSpecies)
                            decompAboveSpecies.pop(0)
                            # this is used for the koff rates which are 
                            # multiplied by the number of occupied sites  
                            indicesAbove = myDictLoc[aboveSpecies]
                            # print("indices current:")
                            # print(indicesCurrent)
                            # time.sleep(2)                 
                            # print("indices above:")
                            # print(indicesAbove)
                            # time.sleep(1)
                            if adairBool:
                                differencesSites = []
                                for index, (first, second) in enumerate(zip(indicesCurrent, indicesAbove)):
                                    if first != second:
                                        differencesSites.append(index)
                            else:
                                elementsOR = (set(indicesCurrent) | set(indicesAbove))
                                elementsAND = (set(indicesCurrent) & set(indicesAbove))
                                differencesSites = list(elementsOR - elementsAND)   

                            if adairBool:
                                # the number of occupied sites in the above species 
                                # would be how many are already taken 
                                # (given by the index of the below species)
                                adairNumber = decompAboveSpecies
                                # in the Adair case differencesSites stores the index of 
                                # the subscript that we are updating
                                adairNumber = adairNumber[differencesSites[0]]
                                adairNumber = str(int(adairNumber));                                                    
                            # print("elementsOR:")
                            # print(elementsOR)
                            # time.sleep(1)
                            # print("elementsAND:")
                            # print(elementsAND)
                            # time.sleep(1)
                            multiBool = True
                            if not adairBool:
                                if len(differencesSites) == 1:
                                    # placeholder that changes if any site besides 
                                    # the one to be updated is different  
                                    multiBool = True
                                    # comparing leach site except the one that will be 
                                    # updated.
                                    for siteIndex in range(0,len(decompCurrentSpecies)):
                                        if siteIndex != differencesSites[0]:
                                            if decompCurrentSpecies[siteIndex] != \
                                                 decompAboveSpecies[siteIndex]:
                                                 multiBool = False
                                        
                            # this implies that all TFs are in the same location as the currentSpecies
                            # except for a single difference (a single displacement of the arrow in the 
                            # stoichiometric 2D plane)
                            if len(differencesSites) == 1 and (multiBool or adairBool):
                                if multiBinding:
                                    # this also works in the general model. Remember that here species just 
                                    # look like A_112233...
                                    relevantTF = "T" + decompAboveSpecies[differencesSites[0]]
                                elif adairBool:
                                    # in the Adair case differencesSites stores the index of 
                                    # the subscript that we are updating
                                    relevantTF = currentTFs[differencesSites[0]]
                                else:
                                    reactionTFIndex = [(i, currentSite.index(differencesSites[0])) \
                                                        for i, currentSite in enumerate(currentSitesOfTFs)  \
                                                        if differencesSites[0] in currentSite]
                                    # getting top level index so we know which TF owns this site
                                    reactionTFIndex = reactionTFIndex[0][0]
                                    relevantTF = currentTFs[reactionTFIndex]

                                rateNumber = list(relevantTF)
                                rateNumber.pop(0)
                                rateNumber = rateNumber[0]
                                rateNumber = str(int(rateNumber));

                                file.write('System.reaction(' + str(reactionID) + ').educt = ' \
                                    + aboveSpecies + ';\n') 
                                # the relevant TF is released this time in the product reaction 
                                file.write('System.reaction(' + str(reactionID) + ').product = [' \
                                    + currentSpecies + ',' + relevantTF + '];\n')  
                                if adairBool:  
                                    file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                                    + adairNumber + ' * ' + 'koff' + rateNumber + ' * ' + \
                                             aboveSpecies + ';\n') 
                                else:
                                    file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                                    + 'koff' + rateNumber + ' * ' + aboveSpecies + ';\n')                                     
                                file.write('\n')
                                reactionID += 1 

                    # adding reactions for mRNA production
                    if rank != 0:
                        #  adding regular additive reactions
                        file.write('System.reaction(' + str(reactionID) + ').educt = ' \
                            + currentSpecies + ';\n')                   
                        file.write('System.reaction(' + str(reactionID) + ').product = [' \
                            + currentSpecies + ',R];\n')                    
                        # we assume that transcription rates are proportional
                        # to rank normalized
                        # file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                        #   + str(rank) + '/' + str(largestRank) + ' * r * ' + currentSpecies + ';\n')
                        # new edit, just any bound enhancer state makes mRNA at rate r
                        # this is the case where the enhancer makes mRNA at a rate 
                        # that is additive to all the TFs that are bound e.g. given
                        # A112203, the rate would be r1 + r1 + r2 + r2 + r3
                        bigR = '('
                        if adairBool:
                            # the rate at which this species will make mRNA
                            # is the number of TF times its corresponding rate.
                            # e.g. T1num*r0 + T2num*r1 + ....
                            # note that if you are splitting then it does not have
                            # to start at T1. e.g. T2num * r1 + T3num*r2 ....
                            # going over each TF
                            subscriptTracker = 0
                            for tf in currentTFs:
                                decompCurrentTF = list(tf) 
                                decompCurrentTF.pop(0)
                                rateNum = str(int(decompCurrentTF[0]))

                                tfNum = decompCurrentSpecies[subscriptTracker]
                                
                                subscriptTracker += 1
                                
                                bigR = bigR + tfNum + ' * ' + 'r' + rateNum + ' + '  
                                
                        else:
                            for tfBound in decompCurrentSpecies:
                                if tfBound != '0':
                                    tfBoundHolder = str(int(tfBound));
                                    bigR = bigR + 'r' + tfBoundHolder + ' + ' 
                                    # getting rid of the last + and space
                        
                        bigR = bigR[:-3]
                        bigR = bigR + ')'

                        file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                        + bigR + ' * ' + currentSpecies + ';\n')                                               

                        file.write('\n')
                        reactionID += 1
        #  adding TF creation and degration 
        # going over all TFs T_i
        # keeps track of the tf loop in case we are using different rates
        trackerTfRate = 1
        for tf in tfNames:
            # degradation
            file.write('System.reaction(' + str(reactionID) + ').educt = ' \
            + tf + ';\n')                   
            file.write('System.reaction(' + str(reactionID) + ').product = [];\n')      

            if bsT1 == 0:
                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                        + 'betam' + str(2) + ' * ' + tf + ';\n')    
            elif bsT2 == 0: 
                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                        + 'betam' + str(1) + ' * ' + tf + ';\n')   
            else:
                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                    + 'betam' + str(trackerTfRate) + ' * ' + tf + ';\n')    
            file.write('\n')
            reactionID += 1

            # creation
            file.write('System.reaction(' + str(reactionID) + ').educt = ' \
            + '[]' + ';\n')                 


            numTFsBurstArray = [4,12]
            # making the TF string [T,T,....,T] so that TFs come in a burst if required
            tfBurst = '[' + tf;
            if bsT1 == 0:
                for burstIndex in range(1,12):
                    # tfs are added to tfNames in numerical order so [T1, T2, T3,....,Tn]
                    tfBurst = tfBurst + ', ' + tf 
                tfBurst = tfBurst + ']'

                file.write('System.reaction(' + str(reactionID) + ').product = ' + tfBurst +  ';\n')    

                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                    + 'beta' + str(2) + ';\n')  

                file.write('\n')
            
            elif bsT2 == 0:
                for burstIndex in range(1,4):
                    # tfs are added to tfNames in numerical order so [T1, T2, T3,....,Tn]
                    tfBurst = tfBurst + ', ' + tf 
                tfBurst = tfBurst + ']'

                file.write('System.reaction(' + str(reactionID) + ').product = ' + tfBurst +  ';\n')    

                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                    + 'beta' + str(1) + ';\n')  

                file.write('\n')

            else:
                for burstIndex in range(1,numTFsBurstArray[trackerTfRate - 1]):
                    # tfs are added to tfNames in numerical order so [T1, T2, T3,....,Tn]
                    tfBurst = tfBurst + ', ' + tf 
                tfBurst = tfBurst + ']'


                file.write('System.reaction(' + str(reactionID) + ').product = ' + tfBurst +  ';\n')    

                file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
                    + 'beta' + str(trackerTfRate) + ';\n')  
                file.write('\n')
                trackerTfRate = trackerTfRate + 1

            reactionID += 1

        # adding reaction for mRNA degradation
        file.write('System.reaction(' + str(reactionID) + ').educt = ' \
            + 'R' + ';\n')                  
        file.write('System.reaction(' + str(reactionID) + ').product = [];\n')                  
        file.write('System.reaction(' + str(reactionID) + ').propensity = ' \
            + 'myalpha * R;\n')
        file.write('\n')

        # printing output variable (not necessary for now so we just set it equal 
        # to mRNA but a lot of the methods require it
        file.write('System.output.variable = [R];\n')
        file.write('System.output.function = [R];\n')

        file.close()

# single file generation
# dupBool = False
# bsT1 = 1
# bsT2 = 1
# numEnhancers = 2
# fileTag = 1
# adairBool = False
# multiBinding = False
# FileGenerator(dupBool,bsT1,bsT2,numEnhancers,fileTag, \
#             adairBool,multiBinding, extraModel)

ready = True
if ready:
    # file generation 
    myDup = True;
    # not implemented yet
    multiBinding = False

    if not myDup:
        # DupBool = False section: bsTi is
        # total binding sites for Ti
        dupBool = False
        # iterating over 0 <= bsT1,bsT2 <= 4 (1 <= bsTotal <= 4), 1 <= numEnhancers <= 4
        # and also generating an Adair version of every model starting at +200
        for bsT1 in range(5):
            for bsT2 in range(5):
                totalBs = bsT1 + bsT2
                if totalBs < 5 and totalBs > 0:
                    for numEnhancers in range(1,5):
                        # cannot have empty enhancers
                        if totalBs >= numEnhancers:
                            # doing 1 Adair per model
                            for multiAdairInd in range(2):
                            # for generating only Adair
                            # for multiAdairInd in [1]:
                                # reading the current file tag
                                tagFile = open('tagTracker.txt','r') 
                                fileTag = int(tagFile.read())
                                tagFile.close()

                                if multiAdairInd == 1:
                                    adairBool = True
                                    tempFileTag = fileTag + 200 
                                else:
                                    adairBool = False
                                    tempFileTag = fileTag

                                if adairBool:
                                    extraModel = "Adair" 
                                else:
                                    # directory path 
                                    extraModel = "" 

                                FileGenerator(dupBool,bsT1,bsT2,numEnhancers,tempFileTag, \
                                            adairBool,multiBinding, extraModel)

                            # fileTag is increased each time the program is run
                            tagFile = open('tagTracker.txt','w+') 
                            if numEnhancers > 1 and numEnhancers < totalBs and \
                                    bsT1 > 0 and bsT2 > 0:
                                    # we increase by 2 since the sort and pepper
                                    # generates two model files
                                    tagFile.write(str(fileTag + 2))
                            else:
                                    tagFile.write(str(fileTag + 1))

                            tagFile.close()

    else:  
        # DupBool = True section: bsTi is
        # binding sites per enhancer for Ti
        dupBool = True
        # iterating over 0 <= bsT1,bsT2 <= 4 (1 <= bsTotal <= 4), 1 <= numEnhancers <= 4
        # and also generating an Adair version of every model starting at +200
        for bsT1 in range(5):
            for bsT2 in range(5):
                totalBs = bsT1 + bsT2
                if totalBs < 5 and totalBs > 0:
                    # the upper bound of enhancers depends on the 
                    # total number of TF binding sites to keep
                    # the number of reactions at bay. 1 total bs
                    # can go up to 4 enhancers, 2 to 3, 3 to 2, ... 
                    enhancersUpper = 5 - totalBs + 1
                    for numEnhancers in range(1,enhancersUpper):
                        # doing 1 Adair per model
                        for multiAdairInd in range(2):
                            # reading the current file tag
                            tagFile = open('tagTracker.txt','r') 
                            fileTag = int(tagFile.read())
                            tagFile.close()

                            if multiAdairInd == 1:
                                adairBool = True
                                tempFileTag = fileTag + 200 
                            else:
                                adairBool = False
                                tempFileTag = fileTag                                  

                            if adairBool:
                                extraModel = "DupAdair" 
                            else:
                                # directory path 
                                extraModel = "Dup" 

                            FileGenerator(dupBool,bsT1,bsT2,numEnhancers,tempFileTag,\
                                        adairBool,multiBinding, extraModel)

                        # fileTag is increased each time the program is run
                        tagFile = open('tagTracker.txt','w+') 
                        if numEnhancers > 1 and numEnhancers < totalBs and \
                                bsT1 > 0 and bsT2 > 0:
                                # we increase by 2 since the sort and pepper
                                # generates two model files
                                tagFile.write(str(fileTag + 2))
                        else:
                                tagFile.write(str(fileTag + 1))

                        tagFile.close()

