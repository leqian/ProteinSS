from Settings import *
import math
from operator import add, mul

def fillPssm(pssmFile, pssm):
    for i in range(0, 3):
        pssmFile.readline() #skip first three lines
    rIndex = 0  #index for residues
    for line in pssmFile:
        line = line.strip()
        if len(line) ==0:
            return
        pssmi = []
        llist = line.split()
        index = 2 #index for
        for res in MySetting.residues:
            pssmi.append(1/(1+math.exp(-int(llist[index]))))
            index+=1
        pssm[rIndex] = pssmi
        rIndex+=1

def tripProb(probList, freqList, pL, pM, pR):#this one does not apply threshold
    #this function set a 3x3x3 list based on the combination of three targets
    #pL, pM and pR are three triple element lists represents probablity of H, E, C
    index = 0
    if max(pL)>=MySetting.L3Threshold:
        for i in range(0,3):
            if pL[i]>=MySetting.L3Threshold:
                break
        index=i*(MySetting.numOfStates*MySetting.numOfStates)
        if max(pR)>=MySetting.L3Threshold:
            for iii in range(0,3):
                if pR[iii]>=MySetting.L3Threshold:
                    break
            index += iii
            for ii in range(0,3):
                probList[ii] += freqList[index]
                index += MySetting.numOfStates
        else:
            for ii in range(0,3):
                for iii in range(0,3):
                    probList[ii] += freqList[index]
                    index += 1
        lSum = sum(probList)
        if lSum>0:
            for i in range(0,3):
                probList[i] /= lSum

    else:
        if max(pR)>=MySetting.L3Threshold:
            for iii in range(0,3):
                if pR[iii]>=MySetting.L3Threshold:
                    break
            index += iii
            for i in range(0,3):
                for ii in range(0,3):
                    probList[ii] += freqList[index]
                    index += MySetting.numOfStates
            lSum = sum(probList)
            if lSum>0:
                for i in range(0,3):
                    probList[i] /= lSum
        else:
            for i in range(0,3):
                for ii in range(0,3):
                    for iii in range(0,3):
                        probList[ii] += freqList[index]
                        index += 1

def tripProbB(probList, freqList, pL, pM, pR):#this one does not apply threshold
    #this function set a 3x3x3 list based on the combination of three targets
    #pL, pM and pR are three triple element lists represents probablity of H, E, C
    index = 0
    for i in range(0,3):
        for ii in range(0,3):
            for iii in range(0,3):
                probList[ii] += pL[i]*pM[ii]*pR[iii]*freqList[index]
                index += 1

def averageList(llist):
    s = sum(llist)
    if s==0:
        return llist
    else:
        return map(lambda x:x/s, llist)

def formatFloat(number):
    return "%5.3f" % number

def encodeScr(inputFile, pos, fasta):
    #write scores into input file
    singletScr = MySetting.singletScr
    doubletScr = MySetting.doubletScr
    tripletScr = MySetting.tripletScr

    energy_d = [0.0]*MySetting.numOfStates
    energy_t = [0.0]*MySetting.numOfStates

    #process middle score
    currRes = fasta[pos]
    for tripType in range(0, 9):
        t = tripType/3 +1
        tt = tripType%3 + 1
        if pos-t>=0 and pos+tt<len(fasta):
            prev_aa = fasta[pos-t]
            next_aa = fasta[pos+tt]
            res = prev_aa+currRes+next_aa;
            energy_t = map(add, energy_t, tripletScr.values(tripType, TripletScr.mid, res))

    #process left and right score
    for tripType in range(0,3):
        t = 1 if tripType <2 else 2
        tt = 3 if tripType else 2
        #left
        if  pos+tt<len(fasta):
            next_aa =fasta[pos+t]
            nextnext_aa=fasta[pos+tt]
            res=currRes+next_aa+nextnext_aa
            energy_t = map(add, energy_t, tripletScr.values(tripType,TripletScr.left, res))

        #right
        if  pos-tt>=0:
            next_aa =fasta[pos-tt+t]
            nextnext_aa=fasta[pos-tt]
            res=nextnext_aa+next_aa+currRes
            energy_t = map(add, energy_t, tripletScr.values(tripType,TripletScr.right, res))

    #process doublet
    for doubleType in range(0, 3):
        t = doubleType+1
        #right
        if pos-t>=0:
            next_aa = fasta[pos-t]
            res = next_aa+currRes
            energy_d = map(add, energy_d, doubletScr.values(doubleType, DoubletScr.right, res))

        #left
        if pos+t<len(fasta):
            next_aa = fasta[pos+t]
            res = currRes+next_aa
            energy_d = map(add, energy_d, doubletScr.values(doubleType, DoubletScr.left, res))

    #process singlet
    energy_s = singletScr.values(currRes)
    energy_all = map(add, map(add, energy_t, energy_s), [-4*item for item in energy_d])


    #final processing
    if MySetting.encodingMode==EncodingMode.ss:
        if pos==0 or pos == len(fasta)-1:
            energy_all[MySetting.numOfStates-1] -=2;
    inputFile.write(' '.join(map(formatFloat, energy_all)))



def readStats():
    print "Reading singlets scores"
    if MySetting.numOfStates==3:
        statsDir = MySetting.statDir+MySetting.statDir3Suffix
    elif MySetting.numOfStates==8:
        statsDir = MySetting.statDir+MySetting.statDir8Suffix
    else:
        statsDir = MySetting.statDir+MySetting.statDir2Suffix

    singletScr = MySetting.singletScr
    doubletScr = MySetting.doubletScr
    tripletScr = MySetting.tripletScr

    singletFileName = statsDir+"0.score"
    try:
        singletFile = open(singletFileName)
        #sscomb = (singletFile.readline().strip().split())[1:]
        singletFile.readline()
        for line in singletFile:
            llist = line.strip().split()
            if not llist:
                break
            singletScr.setValues(llist[0], map(float, llist[1:]))
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            singletFile.close()
        except Exception:
            pass
    ###################### get doublets
    print "Reading doublets scores"
    try:
        for doubleType in range(0, 3):
            t=doubleType+1
            doubletFileName = statsDir+"0_"+str(t)+".score"
            doubletFile = open(doubletFileName)
            #sscomb = (singletFile.readline().strip().split())[1:]
            doubletFile.readline()
            for line in doubletFile:
                llist = line.strip().split()
                if not llist:
                    break
                doubletScr.setValues(doubleType, DoubletScr.left, llist[0],map(float, llist[1:(MySetting.numOfStates+1)])) #left
                doubletScr.setValues(doubleType, DoubletScr.right, llist[0],map(float, llist[(MySetting.numOfStates+1):])) #right
            doubletFile.close()
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            doubletFile.close()
        except Exception:
            pass

    # get triplet scores
    print "Reading triplets scores"
    try:
        ###################### get mid triplets
        for index in range(0, 9):
            t = index/3 +1
            tt = index%3 + t +1
            tripletFileName = statsDir+"0_"+str(t)+"_"+str(tt)+"_mid.score"
            tripletFile = open(tripletFileName)
            #sscomb = (singletFile.readline().strip().split())[1:]
            tripletFile.readline()
            for line in tripletFile:
                llist = line.strip().split()
                if not llist:
                    break
                tripletScr.setValues(index, TripletScr.mid, llist[0] ,map(float, llist[1:]))
            tripletFile.close()
        # Get left triplets
        for index in range(0, 3):
            t = 1 if index <2 else 2  #index = 0: (1, 2), 1: (1, 3), 2: (2, 3)
            tt = 3 if index else 2
            tripletFileName = statsDir+"0_"+str(t)+"_"+str(tt)+"_left.score"
            tripletFile = open(tripletFileName)
            #sscomb = (singletFile.readline().strip().split())[1:]
            tripletFile.readline()
            for line in tripletFile:
                llist = line.strip().split()
                if not llist:
                    break
                tripletScr.setValues(index, TripletScr.left, llist[0],map(float, llist[1:]))
            tripletFile.close()
        # Get right triplets
        for index in range(0, 3):
            t = 1 if index <2 else 2
            tt = 3 if index else 2
            tripletFileName = statsDir+"0_"+str(t)+"_"+str(tt)+"_right.score"
            tripletFile = open(tripletFileName)
            #sscomb = (singletFile.readline().strip().split())[1:]
            tripletFile.readline()
            for line in tripletFile:
                llist = line.strip().split()
                if not llist:
                    break
                tripletScr.setValues(index, TripletScr.right, llist[0], map(float, llist[1:]))
            tripletFile.close()
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            tripletFile.close()
        except Exception:
            pass


def readStatsL3():
    print "Reading singlets freqs"
    if MySetting.numOfStates==3:
        freqDir = MySetting.statDir+MySetting.freqDir3Suffix
    elif MySetting.numOfStates==8:
        freqDir = MySetting.statDir+MySetting.freqDir8Suffix

    singletFreq = MySetting.singletFreq
    doubletFreq = MySetting.doubletFreq
    tripletFreq = MySetting.tripletFreq

    singletFileName = freqDir+"0.freq"
    try:
        singletFile = open(singletFileName)
        singletFile.readline() #skip the first line
        for line in singletFile:
            llist = line.strip().split()
            if not llist:
                break
            singletFreq.setValues(llist[0], averageList(map(float, llist[1:])))
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            singletFile.close()
        except Exception:
            pass

    '''
    ###################### get doublets
    print "Reading doublets freqs"
    try:
        for doubleType in range(0, 3):
            t=doubleType+1
            doubletFileName = statsDir+"0_"+str(t)+".freq"
            doubletFile = open(doubletFileName)
            doubletFile.readline()
            for line in doubletFile:
                llist = line.strip().split()
                if not llist:
                    break
                doubletFreq.setValues(doubleType, DoubletFreq.left, llist[0], averageList(map(float, llist[1:]))) #left
            doubletFile.close()
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            doubletFile.close()
        except Exception:
            pass
    '''
    # get triplet scores
    print "Reading triplets Freqs"
    try:
        ###################### get mid triplets
        for index in range(0, 3):
            t = index+1
            tt = t*2
            tripletFileName = freqDir+"0_"+str(t)+"_"+str(tt)+".freq"
            tripletFile = open(tripletFileName)
            #sscomb = (singletFile.readline().strip().split())[1:]
            tripletFile.readline()
            for line in tripletFile:
                llist = line.strip().split()
                if not llist:
                    break
                tripletFreq.setValues(index, TripletFreq.mid, llist[0] , averageList(map(float, llist[1:])))
            tripletFile.close()
    except Exception as e:
        errorLog(e)
        exit()
    finally:
        try:
            tripletFile.close()
        except Exception:
            pass


def encodeScrL3b(inputFile, pos, fasta, pp):
    #inputFile: file to write output
    #pos: current residue position
    #fasta: string of the protein
    #pp: probability output from previous level of training.

    #write score
    inputFile.write(" ".join(map(formatFloat, pp[pos])))
    inputFile.write(" ")

    '''
    #check threshold
    pM = pp[pos]  #probability of the middle residue
    if  (pM[0]>=MySetting.L3Threshold):
        inputFile.write('1 0 0 '*3)
        return
    elif (scores[1]>=MySetting.L3Threshold):
        inputFile.write('0 1 0 '*3)
        return
    elif (scores[0]>=MySetting.L3Threshold):
        inputFile.write('0 0 1 '*3)
        return
    '''

    #process score
    currRes = fasta[pos]
    for tripType in range(0, 3):
        t = tripType+1
        if pos-t>=0 and pos+t<len(fasta):
            freqs = [0.0]*(MySetting.numOfStates)
            prev_aa = fasta[pos-t]
            next_aa = fasta[pos+t]
            res = prev_aa+currRes+next_aa;
            tripProb(freqs, MySetting.tripletFreq.values(tripType, TripletFreq.mid, res), pp[pos-t], pp[pos], pp[pos+t])
        else: #end points
            freqs=map(mul, MySetting.singletFreq.values(currRes), pp[pos])
        inputFile.write(' '.join(map(formatFloat, freqs)))
        inputFile.write(' ')

def encodeScrL3Tradition(inputFile, pos, fasta, pp):
    #inputFile: file to write output
    #pos: current residue position
    #fasta: string of the protein
    #pp: probability output from previous level of training.

    #check threshold
    pM = pp[pos]  #probability of the middle residue
    if  (pM[0]>=MySetting.L3Threshold):
        inputFile.write('1 0 0 '*4)
        return
    elif (pM[1]>=MySetting.L3Threshold):
        inputFile.write('0 1 0 '*4)
        return
    elif (pM[2]>=MySetting.L3Threshold):
        inputFile.write('0 0 1 '*4)
        return

    #write score
    inputFile.write(" ".join(map(formatFloat, pp[pos])))
    inputFile.write(" ")

    #process score
    currRes = fasta[pos]
    for tripType in range(0, 3):
        t = tripType+1
        if pos-t>=0 and pos+t<len(fasta):
            freqs = [0.0]*(MySetting.numOfStates)
            prev_aa = fasta[pos-t]
            next_aa = fasta[pos+t]
            res = prev_aa+currRes+next_aa;
            tripProb(freqs, MySetting.tripletFreq.values(tripType, TripletFreq.mid, res), pp[pos-t], pp[pos], pp[pos+t])
        else: #end points
            freqs=MySetting.singletFreq.values(currRes)
        inputFile.write(' '.join(map(formatFloat, freqs)))
        inputFile.write(' ')

def encodeScrL3(inputFile, pos, fasta, pp):
    #inputFile: file to write output
    #pos: current residue position
    #fasta: string of the protein
    #pp: probability output from previous level of training.

    #check threshold
    pM = pp[pos]  #probability of the middle residue

    #write score
    inputFile.write(" ".join(map(formatFloat, pp[pos])))
    inputFile.write(" ")

    #process score
    currRes = fasta[pos]
    for tripType in range(0, 3):
        t = tripType+1
        if pos-t>=0 and pos+t<len(fasta):
            freqs = [0.0]*(MySetting.numOfStates)
            prev_aa = fasta[pos-t]
            next_aa = fasta[pos+t]
            res = prev_aa+currRes+next_aa;
            tripProb(freqs, MySetting.tripletFreq.values(tripType, TripletFreq.mid, res), pp[pos-t], pp[pos], pp[pos+t])
        else: #end points
            freqs=MySetting.singletFreq.values(currRes)
        inputFile.write(' '.join(map(formatFloat, freqs)))
        inputFile.write(' ')

def errorLog(str):
    print str
    try:
        MySetting.errorLogFile.write(str+'\n')
    except:
        pass
