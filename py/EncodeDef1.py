from Settings import *
from Helpers import *
import re
import os

def closeFiles():
    try:
        listFile.close()
        if MySetting.level>1:
            resultFile.close()
        inputFile.close()
        targetFile.close()
    except Exception:
        pass

def closeFile(myfile):
    try:
        myfile.close()
    except Exception:
        pass

def encode():
    level = MySetting.level
    #filename = MySetting.getFilename()
    proteinCount = 0
    residueCount = 0
#    newTagFlag = 0
#    nohomolgFlag = 0
#    noindentityFlag = 0
#    if level == 1:
#        newTagFlag = 1

    #101
    if level == 1 and MySetting.method.val==EncodingMethod.scr:
        readStats()
    elif level == 3:
        #102
        readStatsL3()

    FETURS = MySetting.getNoFeatures()    #104-111

    #code line 115-139 was ignored.
    #No nohomolgflag and noidentityflag processing

    #140
    #Prepare output file
    try:
        listFile = open(MySetting.listFileName)
        if MySetting.level>1:
            resultFile = open(MySetting.getResultFilename())
        inputFile = open(MySetting.getInputFilename(), 'w')
        targetFile = open(MySetting.getTargetFilename(), 'w')
    except IOError as e:
        print e
        try:
            closeFiles()
        except Exception:
            pass
        exit()

    listParser = re.compile(r"(\S+)")
    listFile.readline()
    count = 0
    faFileNameBase = MySetting.inputDir+"/fasta/"
    tgtFileNameBase = MySetting.inputDir+"/"+ str(MySetting.encodingMode)+"/"
    if MySetting.numOfstates ==2:
        tgtFileSuffix = ".sa25"
    elif MySetting.numOfstates == 3:
        tgtFileSuffix = ".ss"
    else:
        tgtFileSuffix = ".ss8"
    pssmFileNameBase = MySetting.inputDir+"/pssm/"
    pssm = {}

    #149
    for listLine in listFile:
        if count>=MySetting.MAXproteins:
            break
        id = listParser.match(listLine).group(0)
        pid = id[0:4]
        if len(id)==4:
            cid = '_'
        else:
            cid = id[4]
        id5 = pid+cid  # 5-letter id

        #process fasta sequence
        faFileName = faFileNameBase+id5+".fa"
        try:
            if os.path.getsize(faFileName) ==0:
                raise WindowsError("Fasta file "+faFileName+" is empty!")
            faFile = open(faFileName)
            faFile.readline()
            fasta = faFile.readline().strip()
            if len(fasta)<40 or len(fasta)>1000 or fasta.find('X') >= 0:
                print "Skip ", id5
                continue # if sequence is too short, too long or include X, discard them
        except WindowsError as e:
            print "Error: "+str(e)
        finally:
            try:
                faFile.close()
            except Exception:
                pass
        #end processing fasta sequence
        count+=1
        #process target file
        tgtFileName = tgtFileNameBase+id5+tgtFileSuffix
        try:
            if os.path.getsize(tgtFileName) ==0:
                raise WindowsError("Target file "+tgtFileName+" is empty!")
            tgtFile = open(tgtFileName)
            target = tgtFile.readline().strip()
            if MySetting.encodingMode == EncodingMode.ss:
                target = target.replace('X', 'C')
        except WindowsError as e:
            print "Error: "+str(e)
        finally:
            try:
                tgtFile.close()
            except Exception:
                pass
        #end processing target file
        print str(count), "Processing ", id5
        #200 process pssm files
        if level == 1 or level == 3:
            pssmFileName = pssmFileNameBase+id5+".pssm"
        try:
            if os.path.getsize(pssmFileName) ==0:
                raise WindowsError("PSSM file "+pssmFileName+" is empty!")
            pssmFile = open(pssmFileName)
            fillPssm(pssmFile, pssm)
        except WindowsError as e:
            print "Error: "+str(e)
        finally:
            try:
                pssmFile.close()
            except Exception:
                pass
        #end processing fasta sequence

        #201: read the result file
        if level>1:
            pp = []
            ppCount = 0
            for tr_a in target:
                if tr_a != '.':
                    ppItem = {}
                    trline = resultFile.readline().strip()
                    trlist = trline.split()
                    trIndex = 0
                    for itr in MySetting.tr:
                        ppItem[itr] = float(trlist[trIndex])
                        trIndex+=1
                    pp.append(ppItem)
#                else:
#                    for itr in MySetting.tr:
#                        if itr=='C':
#                            ppItem[itr] = 1.0
#                        else:
#                            pp[itr] = 0.0


        for i in range(0, len(fasta)): #223
            inputStr = ""
            outputRes = target[i]
            if outputRes==".":
                continue
            for pos in range(i-MySetting.halfWINsize, i+MySetting.halfWINsize+1):
                if pos<0 or pos>=len(fasta):
                    for iFETURS in range(1,FETURS):
                        #inputStr += "0.00 "
                        inputFile.write("0.00 ")
                    inputFile.write("1.00 ")
                else:
                    if level == 1 or level == 3:
                        inputFile.write(" ".join(map(formatFloat, pssm[pos])))
                        inputFile.write(" ")
                    if level == 1:
                        encodeScr(inputFile, pos, fasta)
                    else: #Copy results
                        inputFile.write(" ".join(map(formatFloat, pp)))

            inputFile.write('\n')

    closeFiles()



