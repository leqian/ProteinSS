from Settings import *
from Helpers import *
import re
import os

#targetFile is for the output
#tgtFile is input result.

def closeFiles():
    try:
        listFile.close()
    except Exception:
        pass

    if MySetting.level>1:
        try:
            resultFile.close()
        except Exception:
            pass
    try:
       inputFile.close()
    except Exception:
        pass

    if MySetting.level==1:
        try:
            targetFile.close()
        except Exception:
            pass

    try:
        MySetting.errorLogFile.close()
    except:
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

    trEmpty = [0.0]*MySetting.numOfStates
    trEmpty[-1] = 1.0
    #101
    if level == 1 and MySetting.method==EncodingMethod.scr:
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
        if MySetting.level==1:
            targetFile = open(MySetting.getTargetFilename(), 'w')
    except IOError as e:
        errorLog(e)
        try:
            closeFiles()
        except Exception:
            pass
        exit()

    listParser = re.compile(r"(\S+)")
    listFile.readline()
    count = 0
    faFileNameBase = MySetting.inputDir+"/fasta/"
    tgtFileNameBase = MySetting.inputDir+"/"+ EncodingMode.toStr(MySetting.encodingMode)+"/"
    if MySetting.numOfStates ==2:
        tgtFileSuffix = ".sa25"
    elif MySetting.numOfStates == 3:
        tgtFileSuffix = ".ss"
    else:
        tgtFileSuffix = ".ss8"
    pssmFileNameBase = MySetting.inputDir+"/pssm/"
    pssm = {}
    ppEmpty = [0.0]*(MySetting.numOfStates-1)
    ppEmpty.append(1.0)

    #149
    offsetCount = 0
    for listLine in listFile:
        if offsetCount<MySetting.offset:
            offsetCount+=1;
            continue
        if (count>=MySetting.MAXproteins) and (MySetting.MAXproteins>0):
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
                errorLog("Skip "+str(id5))
                continue # if sequence is too short, too long or include X, discard them
        except WindowsError as e:
            errorLog("Error: "+str(e))
        finally:
            try:
                faFile.close()
            except Exception:
                pass
        #end processing fasta sequence
        count+=1
        #process tgt file
        tgtFileName = tgtFileNameBase+id5+tgtFileSuffix
        try:
            if os.path.getsize(tgtFileName) ==0:
                raise WindowsError("Target file "+tgtFileName+" is empty!")
            tgtFile = open(tgtFileName)
            tgt = tgtFile.readline().strip()
            if MySetting.encodingMode == EncodingMode.ss:
                tgt = tgt.replace('X', 'C')
        except WindowsError as e:
            errorLog("Error: "+str(e))
            continue
        finally:
            try:
                tgtFile.close()
            except Exception:
                pass
        #end processing target file
        errorLog(str(count)+ "Processing "+str(id5))
        #200 process pssm files
        if level == 1 or level == 3:
            pssmFileName = pssmFileNameBase+id5+".pssm"
            try:
                if os.path.getsize(pssmFileName) ==0:
                    raise WindowsError("PSSM file "+pssmFileName+" is empty!")
                pssmFile = open(pssmFileName)
                fillPssm(pssmFile, pssm)
            except WindowsError as e:
                errorLog("Error: "+str(e))
            finally:
                try:
                    pssmFile.close()
                except Exception:
                    pass
        #end processing fasta sequence

        #201: read the result file
        #structure of pp: pp = [tr0_Prob, tr1_prob, ...]
        # tri_prob = [p0, p1, ...]
        if level>1:
            pp = []
            ppCount = 0
            for tr_a in tgt:
                if tr_a != '.':
                    trline = resultFile.readline().strip()
                    trlist = trline.split()
                    pp.append(map(float, trlist))
                else:
                    pp.append(trEmpty)

        #223: encode
        terminalStr = '0 '*(FETURS-1)
        terminalStr+='1 '
        terminalStrL3 = '0 '*(MySetting.getNoFeaturesL3()-1)
        terminalStrL3+='1 '

        for i in range(0, len(fasta)):
            outputRes = tgt[i]
            if outputRes==".":
                continue
            for pos in range(i-MySetting.halfWINsize, i+MySetting.halfWINsize+1):
                if pos<0 or pos>=len(fasta):
                    if level==3 and abs(i-pos)>MySetting.halfWINsizeL3:
                        inputFile.write(terminalStrL3)  #terminal
                    else:
                        inputFile.write(terminalStr)
                else:
                    if level == 1 or (level == 3 and abs(i-pos)<=MySetting.halfWINsizeL3) :
                        inputFile.write(" ".join(map(formatFloat, pssm[pos])))
                        inputFile.write(" ")
                    if level == 1:
                        encodeScr(inputFile, pos, fasta)
                    if level==2:
                        inputFile.write(" ".join(map(formatFloat, pp[pos])))
                        inputFile.write(" ")
                    if level==3:
                        try:
                            encodeScrL3(inputFile, pos, fasta, pp)
                        except:
                            print i, ' ', pos
                    inputFile.write(' 0 ') #non-terminal

            inputFile.write('\n')
            #263: Write target file
            if MySetting.level==1:
                targetFile.write(MySetting.TGTCode[outputRes])

            residueCount+=1
        proteinCount+=1
    closeFiles()
    errorLog("Total Proteins: "+str(proteinCount))
    errorLog("Total Residues: "+str(residueCount))




