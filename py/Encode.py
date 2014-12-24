from Settings import *
from EncodeDef import *
from DirSettings import *
import sys

def main():
    #Parse arguments
    try:
        if len(sys.argv)<2:
            raise ValueError("Too Few arguments!")
        MySetting.numOfStates = int(sys.argv[1])
        if len(sys.argv)>2:
            MySetting.MAXproteins = int(sys.argv[2])
        else:
            MySetting.MAXproteins = 0
        if MySetting.MAXproteins<=10 and MySetting.MAXproteins!=0:
            raise ValueError("Max Proteins is too small (must be > 10)")
        if len(sys.argv)>3:
            MySetting.level = int(sys.argv[3]);
        else:
            MySetting.level = 1
        if MySetting.level<1 or MySetting.level>3:
            raise ValueError("Level must be 1, 2 or 3")
        if len(sys.argv)>4:
            MySetting.offset = int(sys.argv[4])
        else:
            MySetting.offset = 0
        if MySetting.offset<0:
            raise ValueError("Offset must be non-negative")
    except ValueError as detail:
        print detail
        print "Usage: Python encode.py numOfStates[2/3/8] [MAXProteins] [Level[1/2/3]] [Offset]"
        exit()

    #Initiallization:
    setDirectory()
    if MySetting.numOfStates == 2:
        MySetting.encodingMode = EncodingMode.sa
    else:
        MySetting.encodingMode = EncodingMode.ss


    if MySetting.numOfStates == 2:
        MySetting.tr = MySetting.sa2TR
    elif MySetting.numOfStates == 3:
        MySetting.tr = MySetting.ss3TR
    else:
        MySetting.tr = MySetting.ss8TR

    pos = 0
    for c in MySetting.tr:
        cc = ''
        for i in range(0, MySetting.numOfStates):
            if i==pos:
                cc+="1 "
            else:
                cc += '0 '
        cc += '\n'
        pos += 1
        MySetting.TGTCode[c] = cc

    if MySetting.numOfStates == 3:
        for r in MySetting.ss3TR:
            for rr in MySetting.ss3TR:
                MySetting.dTr.append(r+rr)

        for r in MySetting.ss3TR:
            for rr in MySetting.ss3TR:
                for rrr in MySetting.ss3TR:
                    MySetting.tTr.append(r+rr+rrr)

    encode()

if __name__ == '__main__':
    main()
