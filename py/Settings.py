#constants
#from enum import Enum

class EncodingMethod:
    scr = 1
    pssm = 0
    @staticmethod
    def toStr(val):
        if val==EncodingMethod.scr:
            return "scr"
        elif val==EncodingMethod.pssm:
            return "pssm"

class EncodingMode:
    ss = 0
    sa = 1
    @staticmethod
    def toStr(val):
        if val==EncodingMode.ss:
            return "ss"
        elif val==EncodingMode.sa:
            return "sa"

class SingletScr:
    #format of singleScr:
        #singletScr['A'] = [score_of_H,...]
    def __init__(self):
        self.scrs = {}

    def values(self, res):
        try:
            return self.scrs[res]
        except:
            return [0.0]*MySetting.numOfStates

    def value(self, res, idx):
        try:
            return self.scrs[res][idx]
        except:
            return 0.0

    def setValues(self, res, values):
        self.scrs[res] = values

class DoubletScr:
        #format of doubletScr:
        #doubletScr=[doublet_offset_1_scr, ..., doublet_offset_3_scr]
        #doublet_offset_i = [doublet_offset_i_left_scr, doublet_offset_i_right_scr]
        #doublet_offset_i_left_scr = {"AA": [score_of_H, ...], "AC":[....]}
        #doublet_offset_i_right_scr = {"AA": [score_of_H, ...], "AC":[....]}
    left = 0
    right = 1
    def __init__(self):
        self.scrs = []
        for i in range(0,3):
            left_scr = {}
            right_scr = {}
            scr = [left_scr, right_scr]
            self.scrs.append(scr)

    def value(self, offset, lr, res, idx):
        try:
            return self.scrs[offset][lr][res][idx]
        except:
            return 0.0

    def values(self, offset, lr, res):
        try:
            return self.scrs[offset][lr][res]
        except:
            return [0.0]*MySetting.numOfStates

    def setValues(self, offset, lr, res, values):
        self.scrs[offset][lr][res] = values

class TripletScr:
    #format of tripletScr:
    #tripletScr = [TripletMidScr, TripetLeftScr, TripletRightScr]
    #TripletMidScr = [tripletMidScr_0, ...,tripletMidScr_8]
    # index/3+1 is t
    # index%3+index/3+2 is tt
    #for example: 0-> 0-1-2. 5->0-2-5
    #tripletMidScr_i = {"AAA":[score_of_H, score_of_E, score_of_C],...} for 3-state
    #tripletLeftScr = [TripletLeftScr_0, ..., TripletLeftScr_2]
    # i = 0->(1,2), 1->(1,3), 2->(2,3)
    # TripletLeftScr_i ={"AAA":[score_of_H, score_of_E, score_of_C],...} for 3-state
    # offset = 0~8 (mid) or 0~3 (left or right).
    #  res =
    left = 0
    right = 1
    mid = 2
    def __init__(self):
        midscrs = []
        leftscrs = []
        rightscrs = []
        self.scrs = [leftscrs, rightscrs, midscrs]
        for i in range(0,9):
            scri = {}
            midscrs.append(scri)
        for i in range(0, 3):
            scrli = {}
            leftscrs.append(scrli)
            scrri = {}
            rightscrs.append(scrri)

    def values(self, offset, lmr, res):
        try:
            return self.scrs[lmr][offset][res]
        except:
            return [0.0]*MySetting.numOfStates

    def value(self, offset, lmr, res, idx):
        try:
            return self.scrs[lmr][offset][res][idx]
        except:
            return 0.0

    def setValues(self, offset, lmr, res, values):
        self.scrs[lmr][offset][res] = values

class SingletFreq:
    #format of singleScr:
        #singletScr['A'] = {seq:value ...}
    def __init__(self):
        self.scrs = {}

    def value(self, seq):
        try:
            return self.scrs[seq]
        except:
            return 0.0

    def values(self, seq):
        try:
            return self.scrs[seq]
        except:
            return [0.0]*MySetting.numOfStates

    def setValues(self, res, values):
        self.scrs[res] = values

class DoubletFreq:
        #format of doubletScr:
        #doubletScr=[doublet_offset_1_scr, ..., doublet_offset_3_scr]
        #doublet_offset_i = [doublet_offset_i_left_scr, doublet_offset_i_right_scr]
        #doublet_offset_i_left_scr = {"AA": [score_of_H, ...], "AC":[....]}
        #doublet_offset_i_right_scr = {"AA": [score_of_H, ...], "AC":[....]}
    left = 0
    right = 1
    def __init__(self):
        self.scrs = []
        for i in range(0,3):
            left_scr = {}
            right_scr = {}
            scr = [left_scr, right_scr]
            self.scrs.append(scr)

    def value(self, offset, lr, seq):
        try:
            return self.scrs[offset][lr][seq]
        except:
            return 0.0

    def setValues(self, offset, lr, res, values):
        for i in range(0, len(values)):
            self.scrs[offset][lr][res+MySetting.dTr[i]] = values[i]

class TripletFreq:
    #format of tripletScr:
    #tripletScr = [TripletMidScr, TripetLeftScr, TripletRightScr]
    #TripletMidScr = [tripletMidScr_0, ...,tripletMidScr_8]
    # index/3+1 is t
    # index%3+index/3+2 is tt
    #for example: 0-> 0-1-2. 5->0-2-5
    #tripletMidScr_i = {"AAA":[score_of_H, score_of_E, score_of_C],...} for 3-state
    #tripletLeftScr = [TripletLeftScr_0, ..., TripletLeftScr_2]
    # i = 0->(1,2), 1->(1,3), 2->(2,3)
    # TripletLeftScr_i ={"AAA":[score_of_H, score_of_E, score_of_C],...} for 3-state
    # offset = 0~8 (mid) or 0~3 (left or right).
    #  res =
    left = 0
    right = 1
    mid = 2
    def __init__(self):
        midscrs = []
        leftscrs = []
        rightscrs = []
        self.scrs = [leftscrs, rightscrs, midscrs]
        for i in range(0,9):
            scri = {}
            midscrs.append(scri)
        for i in range(0, 3):
            scrli = {}
            leftscrs.append(scrli)
            scrri = {}
            rightscrs.append(scrri)

    #def value(self, offset, lmr, seq):
    #    try:
    #        return self.scrs[lmr][offset][seq]
    #    except:
    #        return 0.0

    def values(self, offset, lmr, res):
        try:
            return self.scrs[lmr][offset][res]
        except:
            return [0.0]*(MySetting.numOfStates*MySetting.numOfStates*MySetting.numOfStates)

    def setValues(self, offset, lmr, res, values):
        self.scrs[lmr][offset][res] = values

class MySetting:
    inputDir = ""  #inputs
    outputDir = ""  #output
    resultDir = ""  #trained result
    statDir = ""  #conext-based score
    statDir3Suffix = ""
    statDir8Suffix = ""
    freqDir3Suffix = ""
    freqDir8Suffix = ""
    statDir2Suffix = ""
    residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    sa2TR = ['B', 'E']
    ss3TR = ['H','E','C']
    ss8TR = ['H','G','I','E','B','T','S','C']
    TGTCode = {}
    tr = []
    dTr = []
    tTr = []
    halfWINsize	= 7
    halfWINsizeL3	= 2
    offset = 0
    L3Threshold = 0.9
    WINsize	= halfWINsize*2+1
    #fileList =
    numOfStates =3 	#3, 8 ,2
    encodingMode = EncodingMode.ss if numOfStates>2 else EncodingMode.sa
    method = EncodingMethod.scr 		#scr = 1 or pssm = 0
    level = 1
    limited = 0
    MAXproteins = 1
    noPSSM = 20 #number of PSSM
    # argv>4
    listFileName = ""
    errorLogFileName = ""
    errorLogFile = ""
    title = 'cull'
    singletScr = SingletScr()
    doubletScr = DoubletScr()
    tripletScr = TripletScr()
    singletFreq = SingletFreq()
    doubletFreq = DoubletFreq()
    tripletFreq = TripletFreq()

    @classmethod
    def getFilename(cls):
        base = ""+EncodingMode.toStr(cls.encodingMode)+str(cls.numOfStates)+"_"+EncodingMethod.toStr(cls.method)+"M"+str(cls.MAXproteins)
        if cls.offset>0:
            base+="F"+str(cls.offset)
        return base

    @classmethod
    def getNoFeatures(cls):
        if cls.level == 1:
            return cls.noPSSM+cls.numOfStates+1
        elif cls.level == 2:
            return cls.numOfStates+1
        else:
            return cls.noPSSM+4*cls.numOfStates+1
            #return cls.noPSSM+4*cls.numOfStates+1

    @classmethod
    def getNoFeaturesL3(cls):
        return 4*cls.numOfStates+1
        #return 4*cls.numOfStates+1

    @classmethod
    def getResultFilename(cls):
        return cls.resultDir+"/output"+str(cls.level)+"_"+cls.getFilename()+".dat"

    @classmethod
    def getInputFilename(cls):
        return cls.outputDir+"/input"+str(cls.level)+"_"+cls.getFilename()+".dat"

    @classmethod
    def getTargetFilename(cls):
        return cls.outputDir+"/target"+"_"+cls.getFilename()+".dat"

