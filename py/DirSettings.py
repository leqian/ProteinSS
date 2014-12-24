from Settings import *
from datetime import datetime

def setDirectory():
    MySetting.inputDir = "L:/odu/Cull_all"  #inputs
    MySetting.outputDir = "C:/odu/Output"  #output
    MySetting.resultDir = "C:/odu/Output"  #trained result
    MySetting.statDir = "L:/odu/TripletFreq"  #conext-based score
    MySetting.statDir3Suffix = "/SS3_stat/cull16633freq_ss3/"
    MySetting.freqDir3Suffix = "/SS3_stat/cull16633freq_ss3x/"
    MySetting.statDir8Suffix = "/SS8_stat/cull16633freq_ss8/"
    MySetting.freqDir8Suffix = "/SS3_stat/cull16633freq_ss8x/"
    MySetting.statDir2Suffix = "/SA_stat/cull16633freq_sa/"
    MySetting.listFileName = "E:/gdrive/Documents/Research/ODU/Lists/cullpdb_pc25_res3.0_R1.0_d120428_chains7987"
    MySetting.errorLogFileName = MySetting.outputDir+'/errorlog.txt'
    try:
        MySetting.errorLogFile = file(MySetting.errorLogFile, 'a')
        MySetting.errorLogFile.write('\n')
        MySetting.errorLogFile.write(str(datetime.now()))
        MySetting.errorLogFile.write('\n')
    except:
        print('Cannot open errorlog file '+MySetting.errorLogFileName)

