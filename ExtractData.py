#!/usr/bin/python

import os, sys, fnmatch

PlotDirectory = r'Plot'
SaveExtractedData = 'ExtractedData'

TemperatureList = ['0C','-20C', 'RT']
AnnealingSteps = ['0min', '20min', '40min', '80min', '156min', '312min', '624min', '1248min']
SensorNames = {# in increasing fluence order
               'FTH200N': ("FTH200N_25_DiodeL_8","FTH200N_25_DiodeL_3","FTH200N_25_DiodeL_2","FTH200N_24_DiodeL_8"), 
               'FTH200Y': ("FTH200Y_04_DiodeL_9","FTH200Y_03_DiodeL_9","FTH200Y_04_DiodeL_8","FTH200Y_04_DiodeL_3"),
               'MCZ200N': ("MCZ200N_06_DiodeL_9","MCZ200N_11_DiodeL_9","MCZ200N_09_DiodeL_8","MCZ200N_11_DiodeL_8"),
               'MCZ200Y': ("MCZ200Y_06_DiodeL_5","MCZ200Y_04_DiodeL_3","MCZ200Y_07_DiodeL_9","MCZ200Y_05_DiodeL_2"),
               }
Fluences = [1.7, 10.0, 30.0, 150.0] # in 1E13cm-2
Frequencies = ['455', '1000']


def findfiles(path, filtre):
    for root, dirs, files in os.walk(path):
        for f in fnmatch.filter(files, filtre):
            yield os.path.join(root, f)

CurrentUnitDict = {'(fA)': (1.e-9), '(pA)': (1.e-6), '(nA)': (1.e-3), '(uA)': (1.), '(mA)': (1.e3), '(A)': (1.e6)}

def readfile(textfile):
    f = open(textfile, 'rU')
    lines = f.readlines()
    f.close()

    i = 0
    for line in lines:
        i+=1
        data = line.strip().split()
        if line.startswith("V_dep:"):
            if i<18:
                Vdepf1log = float(data[1])
                Vdepf1logEr = float(data[3])
            if i>18 and i<35:
                Vdepf2log = float(data[1])
                Vdepf2logEr = float(data[3])
            if i>35 and i<52:
                Vdepf1IC = float(data[1])
                Vdepf1ICEr = float(data[3])
            if i>52:
                Vdepf2IC = float(data[1])
                Vdepf2ICEr = float(data[3])
        elif line.startswith("I_(V_dep)("):
            if i<18:
                if len(data)>3:
                    Idepf1log = float(data[1])*CurrentUnitDict[data[4]]
                    Idepf1logEr = float(data[3])*CurrentUnitDict[data[4]]
                else:
                    Idepf1log = float(data[1])*CurrentUnitDict[data[2]]
                    Idepf1logEr = -1
            if i>18 and i<35:
                if len(data)>3:
                    Idepf2log = float(data[1])*CurrentUnitDict[data[4]]
                    Idepf2logEr = float(data[3])*CurrentUnitDict[data[4]]
                else:
                    Idepf2log = float(data[1])*CurrentUnitDict[data[2]]
                    Idepf2logEr = -1
            if i>35 and i<52:
                if len(data)>3:
                    Idepf1IC = float(data[1])*CurrentUnitDict[data[4]]
                    Idepf1ICEr = float(data[3])*CurrentUnitDict[data[4]]
                else:
                    Idepf1IC = float(data[1])*CurrentUnitDict[data[2]]
                    Idepf1ICEr = -1
            if i>52:
                if len(data)>3:
                    Idepf2IC = float(data[1])*CurrentUnitDict[data[4]]
                    Idepf2ICEr = float(data[3])*CurrentUnitDict[data[4]]
                else:
                    Idepf2IC = float(data[1])*CurrentUnitDict[data[2]]
                    Idepf2ICEr = -1
        else:
            pass
    Liste = [Vdepf1log,Vdepf1logEr,Vdepf2log,Vdepf2logEr,Vdepf1IC,Vdepf1ICEr,Vdepf2IC,Vdepf2ICEr,Idepf1log,Idepf1logEr,Idepf2log,Idepf2logEr,Idepf1IC,Idepf1ICEr,Idepf2IC,Idepf2ICEr]
    return Liste

textfilelist = []
for textfile in findfiles(PlotDirectory, '*.dat'):
    textfilelist.append(textfile)

for item1 in SensorNames:
    SortedSensor = fnmatch.filter(textfilelist, '*'+item1+'*')
    for item2 in TemperatureList:
        SortedTemp = fnmatch.filter(SortedSensor, '*_'+item2+'*')
        #print item1, item2, SortedTemp
        with open(SaveExtractedData+'/'+item1+'_'+item2+'.txt','w') as fout:
            VdepAll = []
            VdepErAll = []
            IdepAll = []
            IdepErAll = []
            for item3 in AnnealingSteps:
                Sorted = fnmatch.filter(fnmatch.filter(SortedTemp, '*CONNECTED*'), '*_'+item3+'*')
                ExtractedData = {}
                for sensor in SensorNames[item1]:
                    sensortextfile = fnmatch.filter(Sorted, '*'+sensor+'*')
                    if sensortextfile == []:
                        ExtractedData[sensor] = [-1]*16
                    else:
                        ExtractedData[sensor] = readfile(sensortextfile[-1])
                Vdep = [-1]*(4*len(SensorNames[item1])+1)
                VdepEr = [-1]*(4*len(SensorNames[item1])+1)
                Idep = [-1]*(4*len(SensorNames[item1])+1)
                IdepEr = [-1]*(4*len(SensorNames[item1])+1)
                Vdep[0] = float(item3[:-3])
                VdepEr[0] = float(item3[:-3])
                Idep[0] = float(item3[:-3])
                IdepEr[0] = float(item3[:-3])
                for i in range(len(SensorNames[item1])):
                    Vdep[i+1] = ExtractedData[SensorNames[item1][i]][0]
                    Vdep[i+1+len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][4]
                    Vdep[i+1+2*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][2]
                    Vdep[i+1+3*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][6]

                    VdepEr[i+1] = ExtractedData[SensorNames[item1][i]][1]
                    VdepEr[i+1+len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][5]
                    VdepEr[i+1+2*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][3]
                    VdepEr[i+1+3*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][7]

                    Idep[i+1] = ExtractedData[SensorNames[item1][i]][8]
                    Idep[i+1+len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][12]
                    Idep[i+1+2*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][10]
                    Idep[i+1+3*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][14]

                    IdepEr[i+1] = ExtractedData[SensorNames[item1][i]][9]
                    IdepEr[i+1+len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][13]
                    IdepEr[i+1+2*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][11]
                    IdepEr[i+1+3*len(SensorNames[item1])] = ExtractedData[SensorNames[item1][i]][15]
                VdepAll.append(Vdep)
                VdepErAll.append(VdepEr)
                IdepAll.append(Idep)
                IdepErAll.append(IdepEr)
            fout.write("Sample: ")
            fout.write(item1)
            fout.write("\nTemperature(C): ")
            fout.write(item2[:-1])
            fout.write("\nFluences(1E13cm-2): ")
            [fout.write(str(item)+'\t') for item in Fluences]
            fout.write("\nSensors: ")
            [fout.write(str(item)+'\t') for item in SensorNames[item1]]
            #[fout.write(str(item)+'\t') for item in SensorNames[item1]]
            #[fout.write(str(item)+'\t') for item in SensorNames[item1]]
            #[fout.write(str(item)+'\t') for item in SensorNames[item1]]
            fout.write("\n######################################################################################################################################\n")
            fout.write("Freq.: ")
            [fout.write(Frequencies[0]+'\t') for item in SensorNames[item1]]
            [fout.write(Frequencies[0]+'\t') for item in SensorNames[item1]]
            [fout.write(Frequencies[1]+'\t') for item in SensorNames[item1]]
            [fout.write(Frequencies[1]+'\t') for item in SensorNames[item1]]
            fout.write("\nPlot: ")
            [fout.write('loglog\t') for item in SensorNames[item1]]
            [fout.write('1C2\t') for item in SensorNames[item1]]
            [fout.write('loglog\t') for item in SensorNames[item1]]
            [fout.write('1C2\t') for item in SensorNames[item1]]
            fout.write("\n######################################################################################################################################\n")
            fout.write("MEASUREMENT : DEPLETION VOLTAGE [V]\n\n")
            for i in range(len(AnnealingSteps)):
                [fout.write(str(item)+'\t') for item in VdepAll[i]]
                fout.write("\n")
            fout.write("Depletion Voltage Errors:\n")
            for i in range(len(AnnealingSteps)):
                [fout.write(str(item)+'\t') for item in VdepErAll[i]]
                fout.write("\n")
            fout.write("\n\nMEASUREMENT : DEPLETION CURRENT (0C) [uA]\n\n")
            for i in range(len(AnnealingSteps)):
                [fout.write(str(item)+'\t') for item in IdepAll[i]]
                fout.write("\n")
            fout.write("Depletion Current Errors:\n")
            for i in range(len(AnnealingSteps)):
                [fout.write(str(item)+'\t') for item in IdepErAll[i]]
                fout.write("\n")
            fout.write("\n\nEnd")


