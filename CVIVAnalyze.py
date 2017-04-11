#!/usr/bin/python

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""  File : Source code for CV-IV analysis   """""""""""""""""""""
"""""""""""""""""""""             Version : 3.17               """""""""""""""""""""
"""""""""""""""""""""          Author : Sinan Sagir            """""""""""""""""""""
"""""""""""""""""""""        Date : December 15, 2013          """""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

import os
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

def ConfigurationParameters(CVFileLocation, IVFileLocation, PlotSavingLocation, d, A, Irradiation, SensorType, ScaleTemp, FlatFittingMethod1, FlatFittingMethod2, CapManualf1, CapManualf1Err, CapManualf2, CapManualf2Err, slopef1min, slopef1max, flatf1min, flatf1max, slopef2min, slopef2max, flatf2min, flatf2max, slopef1min_1C, slopef1max_1C, flatf1min_1C, flatf1max_1C, slopef2min_1C, slopef2max_1C, flatf2min_1C, flatf2max_1C, Vmin, Vmax, CapMode, CapOffset, VdepCoeff, PlotDopingProfile, PlotDoubleDerivative, SavePopUpPlots, LCRAccuracy,OpenMeasError, CurrentAccuracy):

    if SensorType == 1:
        print "Sensor Type: P"
    if SensorType == 0:
        print "Sensor Type: N"
    if Irradiation == 1:
        print "Irradiated?: Yes\nNote: Change irradiation status to 0 if sensor thickness is wanted to be calculated from measurement!"
    if Irradiation == 0:
        print "Irradiated?: No"
    print "Current Scaling reference temperature(C): ",ScaleTemp
    if CapMode == 0:
        print "Capacitance Measurement Mode: Parallel\n"
    if CapMode == 1:
        print "Capacitance Measurement Mode: Serial\n"

    if os.path.exists(CVFileLocation):
        fCV = open(CVFileLocation, 'rU')
        linesCV = fCV.readlines()
        fCV.close()
    else:
        print "Error: CV data file path does not exist. Please check the file location.\nAborting..."
        os._exit(1)
                
    if os.path.exists(IVFileLocation):
        fIV = open(IVFileLocation, 'rU')
        linesIV = fIV.readlines()
        fIV.close()
    else:
        print "Error: IV data file path does not exist. Please check the file location.\nAborting..."
        os._exit(1)

    count = 1
    LCV1 = 0
    for line in linesCV:

        if line.startswith("Sensor Name:"):
            SenName = line[13:-1]
        if line.startswith("Annealing Status:"):
            Annealing = line[18:-1]
        if line.startswith("Environment:"):
            Environment = line[13:-1]
        if line.startswith("Note: first column reads voltage"):
            LCV1 = count
        if line.startswith("Temperature"):
            LCV2 = count
            break
        count +=1
    if LCV1 == 0:
        print "Error: Unexpected CV data file. Please check CV data file.\nAborting..."
        os._exit(1) 
    LCV1 = LCV1 + 1
    LCV2 = LCV2 - 2

    count = 1
    LIV1 = 0
    for line in linesIV:

        if line.startswith("Sensor Name:"):
            if SenName != line[13:-1]:
                print "Error: Sensor names do not match. Please check data files.\nAborting..."
                os._exit(1)            
        if line.startswith("Voltage	Current"):
            LIV1 = count
        if line.startswith("Individual IV Measurements"):
            LIV2 = count
            break
        count +=1
    if LIV1 == 0:
        print "Error: Unexpected IV data file. Please check IV data file.\nAborting..."
        os._exit(1)
    LIV1 = LIV1
    LIV2 = LIV2 - 2

    count = 0
    BVolCV = []
    Capf1paral = []
    Capf2paral = []
    Gf1 = []
    Gf2 = []
    TempCV = []
    DewTempCV = []
    for line in linesCV:
        count +=1
        if count == LCV1:
             dataCV = line.strip().split()
             freq1 = int(dataCV[1])
             freq2 = int(dataCV[2])
        if count > LCV1 and count < LCV2:
             dataCV = line.strip().split()
             if SensorType == 1:
                 BVolCV.append(-float(dataCV[0]))
             if SensorType == 0:
                 BVolCV.append(float(dataCV[0]))
             Capf1paral.append(float(dataCV[1]) + CapOffset)
             Capf2paral.append(float(dataCV[2]) + CapOffset)
             if len(dataCV) > 4:
                 Gf1.append(float(dataCV[3]))
                 Gf2.append(float(dataCV[4]))
             TempCV.append(dataCV[-1])
             if len(dataCV) > 4:
                 DewTempCV.append(dataCV[-2])

    if float(dataCV[0]) > 0 and SensorType == 1 or float(dataCV[0]) < 0 and SensorType == 0:
        print "Error: Sensor type is not correct.\nAborting..."
        os._exit(1)
        
    count = 0
    BVolIV = []
    LCur = []
    GCur = []
    TempIV = []
    DewTempIV = []
    for line in linesIV:
        count +=1
        if count > LIV1 and count < LIV2:
             dataIV = line.strip().split()
             if SensorType == 1:
                 BVolIV.append(-float(dataIV[0]))
                 LCur.append(-float(dataIV[1]))
                 GCur.append(-float(dataIV[2]))
             if SensorType == 0:
                 BVolIV.append(float(dataIV[0]))
                 LCur.append(float(dataIV[1]))
                 GCur.append(float(dataIV[2]))
             TempIV.append(dataIV[-1])
             if len(dataIV) > 4:
                 DewTempIV.append(dataIV[-2])

    """"""""" Fixing voltage column for CV measurement for P-type sensors """""""""
    """"""""" This is because of a bug in LabView code, but the bug's been fixed """""""""

    if abs(BVolIV[0]) > abs(BVolIV[-1]):
            LCur = [i for i in reversed(LCur)]
            BVolIV = [i for i in reversed(BVolIV)]
    if BVolIV[0] == 0 and BVolCV[-1] == 0:
        for i in range(len(BVolCV)):
            BVolCV[len(BVolCV)-i-1] = BVolCV[len(BVolCV)-i-2]
        BVolCV[0] = -0.0

    """"""""" PARALLEL TO SERIAL MODE CHANGE """""""""

    Capf1 = []
    Capf2 = []
    for i in range(len(Capf1paral)):
         if CapMode == 0:
             Capf1.append(Capf1paral[i])
             Capf2.append(Capf2paral[i])
         if CapMode == 1:
             if len(dataCV) > 4:
                 Capf1.append(Capf1paral[i] + Gf1[i]**2/(Capf1paral[i]*(2.*pi*freq1)**2))
                 Capf2.append(Capf2paral[i] + Gf2[i]**2/(Capf2paral[i]*(2.*pi*freq2)**2))
             else:
                 G = LCur[i]/BVolIV[i]
                 Capf1.append(Capf1paral[i] + G**2/(Capf1paral[i]*(2.*pi*freq1)**2))
                 Capf2.append(Capf2paral[i] + G**2/(Capf2paral[i]*(2.*pi*freq2)**2))

    OneCapf1 = [(1./x**2)*1e-21 for x in Capf1]
    OneCapf2 = [(1./x**2)*1e-21 for x in Capf2]
    TempCV = map(float, TempCV)
    DewTempCV = map(float, DewTempCV)
    TempIV = map(float, TempIV)
    DewTempIV = map(float, DewTempIV)

    """"""""" MODIFIYING DATA FILE FOR P-TYPE SENSORS """""""""

    if SensorType == 1:
        ModifiedFileLocationCV = CVFileLocation[:-4]
        ModifiedFileLocationIV = IVFileLocation[:-4]
        fCVmod = open("%(ModifiedFileLocationCV)s_modified.txt" %locals(), "w")
        fIVmod = open("%(ModifiedFileLocationIV)s_modified.txt" %locals(), "w")
        count = 0
        for line in linesCV:
            count += 1
            if count <= LCV1:
                fCVmod.write( str(line)  )
            if count >= LCV2:
                fCVmod.write( str(line)  )
            if count > LCV1 and count < LCV2:
                if len(dataCV) > 4:
                    colCV0 = BVolCV[count-LCV1-1]
                    colCV1 = Capf1[count-LCV1-1]
                    colCV2 = Capf2[count-LCV1-1]
                    colCV3 = Gf1[count-LCV1-1]
                    colCV4 = Gf2[count-LCV1-1]
                    colCV5 = DewTempCV[count-LCV1-1]
                    colCV6 = TempCV[count-LCV1-1]
                    fCVmod.write( str(colCV0) + "\t" + str(colCV1) + "\t" + str(colCV2) + "\t" + str(colCV3) + "\t" + str(colCV4) + "\t" + str(colCV5) + "\t" + str(colCV6) + "\n"  )     
                else:
                    colCV0 = BVolCV[count-LCV1-1]
                    colCV1 = Capf1[count-LCV1-1]
                    colCV2 = Capf2[count-LCV1-1]
                    colCV3 = DewTempCV[count-LCV1-1]
                    colCV4 = TempCV[count-LCV1-1]
                    fCVmod.write( str(colCV0) + "\t" + str(colCV1) + "\t" + str(colCV2) + "\t" + str(colCV3) + "\t" + str(colCV4) + "\n"  )
        fCVmod.close()
        count = 0
        for line in linesIV:
            count += 1
            if count <= LIV1:
                fIVmod.write( str(line)  )
            if count >= LIV2:
                fIVmod.write( str(line)  )
            if count > LIV1 and count < LIV2:
                colIV0 = BVolIV[count-LIV1-1]
                colIV1 = LCur[count-LIV1-1]
                colIV2 = GCur[count-LIV1-1]
                colIV3 = DewTempIV[count-LIV1-1]
                colIV4 = TempIV[count-LIV1-1]
                fIVmod.write( str(colIV0) + "\t" + str(colIV1) + "\t" + str(colIV2) + "\t" + str(colIV3) + "\t" + str(colIV4) + "\n"  )
        fIVmod.close()

    """"""""" Finding depletion voltage from double derivative (this method is not working """"""
    """"""    properly b/c of low resolution of bias voltage in data) """""""""""""""""""""""""""

    DCapf1 = []
    DCapf2 = []
    for i in range(1,len(BVolCV)-1):
        Df1 = (Capf1[i+1]+Capf1[i-1]-2.*Capf1[i])/(BVolCV[i+1]-BVolCV[i])**2
        Df2 = (Capf2[i+1]+Capf2[i-1]-2.*Capf2[i])/(BVolCV[i+1]-BVolCV[i])**2
        DCapf1.append(Df1)
        DCapf2.append(Df2)

    """"""""""""""""""""" SCALING CURRENT """""""""""""""""""""

    Eg0 = 1.166 # eV # 1.206 eV in http://indico.cern.ch/event/129737/session/3/contribution/24/material/slides/1.pdf
    alpha = 4.74 * 1e-4 # eV/K
    beta = 636.0 # K
    def ScaleCurrent(ScaleTemp, Temp):
        Eg = Eg0 - alpha * Temp**2 / (Temp + beta)
        return pow((ScaleTemp+273.15)/(Temp+273.15),2)*exp((-Eg/(2.*8.6173324*pow(10,-5)))*(1./(ScaleTemp+273.15) - 1./(Temp+273.15)))

    LCurScaled = []
    for i in range(len(TempIV)):
        LCurScaled.append(LCur[i] * ScaleCurrent(ScaleTemp, TempIV[i]))
    GCurScaled = []
    for i in range(len(TempIV)):
        GCurScaled.append(GCur[i] * ScaleCurrent(ScaleTemp, TempIV[i]))

    VminIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vmin)))
    VmaxIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vmax)))

    sf1minIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef1min)))
    sf1maxIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef1max))) + 1
    ff1minIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf1min)))
    ff1maxIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf1max))) + 1

    sf2minIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef2min)))
    sf2maxIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef2max))) + 1
    ff2minIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf2min)))
    ff2maxIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf2max))) + 1

    sf1min1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef1min_1C)))
    sf1max1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef1max_1C))) + 1
    ff1min1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf1min_1C)))
    ff1max1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf1max_1C))) + 1

    sf2min1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef2min_1C)))
    sf2max1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-slopef2max_1C))) + 1
    ff2min1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf2min_1C)))
    ff2max1CIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-flatf2max_1C))) + 1

    """"""""" DEFINING FITTING FUNCTIONS & PERFORMING FITS """""""""

    def LogSlopefunc(voltage, parameter0, parameter1):
        return parameter0 * (voltage**parameter1)

    def Linearfunc(voltage, parameter0, parameter1):
        return parameter0 + parameter1 * voltage

    def Horizontalfunc(voltage, parameter0):
        return parameter0

    LCRerr = LCRAccuracy/100 #%
    OpenErr = OpenMeasError # F
    CurErr = CurrentAccuracy/100 # %

    def FittingRangeCheck(datapoints, mindatapoint):
        if len(datapoints) <= mindatapoint:
            print "Error: Cannot perform fitting. Please change your fitting range.\nAborting..."
            os._exit(1)    

    logvolslopef1 = array(BVolCV[sf1minIdx:sf1maxIdx])
    logcapslopef1 = array(Capf1[sf1minIdx:sf1maxIdx])
    logcapslopef1err = array([sqrt((LCRerr*x)**2 + OpenErr**2) for x in logcapslopef1])
    FittingRangeCheck(logvolslopef1, 2)

    logvolflatf1 = array(BVolCV[ff1minIdx:ff1maxIdx])
    logcapflatf1 = array(Capf1[ff1minIdx:ff1maxIdx])
    logcapflatf1err = array([sqrt((LCRerr*x)**2 + OpenErr**2) for x in logcapflatf1])
    FittingRangeCheck(logvolflatf1, 1)

    logvolslopef2 = array(BVolCV[sf2minIdx:sf2maxIdx])
    logcapslopef2 = array(Capf2[sf2minIdx:sf2maxIdx])
    logcapslopef2err = array([sqrt((LCRerr*x)**2 + OpenErr**2) for x in logcapslopef2])
    FittingRangeCheck(logvolslopef2, 2)

    logvolflatf2 = array(BVolCV[ff2minIdx:ff2maxIdx])
    logcapflatf2 = array(Capf2[ff2minIdx:ff2maxIdx])
    logcapflatf2err = array([sqrt((LCRerr*x)**2 + OpenErr**2) for x in logcapflatf2])
    FittingRangeCheck(logvolflatf2, 1)

    if slopef1min_1C < 0:
        ICvolslopef1 = array(BVolCV[sf1minIdx:sf1maxIdx])
        ICcapslopef1 = array(OneCapf1[sf1minIdx:sf1maxIdx])
        ICcapslopef1err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapslopef1])
        FittingRangeCheck(ICvolslopef1, 2)

        ICvolflatf1 = array(BVolCV[ff1minIdx:ff1maxIdx])
        ICcapflatf1 = array(OneCapf1[ff1minIdx:ff1maxIdx])
        ICcapflatf1err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapflatf1])
        FittingRangeCheck(ICvolflatf1, 2)

        ICvolslopef2 = array(BVolCV[sf2minIdx:sf2maxIdx])
        ICcapslopef2 = array(OneCapf2[sf2minIdx:sf2maxIdx])
        ICcapslopef2err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapslopef2])
        FittingRangeCheck(ICvolslopef2, 2)

        ICvolflatf2 = array(BVolCV[ff2minIdx:ff2maxIdx])
        ICcapflatf2 = array(OneCapf2[ff2minIdx:ff2maxIdx])
        ICcapflatf2err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapflatf2])
        FittingRangeCheck(ICvolflatf2, 2)

    else:
        ICvolslopef1 = array(BVolCV[sf1min1CIdx:sf1max1CIdx])
        ICcapslopef1 = array(OneCapf1[sf1min1CIdx:sf1max1CIdx])
        ICcapslopef1err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapslopef1])
        FittingRangeCheck(ICvolslopef1, 2)

        ICvolflatf1 = array(BVolCV[ff1min1CIdx:ff1max1CIdx])
        ICcapflatf1 = array(OneCapf1[ff1min1CIdx:ff1max1CIdx])
        ICcapflatf1err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapflatf1])
        FittingRangeCheck(ICvolflatf1, 2)

        ICvolslopef2 = array(BVolCV[sf2min1CIdx:sf2max1CIdx])
        ICcapslopef2 = array(OneCapf2[sf2min1CIdx:sf2max1CIdx])
        ICcapslopef2err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapslopef2])
        FittingRangeCheck(ICvolslopef2, 2)

        ICvolflatf2 = array(BVolCV[ff2min1CIdx:ff2max1CIdx])
        ICcapflatf2 = array(OneCapf2[ff2min1CIdx:ff2max1CIdx])
        ICcapflatf2err = array([sqrt((2.*LCRerr*x)**2 + 4.*(OpenErr**2)*(x**3)) for x in ICcapflatf2])
        FittingRangeCheck(ICvolflatf2, 2)

    def general_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
        popt, pcov = curve_fit(f, xdata, ydata, p0, sigma, maxfev=10000)

        if sigma is None:
            chi2 = sum(((f(xdata,*popt)-ydata))**2)
        else:
            chi2 = sum(((f(xdata,*popt)-ydata)/sigma)**2)
        dof = len(ydata) - len(popt)
        rchi2 = chi2/dof
        punc = zeros(len(popt))
        for i in arange(0,len(popt)):
            punc[i] = sqrt(pcov[i,i])
        return popt, punc, rchi2, dof

    logslopepar0f1 = [0.01, 0.5]
    poptlogslopef1, punclogslopef1, rclogslopef1, dlogslopef1 = general_fit(LogSlopefunc,
                                                                            logvolslopef1,
                                                                            logcapslopef1,
                                                                            logslopepar0f1,
                                                                            logcapslopef1err)

    if CapManualf1 < 0:
        if FlatFittingMethod1 == 0:
            logflatpar0f1 = [0.01, 0.5]
            poptlogflatf1, punclogflatf1, rclogflatf1, dlogflatf1 = general_fit(LogSlopefunc,
                                                                                logvolflatf1,
                                                                                logcapflatf1,
                                                                                logflatpar0f1,
                                                                                logcapflatf1err)
        else:
            logflatpar0f1 = [0.01]
            poptlogflatf1, punclogflatf1, rclogflatf1, dlogflatf1 = general_fit(Horizontalfunc,
                                                                                logvolflatf1,
                                                                                logcapflatf1,
                                                                                logflatpar0f1,
                                                                                logcapflatf1err)
    if CapManualf1 > 0:
        poptlogflatf1 = [CapManualf1]
        punclogflatf1 = [CapManualf1Err]

    logslopepar0f2 = [0.01, 0.5]
    poptlogslopef2, punclogslopef2, rclogslopef2, dlogslopef2 = general_fit(LogSlopefunc,
                                                                            logvolslopef2,
                                                                            logcapslopef2,
                                                                            logslopepar0f2,
                                                                            logcapslopef2err)

    if CapManualf2 < 0:
        if FlatFittingMethod1 == 0:
            logflatpar0f2 = [0.01, 0.5]
            poptlogflatf2, punclogflatf2, rclogflatf2, dlogflatf2 = general_fit(LogSlopefunc,
                                                                                logvolflatf2,
                                                                                logcapflatf2,
                                                                                logflatpar0f2,
                                                                                logcapflatf2err)
        else:
            logflatpar0f2 = [0.01]
            poptlogflatf2, punclogflatf2, rclogflatf2, dlogflatf2 = general_fit(Horizontalfunc,
                                                                                logvolflatf2,
                                                                                logcapflatf2,
                                                                                logflatpar0f2,
                                                                                logcapflatf2err)
    if CapManualf2 > 0:
        poptlogflatf2 = [CapManualf2]
        punclogflatf2 = [CapManualf2Err]

    ICslopepar0f1 = [0.01, 0.5]
    poptICslopef1, puncICslopef1, rcICslopef1, dICslopef1 = general_fit(Linearfunc,
                                                                        ICvolslopef1,
                                                                        ICcapslopef1,
                                                                        ICslopepar0f1,
                                                                        ICcapslopef1err)

    if CapManualf1 < 0:
        if FlatFittingMethod2 == 0:
            ICflatpar0f1 = [0.01, 0.5]
            poptICflatf1, puncICflatf1, rcICflatf1, dICflatf1 = general_fit(Linearfunc,
                                                                            ICvolflatf1,
                                                                            ICcapflatf1,
                                                                            ICflatpar0f1,
                                                                            ICcapflatf1err)
        else:
            ICflatpar0f1 = [0.01]
            poptICflatf1, puncICflatf1, rcICflatf1, dICflatf1 = general_fit(Horizontalfunc,
                                                                            ICvolflatf1,
                                                                            ICcapflatf1,
                                                                            ICflatpar0f1,
                                                                            ICcapflatf1err)
    if CapManualf1 > 0:
        poptICflatf1 = [1E-21/CapManualf1**2, 0]
        puncICflatf1 = [CapManualf1Err*2E-21/CapManualf1**3, 0]

    ICslopepar0f2 = [0.01, 0.5]
    poptICslopef2, puncICslopef2, rcICslopef2, dICslopef2 = general_fit(Linearfunc,
                                                                        ICvolslopef2,
                                                                        ICcapslopef2,
                                                                        ICslopepar0f2,
                                                                        ICcapslopef2err)

    if CapManualf2 < 0:
        if FlatFittingMethod2 == 0:
            ICflatpar0f2 = [0.01, 0.5]
            poptICflatf2, puncICflatf2, rcICflatf2, dICflatf2 = general_fit(Linearfunc,
                                                                            ICvolflatf2,
                                                                            ICcapflatf2,
                                                                            ICflatpar0f2,
                                                                            ICcapflatf2err)
        else:
            ICflatpar0f2 = [0.01]
            poptICflatf2, puncICflatf2, rcICflatf2, dICflatf2 = general_fit(Horizontalfunc,
                                                                            ICvolflatf2,
                                                                            ICcapflatf2,
                                                                            ICflatpar0f2,
                                                                            ICcapflatf2err)
    if CapManualf2 > 0:
        poptICflatf2 = [1E-21/CapManualf2**2, 0]
        puncICflatf2 = [CapManualf2Err*2E-21/CapManualf2**3, 0]

    """"""""" CALCULATING SENSOR CHARACTERISTICS """""""""

    if FlatFittingMethod1 == 0:
        Vdepf1Log = pow(poptlogflatf1[0] / poptlogslopef1[0],1. / (poptlogslopef1[1] - poptlogflatf1[1]))
        Vdepf1LogEr = (Vdepf1Log/abs(poptlogslopef1[1]-poptlogflatf1[1]))*sqrt(pow(punclogslopef1[0]/poptlogslopef1[0],2) +
                                                                               pow(punclogflatf1[0]/poptlogflatf1[0],2) +
                                                                               pow(poptlogslopef1[1] - poptlogflatf1[1],-2) *
                                                                               pow(log(poptlogflatf1[0]/poptlogslopef1[0]),2) *
                                                                               (pow(punclogflatf1[1],2)+pow(punclogslopef1[1],2)))

        Vdepf2Log = pow(poptlogflatf2[0] / poptlogslopef2[0],1. / (poptlogslopef2[1] - poptlogflatf2[1]))
        Vdepf2LogEr = (Vdepf2Log/abs(poptlogslopef2[1]-poptlogflatf2[1]))*sqrt(pow(punclogslopef2[0]/poptlogslopef2[0],2) +
                                                                               pow(punclogflatf2[0]/poptlogflatf2[0],2) +
                                                                               pow(poptlogslopef2[1] - poptlogflatf2[1],-2) *
                                                                               pow(log(poptlogflatf2[0]/poptlogslopef2[0]),2) *
                                                                               (pow(punclogflatf2[1],2)+pow(punclogslopef2[1],2)))

    if FlatFittingMethod1 == 1:
        Vdepf1Log = pow(poptlogflatf1[0] / poptlogslopef1[0],1. / poptlogslopef1[1])
        Vdepf1LogEr = (Vdepf1Log/abs(poptlogslopef1[1]))*sqrt(pow(punclogslopef1[0]/poptlogslopef1[0],2) +
                                                              pow(punclogflatf1[0]/poptlogflatf1[0],2) +
                                                              pow(punclogslopef1[1]/poptlogslopef1[1],2) *
                                                              pow(log(poptlogflatf1[0]/poptlogslopef1[0]),2))

        Vdepf2Log = pow(poptlogflatf2[0] / poptlogslopef2[0],1. / poptlogslopef2[1])
        Vdepf2LogEr = (Vdepf2Log/abs(poptlogslopef2[1]))*sqrt(pow(punclogslopef2[0]/poptlogslopef2[0],2) +
                                                              pow(punclogflatf2[0]/poptlogflatf2[0],2) +
                                                              pow(punclogslopef2[1]/poptlogslopef2[1],2) *
                                                              pow(log(poptlogflatf2[0]/poptlogslopef2[0]),2))

    if FlatFittingMethod2 == 0:
        Vdepf1IC = (poptICslopef1[0] - poptICflatf1[0])/(poptICflatf1[1] - poptICslopef1[1])
        Vdepf1ICEr = (1./(poptICflatf1[1] - poptICslopef1[1])**2)*sqrt((puncICslopef1[0]**2 + puncICflatf1[0]**2) *
                                                                       (poptICflatf1[1] - poptICslopef1[1])**2 +
                                                                       (puncICflatf1[1]**2 + puncICslopef1[1]**2) *
                                                                       (poptICslopef1[0] - poptICflatf1[0])**2)

        Vdepf2IC = (poptICslopef2[0] - poptICflatf2[0])/(poptICflatf2[1] - poptICslopef2[1])
        Vdepf2ICEr = (1./(poptICflatf2[1] - poptICslopef2[1])**2)*sqrt((puncICslopef2[0]**2 + puncICflatf2[0]**2) *
                                                                       (poptICflatf2[1] - poptICslopef2[1])**2 +
                                                                       (puncICflatf2[1]**2 + puncICslopef2[1]**2) *
                                                                       (poptICslopef2[0] - poptICflatf2[0])**2)

    if FlatFittingMethod2 == 1:
        Vdepf1IC = (poptICflatf1[0] - poptICslopef1[0])/poptICslopef1[1]
        Vdepf1ICEr = sqrt(pow(puncICslopef1[0]/poptICslopef1[1],2) +
                          pow(puncICflatf1[0]/poptICslopef1[1],2) +
                          pow(puncICslopef1[1]*pow(puncICslopef1[0] - puncICflatf1[0],2)/poptICslopef1[1],2))

        Vdepf2IC = (poptICflatf2[0] - poptICslopef2[0])/poptICslopef2[1]
        Vdepf2ICEr = sqrt(pow(puncICslopef2[0]/poptICslopef2[1],2) +
                          pow(puncICflatf2[0]/poptICslopef2[1],2) +
                          pow(puncICslopef2[1]*pow(puncICslopef2[0] - puncICflatf2[0],2)/poptICslopef2[1],2))

    epsi0 = 8.85418 * pow(10,-12)
    epsi = 11.9 # for Silicon
    qe = 1.60218 * 1e-19

    if Irradiation == 1:
        df1Log = d
        df2Log = d
        df1IC = d
        df2IC = d
        EndCapf1Log = epsi0 * epsi * A / df1Log
        EndCapf2Log = epsi0 * epsi * A / df2Log
        EndCapf1IC = epsi0 * epsi * A / df1IC
        EndCapf2IC = epsi0 * epsi * A / df2IC
        EndCapf1LogErr = 0.0
        EndCapf2LogErr = 0.0
        EndCapf1ICErr = 0.0
        EndCapf2ICErr = 0.0
        df1LogErr = 0.0
        df2LogErr = 0.0
        df1ICErr = 0.0
        df2ICErr = 0.0

    if Irradiation == 0:
        EndCapf1Log = poptlogflatf1[0]
        EndCapf2Log = poptlogflatf2[0]
        EndCapf1IC = 1./sqrt(poptICflatf1[0]*1.e21)
        EndCapf2IC = 1./sqrt(poptICflatf2[0]*1.e21)
        df1Log = epsi0 * epsi * A / EndCapf1Log
        df2Log = epsi0 * epsi * A / EndCapf2Log
        df1IC = epsi0 * epsi * A / EndCapf1IC
        df2IC = epsi0 * epsi * A / EndCapf2IC
        EndCapf1LogErr = punclogflatf1[0]#sqrt((LCRerr*EndCapf1Log)**2 + OpenErr**2)
        EndCapf2LogErr = punclogflatf1[0]#sqrt((LCRerr*EndCapf2Log)**2 + OpenErr**2)
        EndCapf1ICErr = (puncICflatf1[0]/poptICflatf1[0])*(0.5/sqrt(poptICflatf1[0]*1.e21))
        EndCapf2ICErr = (puncICflatf2[0]/poptICflatf2[0])*(0.5/sqrt(poptICflatf2[0]*1.e21))
        df1LogErr = epsi0 * epsi * A * EndCapf1LogErr / EndCapf1Log**2
        df2LogErr = epsi0 * epsi * A * EndCapf2LogErr / EndCapf2Log**2
        df1ICErr = epsi0 * epsi * A * EndCapf1ICErr / EndCapf1IC**2
        df2ICErr = epsi0 * epsi * A * EndCapf2ICErr / EndCapf2IC**2

    Neff1Log = 2.0 * epsi * epsi0 * Vdepf1Log * 1e-17 / (qe * df1Log**2)
    Neff1LogEr = (2.0 * epsi * epsi0 * 1e-17 / (qe * df1Log**2)) * sqrt(Vdepf1LogEr**2 + 4.*((Vdepf1Log/df1Log)**2)*(df1LogErr**2))
    Neff2Log = 2.0 * epsi * epsi0 * Vdepf2Log * 1e-17 / (qe * df1Log**2)
    Neff2LogEr = (2.0 * epsi * epsi0 * 1e-17 / (qe * df2Log**2)) * sqrt(Vdepf2LogEr**2 + 4.*((Vdepf2Log/df2Log)**2)*(df2LogErr**2))
    Neff1IC = 2.0 * epsi * epsi0 * Vdepf1IC * 1e-17 / (qe * df1IC**2)
    Neff1ICEr = (2.0 * epsi * epsi0 * 1e-17 / (qe * df1IC**2)) * sqrt(Vdepf1ICEr**2 + 4.*((Vdepf1IC/df1IC)**2)*(df1ICErr**2))
    Neff2IC = 2.0 * epsi * epsi0 * Vdepf2IC * 1e-17 / (qe * df2IC**2)
    Neff2ICEr = (2.0 * epsi * epsi0 * 1e-17 / (qe * df2IC**2)) * sqrt(Vdepf2ICEr**2 + 4.*((Vdepf2IC/df2IC)**2)*(df2ICErr**2))

    Vdepf1LogIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vdepf1Log))) 
    Vdepf2LogIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vdepf2Log)))
    Vdepf1ICIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vdepf1IC)))
    Vdepf2ICIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-Vdepf2IC)))

    Vpostdepf1LogIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-VdepCoeff*Vdepf1Log)))
    Vpostdepf2LogIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-VdepCoeff*Vdepf2Log)))
    Vpostdepf1ICIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-VdepCoeff*Vdepf1IC)))
    Vpostdepf2ICIdx = BVolCV.index(min(BVolCV, key=lambda x:abs(x-VdepCoeff*Vdepf2IC)))

    if Vdepf1LogIdx == BVolCV.index(BVolCV[-1]):
        Vdepf1LogIdx1 = Vdepf1LogIdx-1
    elif Vdepf1LogIdx == BVolCV.index(BVolCV[0]):
        Vdepf1LogIdx1 = Vdepf1LogIdx+1
    else:
        Vdepf1LogIdx1 = BVolCV.index(min([BVolCV[Vdepf1LogIdx-1],BVolCV[Vdepf1LogIdx+1]], key=lambda x:abs(x-Vdepf1Log)))
        
    if Vdepf2LogIdx == BVolCV.index(BVolCV[-1]):
        Vdepf2LogIdx1 = Vdepf2LogIdx-1
    elif Vdepf2LogIdx == BVolCV.index(BVolCV[0]):
        Vdepf2LogIdx1 = Vdepf2LogIdx+1
    else:
        Vdepf2LogIdx1 = BVolCV.index(min([BVolCV[Vdepf2LogIdx-1],BVolCV[Vdepf2LogIdx+1]], key=lambda x:abs(x-Vdepf2Log)))

    if Vdepf1ICIdx == BVolCV.index(BVolCV[-1]):
        Vdepf1ICIdx1 = Vdepf1ICIdx-1
    elif Vdepf1ICIdx == BVolCV.index(BVolCV[0]):
        Vdepf1ICIdx1 = Vdepf1ICIdx+1
    else:
        Vdepf1ICIdx1 = BVolCV.index(min([BVolCV[Vdepf1ICIdx-1],BVolCV[Vdepf1ICIdx+1]], key=lambda x:abs(x-Vdepf1IC)))

    if Vdepf2ICIdx == BVolCV.index(BVolCV[-1]):
        Vdepf2ICIdx1 = Vdepf2ICIdx-1
    elif Vdepf2ICIdx == BVolCV.index(BVolCV[0]):
        Vdepf2ICIdx1 = Vdepf2ICIdx+1
    else:
        Vdepf2ICIdx1 = BVolCV.index(min([BVolCV[Vdepf2ICIdx-1],BVolCV[Vdepf2ICIdx+1]], key=lambda x:abs(x-Vdepf2IC)))


    if Vpostdepf1LogIdx == BVolCV.index(BVolCV[-1]):
        Vpostdepf1LogIdx1 = Vpostdepf1LogIdx-1
    elif Vpostdepf1LogIdx == BVolCV.index(BVolCV[0]):
        Vpostdepf1LogIdx1 = Vpostdepf1LogIdx+1
    else:
        Vpostdepf1LogIdx1 = BVolCV.index(min([BVolCV[Vpostdepf1LogIdx-1],BVolCV[Vpostdepf1LogIdx+1]], key=lambda x:abs(x-VdepCoeff*Vdepf1Log)))

    if Vpostdepf2LogIdx == BVolCV.index(BVolCV[-1]):
        Vpostdepf2LogIdx1 = Vpostdepf2LogIdx-1
    elif Vpostdepf2LogIdx == BVolCV.index(BVolCV[0]):
        Vpostdepf2LogIdx1 = Vpostdepf2LogIdx+1
    else:
        Vpostdepf2LogIdx1 = BVolCV.index(min([BVolCV[Vpostdepf2LogIdx-1],BVolCV[Vpostdepf2LogIdx+1]], key=lambda x:abs(x-VdepCoeff*Vdepf2Log)))

    if Vpostdepf1ICIdx == BVolCV.index(BVolCV[-1]):
        Vpostdepf1ICIdx1 = Vpostdepf1ICIdx-1
    elif Vpostdepf1ICIdx == BVolCV.index(BVolCV[0]):
        Vpostdepf1ICIdx1 = Vpostdepf1LogIdx+1
    else:
        Vpostdepf1ICIdx1 = BVolCV.index(min([BVolCV[Vpostdepf1ICIdx-1],BVolCV[Vpostdepf1ICIdx+1]], key=lambda x:abs(x-VdepCoeff*Vdepf1IC)))
        
    if Vpostdepf2ICIdx == BVolCV.index(BVolCV[-1]):
        Vpostdepf2ICIdx1 = Vpostdepf2ICIdx-1
    elif Vpostdepf2ICIdx == BVolCV.index(BVolCV[0]):
        Vpostdepf2ICIdx1 = Vpostdepf2LogIdx+1
    else:
        Vpostdepf2ICIdx1 = BVolCV.index(min([BVolCV[Vpostdepf2ICIdx-1],BVolCV[Vpostdepf2ICIdx+1]], key=lambda x:abs(x-VdepCoeff*Vdepf2IC)))

    Vdepf1LogIdx_b = LCur.index(min([LCur[Vdepf1LogIdx],LCur[Vdepf1LogIdx1]]))# b for below
    Vdepf2LogIdx_b = LCur.index(min([LCur[Vdepf2LogIdx],LCur[Vdepf2LogIdx1]]))
    Vdepf1ICIdx_b = LCur.index(min([LCur[Vdepf1ICIdx],LCur[Vdepf1ICIdx1]]))
    Vdepf2ICIdx_b = LCur.index(min([LCur[Vdepf2ICIdx],LCur[Vdepf2ICIdx1]]))

    Vpostdepf1LogIdx_b = LCur.index(min([LCur[Vpostdepf1LogIdx],LCur[Vpostdepf1LogIdx1]]))
    Vpostdepf2LogIdx_b = LCur.index(min([LCur[Vpostdepf2LogIdx],LCur[Vpostdepf2LogIdx1]]))
    Vpostdepf1ICIdx_b = LCur.index(min([LCur[Vpostdepf1ICIdx],LCur[Vpostdepf1ICIdx1]]))
    Vpostdepf2ICIdx_b = LCur.index(min([LCur[Vpostdepf2ICIdx],LCur[Vpostdepf2ICIdx1]]))

    Idepf1Log = LCur[Vdepf1LogIdx_b]+abs((LCur[Vdepf1LogIdx]-LCur[Vdepf1LogIdx1])*(Vdepf1Log-BVolIV[Vdepf1LogIdx_b])/(BVolIV[Vdepf1LogIdx]-BVolIV[Vdepf1LogIdx1]))#LCur[Vdepf1LogIdx]
    Idepf1LogScaled = LCurScaled[Vdepf1LogIdx_b]+abs((LCurScaled[Vdepf1LogIdx]-LCurScaled[Vdepf1LogIdx1])*(Vdepf1Log-BVolIV[Vdepf1LogIdx_b])/(BVolIV[Vdepf1LogIdx]-BVolIV[Vdepf1LogIdx1]))#LCurScaled[Vdepf1LogIdx]
    Idepf2Log = LCur[Vdepf2LogIdx_b]+abs((LCur[Vdepf2LogIdx]-LCur[Vdepf2LogIdx1])*(Vdepf2Log-BVolIV[Vdepf2LogIdx_b])/(BVolIV[Vdepf2LogIdx]-BVolIV[Vdepf2LogIdx1]))#LCur[Vdepf2LogIdx]
    Idepf2LogScaled = LCurScaled[Vdepf2LogIdx_b]+abs((LCurScaled[Vdepf2LogIdx]-LCurScaled[Vdepf2LogIdx1])*(Vdepf2Log-BVolIV[Vdepf2LogIdx_b])/(BVolIV[Vdepf2LogIdx]-BVolIV[Vdepf2LogIdx1]))#LCurScaled[Vdepf2LogIdx]
    Idepf1IC = LCur[Vdepf1ICIdx_b]+abs((LCur[Vdepf1ICIdx]-LCur[Vdepf1ICIdx1])*(Vdepf1IC-BVolIV[Vdepf1ICIdx_b])/(BVolIV[Vdepf1ICIdx]-BVolIV[Vdepf1ICIdx1]))#LCur[Vdepf1ICIdx]
    Idepf1ICScaled = LCurScaled[Vdepf1ICIdx_b]+abs((LCurScaled[Vdepf1ICIdx]-LCurScaled[Vdepf1ICIdx1])*(Vdepf1IC-BVolIV[Vdepf1ICIdx_b])/(BVolIV[Vdepf1ICIdx]-BVolIV[Vdepf1ICIdx1]))#LCurScaled[Vdepf1ICIdx]
    Idepf2IC = LCur[Vdepf2ICIdx_b]+abs((LCur[Vdepf2ICIdx]-LCur[Vdepf2ICIdx1])*(Vdepf2IC-BVolIV[Vdepf2ICIdx_b])/(BVolIV[Vdepf2ICIdx]-BVolIV[Vdepf2ICIdx1]))#LCur[Vdepf2ICIdx]
    Idepf2ICScaled = LCurScaled[Vdepf2ICIdx_b]+abs((LCurScaled[Vdepf2ICIdx]-LCurScaled[Vdepf2ICIdx1])*(Vdepf2IC-BVolIV[Vdepf2ICIdx_b])/(BVolIV[Vdepf2ICIdx]-BVolIV[Vdepf2ICIdx1]))#LCurScaled[Vdepf2ICIdx]

    Ipostdepf1Log = LCur[Vpostdepf1LogIdx_b]+abs((LCur[Vpostdepf1LogIdx]-LCur[Vpostdepf1LogIdx1])*(Vdepf1Log*VdepCoeff-BVolIV[Vpostdepf1LogIdx_b])/(BVolIV[Vpostdepf1LogIdx]-BVolIV[Vpostdepf1LogIdx1]))#LCur[Vpostdepf1LogIdx]
    Ipostdepf1LogScaled = LCurScaled[Vpostdepf1LogIdx_b]+abs((LCurScaled[Vpostdepf1LogIdx]-LCurScaled[Vpostdepf1LogIdx1])*(Vdepf1Log*VdepCoeff-BVolIV[Vpostdepf1LogIdx_b])/(BVolIV[Vpostdepf1LogIdx]-BVolIV[Vpostdepf1LogIdx1]))#LCurScaled[Vpostdepf1LogIdx]
    Ipostdepf2Log = LCur[Vpostdepf2LogIdx_b]+abs((LCur[Vpostdepf2LogIdx]-LCur[Vpostdepf2LogIdx1])*(Vdepf2Log*VdepCoeff-BVolIV[Vpostdepf2LogIdx_b])/(BVolIV[Vpostdepf2LogIdx]-BVolIV[Vpostdepf2LogIdx1]))#LCur[Vpostdepf2LogIdx]
    Ipostdepf2LogScaled = LCurScaled[Vpostdepf2LogIdx_b]+abs((LCurScaled[Vpostdepf2LogIdx]-LCurScaled[Vpostdepf2LogIdx1])*(Vdepf2Log*VdepCoeff-BVolIV[Vpostdepf2LogIdx_b])/(BVolIV[Vpostdepf2LogIdx]-BVolIV[Vpostdepf2LogIdx1]))#LCurScaled[Vpostdepf2LogIdx]
    Ipostdepf1IC = LCur[Vpostdepf1ICIdx_b]+abs((LCur[Vpostdepf1ICIdx]-LCur[Vpostdepf1ICIdx1])*(Vdepf1IC*VdepCoeff-BVolIV[Vpostdepf1ICIdx_b])/(BVolIV[Vpostdepf1ICIdx]-BVolIV[Vpostdepf1ICIdx1]))#LCur[Vpostdepf1ICIdx]
    Ipostdepf1ICScaled = LCurScaled[Vpostdepf1ICIdx_b]+abs((LCurScaled[Vpostdepf1ICIdx]-LCurScaled[Vpostdepf1ICIdx1])*(Vdepf1IC*VdepCoeff-BVolIV[Vpostdepf1ICIdx_b])/(BVolIV[Vpostdepf1ICIdx]-BVolIV[Vpostdepf1ICIdx1]))#LCurScaled[Vpostdepf1ICIdx]
    Ipostdepf2IC = LCur[Vpostdepf2ICIdx_b]+abs((LCur[Vpostdepf2ICIdx]-LCur[Vpostdepf2ICIdx1])*(Vdepf2IC*VdepCoeff-BVolIV[Vpostdepf2ICIdx_b])/(BVolIV[Vpostdepf2ICIdx]-BVolIV[Vpostdepf2ICIdx1]))#LCur[Vpostdepf2ICIdx]
    Ipostdepf2ICScaled = LCurScaled[Vpostdepf2ICIdx_b]+abs((LCurScaled[Vpostdepf2ICIdx]-LCurScaled[Vpostdepf2ICIdx1])*(Vdepf2IC*VdepCoeff-BVolIV[Vpostdepf2ICIdx_b])/(BVolIV[Vpostdepf2ICIdx]-BVolIV[Vpostdepf2ICIdx1]))#LCurScaled[Vpostdepf2ICIdx]

    Cpostdepf1Log = Capf1[Vpostdepf1LogIdx]
    Cpostdepf2Log = Capf2[Vpostdepf2LogIdx]
    Cpostdepf1IC = Capf1[Vpostdepf1ICIdx]
    Cpostdepf2IC = Capf2[Vpostdepf2ICIdx]
    Cpostdepf1LogErr = sqrt((LCRerr*Cpostdepf1Log)**2 + OpenErr**2)
    Cpostdepf2LogErr = sqrt((LCRerr*Cpostdepf2Log)**2 + OpenErr**2)
    Cpostdepf1ICErr = sqrt((LCRerr*Cpostdepf1IC)**2 + OpenErr**2)
    Cpostdepf2ICErr = sqrt((LCRerr*Cpostdepf2IC)**2 + OpenErr**2)

    """Finding the break-down voltage"""
    Vbd = 0
    for i in range(Vdepf2LogIdx+10,len(LCurScaled)):   
        deltaLI = ((LCurScaled[i-1]-LCurScaled[i])/LCurScaled[i]) / (BVolIV[i-1]-BVolIV[i])
        if deltaLI > 0.01: # this is given as 0.1 in HPK measurement specifications V_2.7
            Vbd = int(BVolIV[i-1])
            break
    if Vbd == 0:
        Vbd = int(BVolIV[-1])
        VbdStr1 = ">"
        VbdStr2 = ""
    else:
        VbdStr1 = "\\approx"
        VbdStr2 = "$--yellow$ $line$"
        
    VbdGR = 0
    for i in range(Vdepf2LogIdx+10,len(GCurScaled)):   
        deltaGI = ((GCurScaled[i-1]-GCurScaled[i])/(GCurScaled[i]+1e-21)) / (BVolIV[i-1]-BVolIV[i])
        if deltaGI > 0.009: # this is given as 0.1 in HPK measurement specifications V_2.7
            VbdGR = int(BVolIV[i-3])
            break
    if VbdGR == 0:
        VbdGR = int(BVolIV[-1])
        VbdGRStr1 = ">"
        VbdGRStr2 = ""
    else:
        VbdGRStr1 = "\\approx"
        VbdGRStr2 = "$--lime$ $line$"

    print "Results:"
    print "*******************************************************************************"
    print "Vdep_f1_log = ", Vdepf1Log, "+/-", Vdepf1LogEr
    print "Vdep_f2_log = ", Vdepf2Log, "+/-", Vdepf2LogEr
    print "Vdep_f1_1C = ", Vdepf1IC, "+/-", Vdepf1ICEr
    print "Vdep_f2_1C = ", Vdepf2IC, "+/-", Vdepf2ICEr

    print "Neff_f1_log = ", Neff1Log, "+/-", Neff1LogEr
    print "Neff_f2_log = ", Neff2Log, "+/-", Neff2LogEr
    print "Neff_f1_1C = ", Neff1IC, "+/-", Neff1ICEr
    print "Neff_f2_1C = ", Neff2IC, "+/-", Neff2ICEr

    print "Idep_f1_log = ", Idepf1Log
    print "Idep_f2_log = ", Idepf2Log
    print "Idep_f1_1C = ", Idepf1IC
    print "Idep_f2_1C = ", Idepf2IC

    print "V_break-down = ", Vbd

    print "IdepScaled_f1_log = ", Idepf1LogScaled
    print "IdepScaled_f2_log = ", Idepf2LogScaled
    print "IdepScaled_f1_1C = ", Idepf1ICScaled
    print "IdepScaled_f2_1C = ", Idepf2ICScaled
    print "*******************************************************************************"

    """"""""" CALCULATING DOPING PROFILE """""""""

    Neff1W = []
    Neff2W = []
    Wf1 = []
    Wf2 = []
    for i in range(1,len(BVolCV)):
        Neff1Width = abs((2./(qe*epsi0 * epsi * A**2))*((BVolCV[i]-BVolCV[i-1])/(1./Capf1[i]**2-1./Capf1[i-1]**2+1e-20)))*1e-17 # x10E11cm-3
        Widthf1 = df1Log * sqrt(BVolCV[i]/Vdepf1Log)
        if Widthf1 >= df1Log:
            break
        Neff1W.append(Neff1Width)
        Wf1.append(Widthf1*1e6) # in micron
    for i in range(1,len(BVolCV)):
        Neff2Width = abs((2./(qe*epsi0 * epsi * A**2))*((BVolCV[i]-BVolCV[i-1])/(1./Capf2[i]**2-1./Capf2[i-1]**2+1e-20)))*1e-17 # x10E11cm-3
        Widthf2 = df2Log * sqrt(BVolCV[i]/Vdepf2Log)
        if Widthf2 >= df2Log:
            break
        Neff2W.append(Neff2Width)
        Wf2.append(Widthf2*1e6) # in micron

    """"""""" CREATING LISTS TO PLOT FITTING RESULTS """""""""

    logslopef1x = np.linspace(min(logvolslopef1),BVolCV[Vdepf1LogIdx]+30,40000)
    logslopef1y = LogSlopefunc(logslopef1x, poptlogslopef1[0], poptlogslopef1[1])

    logflatf1x = np.linspace(BVolCV[Vdepf1LogIdx]-30,max(logvolflatf1),40000)
    if FlatFittingMethod1 == 0:
        logflatf1y = LogSlopefunc(logflatf1x, poptlogflatf1[0], poptlogflatf1[1])
    if FlatFittingMethod1 == 1:
        logflatf1y = logflatf1x * 0 + poptlogflatf1[0]

    loghorizf1y = np.linspace(Capf1[VmaxIdx]-Capf1[VmaxIdx]/1.1,Capf1[VminIdx]+10.*Capf1[VminIdx],40000)
    loghorizf1x = loghorizf1y * 0 + Vdepf1Log

    logslopef2x = np.linspace(min(logvolslopef2),BVolCV[Vdepf2LogIdx]+30,40000)
    logslopef2y = LogSlopefunc(logslopef2x, poptlogslopef2[0], poptlogslopef2[1])

    logflatf2x = np.linspace(BVolCV[Vdepf2LogIdx]-30,max(logvolflatf2),40000)
    if FlatFittingMethod1 == 0:
        logflatf2y = LogSlopefunc(logflatf2x, poptlogflatf2[0], poptlogflatf2[1])
    if FlatFittingMethod1 == 1:
        logflatf2y = logflatf2x * 0 + poptlogflatf2[0]

    loghorizf2y = np.linspace(Capf2[VmaxIdx]-Capf2[VmaxIdx]/1.1,Capf2[VminIdx]+10.*Capf2[VminIdx],40000)
    loghorizf2x = loghorizf2y * 0 + Vdepf2Log

    ICslopef1x = np.linspace(min(ICvolslopef1),BVolCV[Vdepf1ICIdx]+30,40000)
    ICslopef1y = Linearfunc(ICslopef1x, poptICslopef1[0], poptICslopef1[1])

    ICflatf1x = np.linspace(BVolCV[Vdepf1ICIdx]-30,max(ICvolflatf1),40000)
    if FlatFittingMethod2 == 0:
        ICflatf1y = Linearfunc(ICflatf1x, poptICflatf1[0], poptICflatf1[1])
    if FlatFittingMethod2 == 1:
        ICflatf1y = ICflatf1x * 0 + poptICflatf1[0]

    IChorizf1y = np.linspace(OneCapf1[VminIdx]-OneCapf1[VminIdx]/1.1,OneCapf1[VmaxIdx]+20.*OneCapf1[VmaxIdx],40000)
    IChorizf1x = IChorizf1y * 0 + Vdepf1IC

    ICslopef2x = np.linspace(min(ICvolslopef2),BVolCV[Vdepf2ICIdx]+30,40000)
    ICslopef2y = Linearfunc(ICslopef2x, poptICslopef2[0], poptICslopef2[1])

    ICflatf2x = np.linspace(BVolCV[Vdepf2ICIdx]-30,max(ICvolflatf2),40000)
    if FlatFittingMethod2 == 0:
        ICflatf2y = Linearfunc(ICflatf2x, poptICflatf2[0], poptICflatf2[1])
    if FlatFittingMethod2 == 1:
        ICflatf2y = ICflatf2x * 0 + poptICflatf2[0]

    IChorizf2y = np.linspace(OneCapf2[VminIdx]-OneCapf2[VminIdx]/1.1,OneCapf2[VmaxIdx]+20.*OneCapf2[VmaxIdx],40000)
    IChorizf2x = IChorizf2y * 0 + Vdepf2IC

    def CurrentUnit(Current):
        if Current < 1e-12:
            CurrentCoeff = 1e15
            print "Warning: Current is very small, make sure your measurement is OK!"
            CurrentUnit = "fA"
            CurrentUnit2 = "fA"
        if Current > 1e-12 and Current < 1e-9:
            CurrentCoeff = 1e12
            CurrentUnit = "pA"
            CurrentUnit2 = "pA"
        if Current > 1e-9 and Current < 1e-6:
            CurrentCoeff = 1e9
            CurrentUnit = "nA"
            CurrentUnit2 = "nA"
        if Current > 1e-6 and Current < 1e-3:
            CurrentCoeff = 1e6
            CurrentUnit = "$\mu$A"
            CurrentUnit2 = "uA"
        if Current > 1e-3 and Current < 1.:
            CurrentCoeff = 1e3
            CurrentUnit = "mA"
            CurrentUnit2 = "mA"
        if Current > 1.:
            CurrentCoeff = 1
            CurrentUnit = "A"
            CurrentUnit2 = "A"
        return CurrentCoeff, CurrentUnit, CurrentUnit2

    Idepf2LogCoeff, Idepf2LogUnit, Idepf2LogUnit2 = CurrentUnit(Idepf2Log)
    Idepf2LogScaledCoeff, Idepf2LogScaledUnit, Idepf2LogScaledUnit2 = CurrentUnit(Idepf2LogScaled)

    if EndCapf2Log < 1e-12:
        CapCoeff = 1e15
        CapUnit = "fF"
        CapUnit2 = "fF"
    if EndCapf2Log > 1e-12 and EndCapf2Log < 1e-9:
        CapCoeff = 1e12
        CapUnit = "pF"
        CapUnit2 = "pF"
    if EndCapf2Log > 1e-9 and EndCapf2Log < 1e-6:
        CapCoeff = 1e9
        CapUnit = "nF"
        CapUnit2 = "nF"
    if EndCapf2Log > 1e-6 and EndCapf2Log < 1e-3:
        CapCoeff = 1e6
        CapUnit = "$\mu$F"
        CapUnit2 = "uF"
    if EndCapf2Log > 1e-3 and EndCapf2Log < 1.:
        CapCoeff = 1e3
        CapUnit = "mF"
        CapUnit = "mF"
    if EndCapf2Log > 1.:
        CapCoeff = 1
        CapUnit = "F"
        CapUnit = "F"

    if df2Log < 1e-12:
        ThickCoeff = 1e15
        ThickUnit = "fm"
        ThickUnit2 = "fm"
    if df2Log > 1e-12 and df2Log < 1e-9:
        ThickCoeff = 1e12
        ThickUnit = "pm"
        ThickUnit2 = "pm"
    if df2Log > 1e-9 and df2Log < 1e-6:
        ThickCoeff = 1e9
        ThickUnit = "nm"
        ThickUnit2 = "nm"
    if df2Log > 1e-6 and df2Log < 1e-3:
        ThickCoeff = 1e6
        ThickUnit = "$\mu$m"
        ThickUnit2 = "um"
    if df2Log > 1e-3 and df2Log < 1.:
        ThickCoeff = 1e3
        ThickUnit = "mm"
        ThickUnit2 = "mm"
    if df2Log > 1.:
        ThickCoeff = 1
        ThickUnit = "m"
        ThickUnit2 = "m"

    if A < 1e-6:
        AreaCoeff = 1e12
        AreaUnit = "$\mu$m$^2$"
        AreaUnit2 = "um2"
    if A >= 1e-6 and A < 1e-4:
        AreaCoeff = 1e6
        AreaUnit = "mm$^2$"
        AreaUnit2 = "mm2"
    if A >= 1e-4 and A < 1.:
        AreaCoeff = 1e4
        AreaUnit = "cm$^2$"
        AreaUnit2 = "cm2"
    if A >= 1.:
        AreaCoeff = 1
        AreaUnit = "m$^2$"
        AreaUnit2 = "m2"
        
    LCur2 = [x*Idepf2LogCoeff for x in LCur]
    LCurScaled2 = [x*Idepf2LogCoeff for x in LCurScaled]
    GCur2 = [x*Idepf2LogCoeff for x in GCur]
    GCurScaled2 = [x*Idepf2LogCoeff for x in GCurScaled]
    BVolIV2 = BVolIV

    horizVbdy = np.linspace(LCur2[0]/10,10.*LCur2[-1],40000)
    horizVbdx = horizVbdy * 0 + Vbd
    horizVbdGRy = np.linspace(GCur2[0]/10,10.*GCur2[-1],40000)
    horizVbdGRx = horizVbdGRy * 0 + VbdGR

    if len(Annealing) > 17:
        Annealing = Annealing[:17]
        
    if CapMode == 0:
        MeasMode = "Parallel"
    if CapMode == 1:
        MeasMode = "Serial"
        
    if os.path.exists(PlotSavingLocation):
        pp = PdfPages('%(PlotSavingLocation)s%(SenName)s_CV_IV_%(Annealing)s_%(Environment)s_%(MeasMode)s.pdf' %locals())
    else:
        print "Error: Plots cannot be saved into given path. Please check if necessary directories are created.\nAborting..."
        os._exit(1)

    """"""""" LOG-LOG SCALED PLOTTING STARTS """""""""

    fig1 = plt.figure(num=1, figsize=(8.27, 11.69), dpi=100, facecolor='w', edgecolor='k')
    plt.suptitle("%(SenName)s_CV-IV_%(Annealing)s_%(Environment)s_%(MeasMode)s" %locals())

    gs1 = GridSpec(5, 1)
    gs1.update(left=0.12, right=0.60, wspace=0.05, hspace=0)

    gs2 = GridSpec(1, 1)
    gs2.update(left=0.70, right=0.98, hspace=0.03)

    axf1log = fig1.add_subplot(gs1[0:2, :])
    axf1log.loglog(BVolIV2,LCur2,'b.')
    axf1log.loglog(BVolIV2,LCurScaled2,'m.')#,color='0.35')
    if Vbd < BVolIV[-1]:
        axf1log.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    if VbdGR < BVolIV[-1]:
        axf1log.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
    axf1log.yaxis.tick_right()
    axf1log.yaxis.set_label_position("right")
    axf1log.set_ylabel("log of Current (%(Idepf2LogUnit)s)" %locals())
    axf1log.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf1log.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    axf1log.set_xticks([])
    axf1log.spines['right'].set_color('blue')
    axf1log.yaxis.label.set_color('blue')
    axf1log.tick_params(axis='y', colors='blue')

    axf1log2 = fig1.add_subplot(gs1[0:2, :], sharex=axf1log, frameon=False)
    axf1log2.loglog(BVolCV,Capf1,'k.')
    axf1log2.loglog(logslopef1x,logslopef1y,'r--', linewidth=2.0)
    axf1log2.loglog(logflatf1x,logflatf1y,'r--',linewidth=2.0)
    axf1log2.plot(loghorizf1x,loghorizf1y,'g--',linewidth=2.0)
    axf1log2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf1log2.set_ylim(min(Capf1[VminIdx:VmaxIdx])-min(Capf1[VminIdx:VmaxIdx])/5,max(Capf1[VminIdx:VmaxIdx])+max(Capf1[VminIdx:VmaxIdx])/6)
    axf1log2.set_xlabel("log of Voltage (V)")
    axf1log2.set_ylabel("log of Capacitance (F)")
    axf1log2.set_xticks([])

    if PlotDopingProfile == 1:
        axf1log3 = fig1.add_subplot(gs1[0:2, :], sharex=axf1log2, frameon=False)
        axf1log3.loglog(BVolIV2,GCur2, color='lime', linestyle='None', marker='.')
        if VbdGR < BVolIV[-1]:
            axf1log3.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
        axf1log3.yaxis.tick_right()
        axf1log3.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
        axf1log3.set_ylim(min([min(GCur2[VminIdx:VmaxIdx]),min(GCurScaled2[VminIdx:VmaxIdx])]),max([max(GCur2[VminIdx:VmaxIdx]),max(GCurScaled2[VminIdx:VmaxIdx])]))
        axf1log3.set_xticks([])
        axf1log3.set_yticks([])

    axf2log = fig1.add_subplot(gs1[2:4, :])
    axf2log.loglog(BVolIV2,LCur2,'b.')
    axf2log.loglog(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        axf2log.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    if VbdGR < BVolIV[-1]:
        axf2log.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
    axf2log.yaxis.tick_right()
    axf2log.yaxis.set_label_position("right")
    axf2log.set_ylabel("log of Current (%(Idepf2LogUnit)s)" %locals())
    axf2log.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf2log.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    axf2log.set_xticks([])
    axf2log.spines['right'].set_color('blue')
    axf2log.yaxis.label.set_color('blue')
    axf2log.tick_params(axis='y', colors='blue')

    axf2log2 = fig1.add_subplot(gs1[2:4, :], sharex=axf2log, frameon=False)
    axf2log2.loglog(BVolCV,Capf2,'k.')
    axf2log2.loglog(logslopef2x,logslopef2y,'r--', linewidth=2.0)
    axf2log2.loglog(logflatf2x,logflatf2y,'r--',linewidth=2.0)
    axf2log2.plot(loghorizf2x,loghorizf2y,'g--',linewidth=2.0)
    axf2log2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf2log2.set_ylim(min(Capf2[VminIdx:VmaxIdx])-min(Capf2[VminIdx:VmaxIdx])/5,max(Capf2[VminIdx:VmaxIdx])+max(Capf2[VminIdx:VmaxIdx])/6)
    axf2log2.set_xlabel("log of Voltage (V)")
    axf2log2.set_ylabel("log of Capacitance (F)")
    axf2log2.set_xticks([])

    axTlog = fig1.add_subplot(gs1[4:5, :])
    axTlog.semilogx(BVolCV,TempCV,'k.')
    axTlog.semilogx(BVolIV,TempIV,'r.')
    axTlog.set_ylabel("Temperature ($^0$C)")
    axTlog.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axTlog.set_ylim(min([min(TempCV[:-1]),min(TempIV[:-1])])-0.05,max([max(TempCV[:-1]),max(TempIV[:-1])])+0.05)
    axTlog.set_xlabel("log of Voltage (V)")
    axTlog.spines['right'].set_color('blue')

    if len(DewTempCV) > 0:
        axTlog2 = fig1.add_subplot(gs1[4:5, :], sharex=axTlog, frameon=False)
        axTlog2.semilogx(BVolCV,DewTempCV,'b.')
        axTlog2.semilogx(BVolIV,DewTempIV,'g.')
        axTlog2.yaxis.tick_right()
        axTlog2.yaxis.set_label_position("right")
        axTlog2.set_ylabel("Dew ($^0$C)")
        axTlog2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
        axTlog2.set_ylim(min([min(DewTempCV[:-1]),min(DewTempIV[:-1])])-0.05,max([max(DewTempCV[:-1]),max(DewTempIV[:-1])])+0.05)
        axTlog2.set_xlabel("log of Voltage (V)")
        axTlog2.spines['right'].set_color('blue')
        axTlog2.yaxis.label.set_color('blue')
        axTlog2.tick_params(axis='y', colors='blue')

    axboxlog = fig1.add_subplot(gs2[:, :])
    axboxlog.axis([0, 10, 0, 10])
    axboxlog.set_xticks([])
    axboxlog.set_yticks([])

    df1Log2 = round(df1Log*ThickCoeff,2)
    df2Log2 = round(df2Log*ThickCoeff,2)
    EndCapf1Log2 = round(EndCapf1Log*CapCoeff,3)
    EndCapf2Log2 = round(EndCapf2Log*CapCoeff,3)

    if Irradiation == 0:
        EndCapf1LogErr2 = round(EndCapf1LogErr*CapCoeff,3)
        EndCapf2LogErr2 = round(EndCapf2LogErr*CapCoeff,3)
        EndCapf1LogErr2 = "%(EndCapf1LogErr2)s" %locals()
        EndCapf2LogErr2 = "%(EndCapf2LogErr2)s" %locals()
        df1LogErr2 = round(df1LogErr*ThickCoeff,2)
        df2LogErr2 = round(df2LogErr*ThickCoeff,2)
        df1LogErr2 = "%(df1LogErr2)s" %locals()
        df2LogErr2 = "%(df2LogErr2)s" %locals()
    if Irradiation == 1:
        EndCapf1LogErr2 = ""
        EndCapf2LogErr2 = ""
        df1LogErr2 = ""
        df2LogErr2 = ""
        
    Cpostdepf1Log2 = round(Cpostdepf1Log * 1e12,3)
    Cpostdepf2Log2 = round(Cpostdepf2Log * 1e12,3)
    Cpostdepf1LogErr2 = round(Cpostdepf1LogErr * 1e12,3)
    Cpostdepf1LogErr2 = ""#"$\pm$%(Cpostdepf1LogErr2)s" %locals()
    Cpostdepf2LogErr2 = round(Cpostdepf2LogErr * 1e12,3)
    Cpostdepf2LogErr2 = ""#"$\pm$%(Cpostdepf2LogErr2)s" %locals()
    A2 = A * AreaCoeff
    v1 = int(BVolCV[VminIdx])
    v2 = int(BVolCV[VmaxIdx])
    vrange = "%(v1)s - %(v2)s" %locals()
    Vdepf1Log2 = round(Vdepf1Log,2)
    Vdepf1LogEr2 = round(Vdepf1LogEr,2)
    Vdepf2Log2 = round(Vdepf2Log,2)
    Vdepf2LogEr2 = round(Vdepf2LogEr,2)

    Neff1Log2 = round(Neff1Log,2)
    Neff1LogEr2 = round(Neff1LogEr,2)
    Neff2Log2 = round(Neff2Log,2)
    Neff2LogEr2 = round(Neff2LogEr,2)

    Idepf1Log2 = round(Idepf1Log*Idepf2LogCoeff,2)
    Idepf2Log2 = round(Idepf2Log*Idepf2LogCoeff,2)
    Idepf1LogScaled2 = round(Idepf1LogScaled*Idepf2LogScaledCoeff,2)
    Idepf2LogScaled2 = round(Idepf2LogScaled*Idepf2LogScaledCoeff,2)
    Ipostdepf1Log2 = round(Ipostdepf1Log*Idepf2LogCoeff,2)
    Ipostdepf2Log2 = round(Ipostdepf2Log*Idepf2LogCoeff,2)
    Ipostdepf1LogScaled2 = round(Ipostdepf1LogScaled*Idepf2LogScaledCoeff,2)
    Ipostdepf2LogScaled2 = round(Ipostdepf2LogScaled*Idepf2LogScaledCoeff,2)
    ScaleTemp2 = int(ScaleTemp)
    import ntpath
    CVFileLocation2 = ntpath.basename(CVFileLocation)
    CVFileLocation21 = CVFileLocation2[:39]
    CVFileLocation22 = CVFileLocation2[39:]
    IVFileLocation2 = ntpath.basename(IVFileLocation)
    IVFileLocation21 = IVFileLocation2[:39]
    IVFileLocation22 = IVFileLocation2[39:]
    if FlatFittingMethod1 == 0:
        Fitting1 = "Linear-Linear"
    if FlatFittingMethod1 == 1:
        Fitting1 = "Linear-Horizontal"

    tlog1 = "Frequency: %(freq1)s (Hz)\nFitting: %(Fitting1)s\nRange: %(vrange)s (V)\nAnnealing:%(Annealing)s\n"\
            "Thickness: %(df1Log2)s$\pm$%(df1LogErr2)s (%(ThickUnit)s)\nActive Area: %(A2)s (%(AreaUnit)s)\nCapacitance:"\
            " %(EndCapf1Log2)s$\pm$%(EndCapf1LogErr2)s (%(CapUnit)s)\nV$_{dep}$: "\
            "%(Vdepf1Log2)s $\pm$ %(Vdepf1LogEr2)s (V)\nV$_{bd}$: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV$_{bd}$(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN$_{eff}$: %(Neff1Log2)s $\pm$ %(Neff1LogEr2)s "\
            "(10$^{11}$cm$^{-3}$)\nI$_{V_{dep}}$: %(Idepf1Log2)s (%(Idepf2LogUnit)s)\nI$_{V_{dep}}$(%(ScaleTemp2)s$^0$C): "\
            "%(Idepf1LogScaled2)s (%(Idepf2LogScaledUnit)s)\nI$_{1.2V_{dep}}$: %(Ipostdepf1Log2)s (%(Idepf2LogUnit)s)\nI$_{1.2V_{dep}}"\
            "$(%(ScaleTemp2)s$^0$C): %(Ipostdepf1LogScaled2)s (%(Idepf2LogScaledUnit)s)\nC$_{1.2V_{dep}}$: %(Cpostdepf1Log2)s%(Cpostdepf1LogErr2)s (%(CapUnit)s)"\
            "\n\n_Configuration Parameters_\n  slopef1min: %(slopef1min)s\n  slopef1max: "\
            "%(slopef1max)s\n  flatf1min: %(flatf1min)s\n  flatf1max: %(flatf1max)s" % locals()

    tlog2 = "Frequency: %(freq2)s (Hz)\nFitting: %(Fitting1)s\nRange: %(vrange)s (V)\nAnnealing:%(Annealing)s\n"\
            "Thickness: %(df2Log2)s$\pm$%(df2LogErr2)s (%(ThickUnit)s)\nActive Area: %(A2)s (%(AreaUnit)s)\nCapacitance: %(EndCapf2Log2)s$\pm$%(EndCapf2LogErr2)s (%(CapUnit)s)\nV$_{dep}$: "\
            "%(Vdepf2Log2)s $\pm$ %(Vdepf2LogEr2)s (V)\nV$_{bd}$: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV$_{bd}$(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN$_{eff}$: %(Neff2Log2)s $\pm$ %(Neff2LogEr2)s "\
            "(10$^{11}$cm$^{-3}$)\nI$_{V_{dep}}$: %(Idepf2Log2)s (%(Idepf2LogUnit)s)\nI$_{V_{dep}}$(%(ScaleTemp2)s$^0$C): "\
            "%(Idepf2LogScaled2)s (%(Idepf2LogScaledUnit)s)\nI$_{1.2V_{dep}}$: %(Ipostdepf2Log2)s (%(Idepf2LogUnit)s)\nI$_{1.2V_{dep}}"\
            "$(%(ScaleTemp2)s$^0$C): %(Ipostdepf2LogScaled2)s (%(Idepf2LogScaledUnit)s)\nC$_{1.2V_{dep}}$: %(Cpostdepf2Log2)s%(Cpostdepf2LogErr2)s (%(CapUnit)s)"\
            "\n\n_Configuration Parameters_\n  slopef2min: %(slopef2min)s\n  slopef2max: "\
            "%(slopef2max)s\n  flatf2min: %(flatf2min)s\n  flatf2max: %(flatf2max)s" % locals()
    tlog3 = "CV Filename:"
    tlog4 = "%(CVFileLocation21)s\n%(CVFileLocation22)s" % locals()
    tlog5 = "IV Filename:"
    tlog6 = "%(IVFileLocation21)s\n%(IVFileLocation22)s" % locals()
    tlog7 = "Note: Cyan is the scaled current on IV curves.\n_Legend for Temperature Plots_\n  Black: CV temperature\n  Blue: CV dew point temperature\n  Red: IV temperature\n  Green: IV dew point temperature"
    axboxlog.text(0.2, 9.9, tlog1, fontsize=10, ha='left', va='top')
    axboxlog.text(0.2, 5.9, tlog2, fontsize=10, ha='left', va='top')
    axboxlog.text(0.2, 1.9, tlog3, fontsize=10, ha='left', va='top')
    axboxlog.text(0.2, 1.75, tlog4, fontsize=7, ha='left', va='top')
    axboxlog.text(0.2, 1.5, tlog5, fontsize=10, ha='left', va='top')
    axboxlog.text(0.2, 1.35, tlog6, fontsize=7, ha='left', va='top')
    axboxlog.text(0.2, 0.8, tlog7, fontsize=7, ha='left', va='top')

    pp.savefig(fig1,papertype='a4',orientation='portrait')
    close()

    """"""""" LOG-LOG SCALED PLOTTING ENDS """""""""

    """"""""" 1/C2 vs. BIAS VOL. PLOTTING STARTS """""""""

    fig2 = plt.figure(num=2, figsize=(8.27, 11.69), dpi=100, facecolor='w', edgecolor='k')
    plt.suptitle("%(SenName)s_CV-IV_%(Annealing)s_%(Environment)s_%(MeasMode)s" %locals())

    gs3 = GridSpec(5, 1)
    gs3.update(left=0.12, right=0.60, wspace=0.05, hspace=0)

    gs4 = GridSpec(1, 1)
    gs4.update(left=0.70, right=0.98, hspace=0.03)

    axf1IC = fig2.add_subplot(gs3[0:2, :])
    axf1IC.plot(BVolIV2,LCur2,'b.')
    axf1IC.plot(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        axf1IC.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    if VbdGR < BVolIV[-1]:
        axf1IC.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
    axf1IC.yaxis.tick_right()
    axf1IC.yaxis.set_label_position("right")
    axf1IC.set_ylabel("Current (%(Idepf2LogUnit)s)" %locals())
    axf1IC.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf1IC.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    axf1IC.set_xticks([])
    axf1IC.spines['right'].set_color('blue')
    axf1IC.yaxis.label.set_color('blue')
    axf1IC.tick_params(axis='y', colors='blue')

    axf1IC2 = fig2.add_subplot(gs3[0:2, :], sharex=axf1IC, frameon=False)
    axf1IC2.plot(BVolCV,OneCapf1,'k.')
    axf1IC2.plot(ICslopef1x,ICslopef1y,'r--', linewidth=2.0)
    axf1IC2.plot(ICflatf1x,ICflatf1y,'r--',linewidth=2.0)
    axf1IC2.plot(IChorizf1x,IChorizf1y,'g--',linewidth=2.0)
    axf1IC2.set_ylabel("1/C$^2$ (x10$^{21}$ F$^{-2}$)")
    axf1IC2.set_xlabel("Voltage (V)")
    axf1IC2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf1IC2.set_ylim(min(OneCapf1[VminIdx:VmaxIdx])-min(OneCapf1[VminIdx:VmaxIdx])/6,max(OneCapf1[VminIdx:VmaxIdx])+max(OneCapf1[VminIdx:VmaxIdx])/5)
    axf1IC2.set_xticks([])

    axf2IC = fig2.add_subplot(gs3[2:4, :])
    axf2IC.plot(BVolIV2,LCur2,'b.')
    axf2IC.plot(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        axf2IC.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    if VbdGR < BVolIV[-1]:
        axf2IC.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
    axf2IC.yaxis.tick_right()
    axf2IC.yaxis.set_label_position("right")
    axf2IC.set_ylabel("Current (%(Idepf2LogUnit)s)" %locals())
    axf2IC.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf2IC.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    axf2IC.set_xticks([])
    axf2IC.spines['right'].set_color('blue')
    axf2IC.yaxis.label.set_color('blue')
    axf2IC.tick_params(axis='y', colors='blue')

    axf2IC2 = fig2.add_subplot(gs3[2:4, :], sharex=axf2IC, frameon=False)
    axf2IC2.plot(BVolCV,OneCapf2,'k.')
    axf2IC2.plot(ICslopef2x,ICslopef2y,'r--', linewidth=2.0)
    axf2IC2.plot(ICflatf2x,ICflatf2y,'r--',linewidth=2.0)
    axf2IC2.plot(IChorizf2x,IChorizf2y,'g--',linewidth=2.0)
    axf2IC2.set_ylabel("1/C$^2$ (x10$^{21}$ F$^{-2}$)")
    axf2IC2.set_xlabel("Voltage (V)")
    axf2IC2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axf2IC2.set_ylim(min(OneCapf2[VminIdx:VmaxIdx])-min(OneCapf2[VminIdx:VmaxIdx])/6,max(OneCapf2[VminIdx:VmaxIdx])+max(OneCapf2[VminIdx:VmaxIdx])/5)
    axf2IC2.set_xticks([])

    axTIC = fig2.add_subplot(gs3[4:5, :])
    axTIC.plot(BVolCV,TempCV,'k.',label="Temperature CV")
    axTIC.plot(BVolIV,TempIV,'r.',label="Temperature IV")
    axTIC.set_ylabel("Temperature ($^0$C)")
    axTIC.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    axTIC.set_ylim(min([min(TempCV[:-1]),min(TempIV[:-1])])-0.05,max([max(TempCV[:-1]),max(TempIV[:-1])])+0.05)
    axTIC.set_xlabel("Voltage (V)")
    axTIC.spines['right'].set_color('blue')

    if len(DewTempCV) > 0:
        axTIC2 = fig2.add_subplot(gs3[4:5, :], sharex=axTIC, frameon=False)
        axTIC2.plot(BVolCV,DewTempCV,'b.',label="Dew CV")
        axTIC2.plot(BVolIV,DewTempIV,'g.',label="Dew IV")
        axTIC2.yaxis.tick_right()
        axTIC2.yaxis.set_label_position("right")
        axTIC2.set_ylabel("Dew ($^0$C)")
        axTIC2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
        axTIC2.set_ylim(min([min(DewTempCV[:-1]),min(DewTempIV[:-1])])-0.05,max([max(DewTempCV[:-1]),max(DewTempIV[:-1])])+0.05)
        axTIC2.set_xlabel("Voltage (V)")
        axTIC2.spines['right'].set_color('blue')
        axTIC2.yaxis.label.set_color('blue')
        axTIC2.tick_params(axis='y', colors='blue')

    axboxIC = fig2.add_subplot(gs4[:, :])
    axboxIC.axis([0, 10, 0, 10])
    axboxIC.set_xticks([])
    axboxIC.set_yticks([])

    df1IC2 = round(df1IC*ThickCoeff,2)
    df2IC2 = round(df2IC*ThickCoeff,2)
    EndCapf1IC2 = round(EndCapf1IC*CapCoeff,3)
    EndCapf2IC2 = round(EndCapf2IC*CapCoeff,3)

    if Irradiation == 0:
        EndCapf1ICErr2 = round(EndCapf1ICErr*CapCoeff,3)
        EndCapf2ICErr2 = round(EndCapf2ICErr*CapCoeff,3)
        EndCapf1ICErr2 = "%(EndCapf1ICErr2)s" %locals()
        EndCapf2ICErr2 = "%(EndCapf2ICErr2)s" %locals()
        df1ICErr2 = round(df1ICErr*ThickCoeff,2)
        df2ICErr2 = round(df2ICErr*ThickCoeff,2)
        df1ICErr2 = "%(df1ICErr2)s" %locals()
        df2ICErr2 = "%(df2ICErr2)s" %locals()
    if Irradiation == 1:
        EndCapf1ICErr2 = ""
        EndCapf2ICErr2 = ""
        df1ICErr2 = ""
        df2ICErr2 = ""

    Cpostdepf1IC2 = round(Cpostdepf1IC * 1e12,3)
    Cpostdepf2IC2 = round(Cpostdepf2IC * 1e12,3)
    Vdepf1IC2 = round(Vdepf1IC,2)
    Vdepf1ICEr2 = round(Vdepf1ICEr,2)
    Vdepf2IC2 = round(Vdepf2IC,2)
    Vdepf2ICEr2 = round(Vdepf2ICEr,2)
    Neff1IC2 = round(Neff1IC,2)
    Neff1ICEr2 = round(Neff1ICEr,2)
    Neff2IC2 = round(Neff2IC,2)
    Neff2ICEr2 = round(Neff2ICEr,2)
    Idepf1IC2 = round(Idepf1IC*Idepf2LogCoeff,2)
    Idepf2IC2 = round(Idepf2IC*Idepf2LogCoeff,2)
    Idepf1ICScaled2 = round(Idepf1ICScaled*Idepf2LogScaledCoeff,2)
    Idepf2ICScaled2 = round(Idepf2ICScaled*Idepf2LogScaledCoeff,2)
    Ipostdepf1IC2 = round(Ipostdepf1IC*Idepf2LogCoeff,2)
    Ipostdepf2IC2 = round(Ipostdepf2IC*Idepf2LogCoeff,2)
    Ipostdepf1ICScaled2 = round(Ipostdepf1ICScaled*Idepf2LogScaledCoeff,2)
    Ipostdepf2ICScaled2 = round(Ipostdepf2ICScaled*Idepf2LogScaledCoeff,2)
    if slopef1min_1C < 0:
        slopef1min2,slopef1max2, flatf1min2, flatf1max2, slopef2min2,slopef2max2, flatf2min2, flatf2max2 = slopef1min,slopef1max, flatf1min, flatf1max, slopef2min,slopef2max, flatf2min, flatf2max
    else:
        slopef1min2,slopef1max2, flatf1min2, flatf1max2, slopef2min2,slopef2max2, flatf2min2, flatf2max2 = slopef1min_1C,slopef1max_1C, flatf1min_1C, flatf1max_1C, slopef2min_1C,slopef2max_1C, flatf2min_1C, flatf2max_1C
    if FlatFittingMethod2 == 0:
        Fitting2 = "Linear-Linear"
    if FlatFittingMethod2 == 1:
        Fitting2 = "Linear-Horizontal"

    tIC1 = "Frequency: %(freq1)s (Hz)\nFitting: %(Fitting2)s\nRange: %(vrange)s (V)\nAnnealing:%(Annealing)s\n"\
           "Thickness: %(df1IC2)s$\pm$%(df1ICErr2)s (%(ThickUnit)s)\nActive Area: %(A2)s (%(AreaUnit)s)\nCapacitance: %(EndCapf1IC2)s$\pm$%(EndCapf1ICErr2)s (%(CapUnit)s)\nV$_{dep}$: "\
           "%(Vdepf1IC2)s $\pm$ %(Vdepf1ICEr2)s (V)\nV$_{bd}$: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV$_{bd}$(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN$_{eff}$: %(Neff1IC2)s $\pm$ %(Neff1ICEr2)s "\
            "(10$^{11}$cm$^{-3}$)\nI$_{V_{dep}}$: %(Idepf1IC2)s (%(Idepf2LogUnit)s)\nI$_{V_{dep}}$(%(ScaleTemp2)s$^0$C): "\
            "%(Idepf1ICScaled2)s (%(Idepf2LogScaledUnit)s)\nI$_{1.2V_{dep}}$: %(Ipostdepf1IC2)s (%(Idepf2LogUnit)s)\nI$_{1.2V_{dep}}"\
            "$(%(ScaleTemp2)s$^0$C): %(Ipostdepf1ICScaled2)s (%(Idepf2LogScaledUnit)s)\nC$_{1.2V_{dep}}$: %(Cpostdepf1IC2)s (%(CapUnit)s)"\
            "\n\n_Configuration Parameters_\n  slopef1min: %(slopef1min2)s\n  slopef1max: "\
            "%(slopef1max2)s\n  flatf1min: %(flatf1min2)s\n  flatf1max: %(flatf1max2)s" % locals()

    tIC2 = "Frequency: %(freq2)s (Hz)\nFitting: %(Fitting2)s\nRange: %(vrange)s (V)\nAnnealing:%(Annealing)s\n"\
           "Thickness: %(df2IC2)s$\pm$%(df2ICErr2)s (%(ThickUnit)s)\nActive Area: %(A2)s (%(AreaUnit)s)\nCapacitance: %(EndCapf2IC2)s$\pm$%(EndCapf2ICErr2)s (%(CapUnit)s)\nV$_{dep}$: "\
           "%(Vdepf2IC2)s $\pm$ %(Vdepf2ICEr2)s (V)\nV$_{bd}$: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV$_{bd}$(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN$_{eff}$: %(Neff2IC2)s $\pm$ %(Neff2ICEr2)s "\
            "(10$^{11}$cm$^{-3}$)\nI$_{V_{dep}}$: %(Idepf2IC2)s (%(Idepf2LogUnit)s)\nI$_{V_{dep}}$(%(ScaleTemp2)s$^0$C): "\
            "%(Idepf2ICScaled2)s (%(Idepf2LogScaledUnit)s)\nI$_{1.2V_{dep}}$: %(Ipostdepf2IC2)s (%(Idepf2LogUnit)s)\nI$_{1.2V_{dep}}"\
            "$(%(ScaleTemp2)s$^0$C): %(Ipostdepf2ICScaled2)s (%(Idepf2LogScaledUnit)s)\nC$_{1.2V_{dep}}$: %(Cpostdepf2IC2)s (%(CapUnit)s)"\
            "\n\n_Configuration Parameters_\n  slopef2min: %(slopef2min2)s\n  slopef2max: "\
            "%(slopef2max2)s\n  flatf2min: %(flatf2min2)s\n  flatf2max: %(flatf2max2)s" % locals()
    tIC3 = "CV Filename:"
    tIC4 = "%(CVFileLocation21)s\n%(CVFileLocation22)s" % locals()
    tIC5 = "IV Filename:"
    tIC6 = "%(IVFileLocation21)s\n%(IVFileLocation22)s" % locals()
    tIC7 = "Note: Cyan is the scaled current on IV curves.\n_Legend for Temperature Plots_\n  Black: CV temperature\n  Blue: CV dew point temperature\n  Red: IV temperature\n  Green: IV dew point temperature"
    axboxIC.text(0.2, 9.9, tIC1, fontsize=10, ha='left', va='top')
    axboxIC.text(0.2, 5.9, tIC2, fontsize=10, ha='left', va='top')
    axboxIC.text(0.2, 1.9, tIC3, fontsize=10, ha='left', va='top')
    axboxIC.text(0.2, 1.75, tIC4, fontsize=7, ha='left', va='top')
    axboxIC.text(0.2, 1.5, tIC5, fontsize=10, ha='left', va='top')
    axboxIC.text(0.2, 1.35, tIC6, fontsize=7, ha='left', va='top')
    axboxIC.text(0.2, 0.8, tIC7, fontsize=7, ha='left', va='top')

    pp.savefig(fig2,papertype='a4',orientation='portrait')
    close()
    pp.close()
    print "\nPlots are succesfully saved in /%(PlotSavingLocation)s%(SenName)s_CV_IV_%(Annealing)s_%(Environment)s_%(MeasMode)s.pdf" %locals()

    """"""""" 1/C2 vs. BIAS VOL. PLOTTING ENDS """""""""

    """"""""" WRITING RESULTS IN A TEXT FILE WITH THE SAME NAME AS THE PDF FILE """""""""

    with open('%(PlotSavingLocation)s%(SenName)s_CV_IV_%(Annealing)s_%(Environment)s_%(MeasMode)s.dat' %locals(), 'w') as textfile:
        textfile.write("Frequency: %(freq1)s (Hz)\nFitting: %(Fitting1)s-loglog\nRange: %(vrange)s (V)\nAnnealing: %(Annealing)s\n"\
            "Thickness: %(df1Log2)s pm %(df1LogErr2)s (%(ThickUnit2)s)\nActive Area: %(A2)s (%(AreaUnit2)s)\nCapacitance:"\
            " %(EndCapf1Log2)s pm %(EndCapf1LogErr2)s (%(CapUnit2)s)\nV_dep: "
            "%(Vdepf1Log2)s pm %(Vdepf1LogEr2)s (V)\nV_bd: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV_bd(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN_eff: %(Neff1Log2)s pm %(Neff1LogEr2)s "\
            "(10E11 cm-3)\nI_(V_dep): %(Idepf1Log2)s (%(Idepf2LogUnit2)s)\nI_(V_dep)(%(ScaleTemp2)sC): "\
            "%(Idepf1LogScaled2)s (%(Idepf2LogScaledUnit2)s)\nI_(1.2V_dep): %(Ipostdepf1Log2)s (%(Idepf2LogUnit2)s)\nI_(1.2V_dep)"\
            "(%(ScaleTemp2)sC): %(Ipostdepf1LogScaled2)s (%(Idepf2LogScaledUnit2)s)\nC_(1.2V_dep): %(Cpostdepf1Log2)s%(Cpostdepf1LogErr2)s (%(CapUnit2)s)" % locals())
        textfile.write('\n')
        textfile.write('\n')
        textfile.write("Frequency: %(freq2)s (Hz)\nFitting: %(Fitting1)s-loglog\nRange: %(vrange)s (V)\nAnnealing: %(Annealing)s\n"\
            "Thickness: %(df2Log2)s pm %(df2LogErr2)s (%(ThickUnit2)s)\nActive Area: %(A2)s (%(AreaUnit2)s)\nCapacitance: %(EndCapf2Log2)s pm %(EndCapf2LogErr2)s (%(CapUnit2)s)\nV_dep: "\
            "%(Vdepf2Log2)s pm %(Vdepf2LogEr2)s (V)\nV_bd: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV_bd(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN_eff: %(Neff2Log2)s pm %(Neff2LogEr2)s "\
            "(10E11 cm-3)\nI_(V_dep): %(Idepf2Log2)s (%(Idepf2LogUnit2)s)\nI_(V_dep)(%(ScaleTemp2)sC): "\
            "%(Idepf2LogScaled2)s (%(Idepf2LogScaledUnit2)s)\nI_(1.2V_dep): %(Ipostdepf2Log2)s (%(Idepf2LogUnit2)s)\nI_(1.2V_dep)"\
            "(%(ScaleTemp2)sC): %(Ipostdepf2LogScaled2)s (%(Idepf2LogScaledUnit2)s)\nC_(1.2V_dep): %(Cpostdepf2Log2)s%(Cpostdepf2LogErr2)s (%(CapUnit2)s)" % locals())
        textfile.write('\n')
        textfile.write('\n')
        textfile.write("Frequency: %(freq1)s (Hz)\nFitting: %(Fitting2)s-1C2\nRange: %(vrange)s (V)\nAnnealing: %(Annealing)s\n"\
           "Thickness: %(df1IC2)s pm %(df1ICErr2)s (%(ThickUnit2)s)\nActive Area: %(A2)s (%(AreaUnit2)s)\nCapacitance: %(EndCapf1IC2)s pm %(EndCapf1ICErr2)s (%(CapUnit2)s)\nV_dep: "\
           "%(Vdepf1IC2)s pm %(Vdepf1ICEr2)s (V)\nV_bd: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV_bd(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN_eff: %(Neff1IC2)s pm %(Neff1ICEr2)s "\
            "(10E11 cm-3)\nI_(V_dep): %(Idepf1IC2)s (%(Idepf2LogUnit2)s)\nI_(V_dep)(%(ScaleTemp2)sC): "\
            "%(Idepf1ICScaled2)s (%(Idepf2LogScaledUnit2)s)\nI_(1.2V_dep): %(Ipostdepf1IC2)s (%(Idepf2LogUnit2)s)\nI_(1.2V_dep)"\
            "(%(ScaleTemp2)sC): %(Ipostdepf1ICScaled2)s (%(Idepf2LogScaledUnit2)s)\nC_(1.2V_dep): %(Cpostdepf1IC2)s (%(CapUnit2)s)" % locals())
        textfile.write('\n')
        textfile.write('\n')
        textfile.write("Frequency: %(freq2)s (Hz)\nFitting: %(Fitting2)s-1C2\nRange: %(vrange)s (V)\nAnnealing: %(Annealing)s\n"\
           "Thickness: %(df2IC2)s pm %(df2ICErr2)s (%(ThickUnit2)s)\nActive Area: %(A2)s (%(AreaUnit2)s)\nCapacitance: %(EndCapf2IC2)s pm %(EndCapf2ICErr2)s (%(CapUnit2)s)\nV_dep: "\
           "%(Vdepf2IC2)s pm %(Vdepf2ICEr2)s (V)\nV_bd: $%(VbdStr1)s$ %(Vbd)s (V)%(VbdStr2)s\nV_bd(GR): $%(VbdGRStr1)s$ %(VbdGR)s (V)%(VbdGRStr2)s\nN_eff: %(Neff2IC2)s pm %(Neff2ICEr2)s "\
            "(10E11 cm-3)\nI_(V_dep): %(Idepf2IC2)s (%(Idepf2LogUnit2)s)\nI_(V_dep)(%(ScaleTemp2)sC): "\
            "%(Idepf2ICScaled2)s (%(Idepf2LogScaledUnit2)s)\nI_(1.2V_dep): %(Ipostdepf2IC2)s (%(Idepf2LogUnit2)s)\nI_(1.2V_dep)"\
            "(%(ScaleTemp2)sC): %(Ipostdepf2ICScaled2)s (%(Idepf2LogScaledUnit2)s)\nC_(1.2V_dep): %(Cpostdepf2IC2)s (%(CapUnit2)s)" % locals())

    """"""""" PLOTTING POP-UP FIGURES TO FINE-TUNE FITTING RANGES """""""""

    fig3 = plt.figure(3)
    ax1 = fig3.add_subplot(111)
    ax1.plot(BVolIV2,LCur2,'b.')
    ax1.plot(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        ax1.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.set_ylabel("Current (%(Idepf2LogUnit)s)" %locals())
    ax1.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax1.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    ax1.spines['right'].set_color('blue')
    ax1.yaxis.label.set_color('blue')
    ax1.tick_params(axis='y', colors='blue')

    ax2 = fig3.add_subplot(111, sharex=ax1, frameon=False)
    ax2.plot(BVolCV,OneCapf1,'k.')
    ax2.plot(ICslopef1x,ICslopef1y,'r--', linewidth=2.0)
    ax2.plot(ICflatf1x,ICflatf1y,'r--',linewidth=2.0)
    ax2.plot(IChorizf1x,IChorizf1y,'g--',linewidth=2.0)
    ax2.set_ylabel("1/C$^2$ (x10$^{21}$ F$^{-2}$)")
    ax2.set_xlabel("Voltage (V)")
    ax2.set_title("%(SenName)s, %(freq1)sHz, %(Annealing)s, %(Environment)s" %locals())
    ax2.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax2.set_ylim(min(OneCapf1[VminIdx:VmaxIdx])-min(OneCapf1[VminIdx:VmaxIdx])/6,max(OneCapf1[VminIdx:VmaxIdx])+max(OneCapf1[VminIdx:VmaxIdx])/6)

    fig4 = plt.figure(4)
    ax3 = fig4.add_subplot(111)
    ax3.plot(BVolIV2,LCur2,'b.')
    ax3.plot(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        ax3.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    ax3.set_ylabel("Current (%(Idepf2LogUnit)s)" %locals())
    ax3.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax3.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    ax3.spines['right'].set_color('blue')
    ax3.yaxis.label.set_color('blue')
    ax3.tick_params(axis='y', colors='blue')

    ax4 = fig4.add_subplot(111, sharex=ax3, frameon=False)
    ax4.plot(BVolCV,OneCapf2,'k.')
    ax4.plot(ICslopef2x,ICslopef2y,'r--', linewidth=2.0)
    ax4.plot(ICflatf2x,ICflatf2y,'r--',linewidth=2.0)
    ax4.plot(IChorizf2x,IChorizf2y,'g--',linewidth=2.0)
    ax4.set_ylabel("1/C$^2$ (x10$^{21}$ F$^{-2}$)")
    ax4.set_xlabel("Voltage (V)")
    ax4.set_title("%(SenName)s, %(freq2)sHz, %(Annealing)s, %(Environment)s" %locals())
    ax4.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax4.set_ylim(min(OneCapf2[VminIdx:VmaxIdx])-min(OneCapf2[VminIdx:VmaxIdx])/6,max(OneCapf2[VminIdx:VmaxIdx])+max(OneCapf2[VminIdx:VmaxIdx])/6)

    ax9 = fig4.add_subplot(111, sharex=ax4, frameon=False)
    ax9.plot(BVolIV2,GCur2, color='lime', linestyle='None', marker='.')
    if VbdGR < BVolIV[-1]:
        ax9.plot(horizVbdGRx,horizVbdGRy, color='lime', linestyle='dashed', linewidth=2.0)
    ax9.yaxis.tick_right()
    ax9.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax9.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(GCur2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(GCur2[VminIdx:VmaxIdx])]))
    ax9.set_yticks([])

    fig5 = plt.figure(5)
    ax5 = fig5.add_subplot(111)
    ax5.loglog(BVolIV2,LCur2,'b.')
    ax5.loglog(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        ax5.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    ax5.set_ylabel("log of Current (%(Idepf2LogUnit)s)" %locals())
    ax5.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax5.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    ax5.spines['right'].set_color('blue')
    ax5.yaxis.label.set_color('blue')
    ax5.tick_params(axis='y', colors='blue')

    ax6 = fig5.add_subplot(111, sharex=ax5, frameon=False)
    ax6.loglog(BVolCV,Capf1,'k.')
    ax6.loglog(logslopef1x,logslopef1y,'r--', linewidth=2.0)
    ax6.loglog(logflatf1x,logflatf1y,'r--',linewidth=2.0)
    ax6.loglog(loghorizf1x,loghorizf1y,'g--',linewidth=2.0)
    ax6.set_ylabel("log of Capacitance (F)")
    ax6.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax6.set_ylim(min(Capf1[VminIdx:VmaxIdx])-min(Capf1[VminIdx:VmaxIdx])/5,max(Capf1[VminIdx:VmaxIdx])+max(Capf1[VminIdx:VmaxIdx])/6)
    ax6.set_xlabel("log of Voltage (V)")
    ax6.set_title("%(SenName)s, %(freq1)sHz, %(Annealing)s, %(Environment)s" %locals())

    fig6 = plt.figure(6)
    ax7 = fig6.add_subplot(111)
    ax7.loglog(BVolIV2,LCur2,'b.')
    ax7.loglog(BVolIV2,LCurScaled2,'m.')
    if Vbd < BVolIV[-1]:
        ax7.plot(horizVbdx,horizVbdy,'y--',linewidth=2.0)
    ax7.yaxis.tick_right()
    ax7.yaxis.set_label_position("right")
    ax7.set_ylabel("log of Current (%(Idepf2LogUnit)s)" %locals())
    ax7.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax7.set_ylim(min([min(LCur2[VminIdx:VmaxIdx]),min(LCurScaled2[VminIdx:VmaxIdx])]),max([max(LCur2[VminIdx:VmaxIdx]),max(LCurScaled2[VminIdx:VmaxIdx])]))
    ax7.spines['right'].set_color('blue')
    ax7.yaxis.label.set_color('blue')
    ax7.tick_params(axis='y', colors='blue')

    ax8 = fig6.add_subplot(111, sharex=ax7, frameon=False)
    ax8.loglog(BVolCV,Capf2,'k.')
    ax8.loglog(logslopef2x,logslopef2y,'r--', linewidth=2.0)
    ax8.loglog(logflatf2x,logflatf2y,'r--',linewidth=2.0)
    ax8.loglog(loghorizf2x,loghorizf2y,'g--',linewidth=2.0)
    ax8.set_ylabel("log of Capacitance (F)")
    ax8.set_xlim(BVolCV[VminIdx],BVolCV[VmaxIdx])
    ax8.set_ylim(min(Capf2[VminIdx:VmaxIdx])-min(Capf2[VminIdx:VmaxIdx])/5,max(Capf2[VminIdx:VmaxIdx])+max(Capf2[VminIdx:VmaxIdx])/6)
    ax8.set_xlabel("log of Voltage (V)")
    ax8.set_title("%(SenName)s, %(freq2)sHz, %(Annealing)s, %(Environment)s" %locals())

    if SavePopUpPlots == 1:
        fig3.savefig('%(PlotSavingLocation)s%(SenName)s_%(freq1)sHz_1-C2_%(Annealing)s_%(Environment)s_%(MeasMode)s.png' %locals())
        fig4.savefig('%(PlotSavingLocation)s%(SenName)s_%(freq2)sHz_1-C2_%(Annealing)s_%(Environment)s_%(MeasMode)s.png' %locals())
        fig5.savefig('%(PlotSavingLocation)s%(SenName)s_%(freq1)sHz_log-log_%(Annealing)s_%(Environment)s_%(MeasMode)s.png' %locals())
        fig6.savefig('%(PlotSavingLocation)s%(SenName)s_%(freq2)sHz_log-log_%(Annealing)s_%(Environment)s_%(MeasMode)s.png' %locals())

    if PlotDoubleDerivative == 1:
        fig7 = plt.figure(7)
        ax9 = fig7.add_subplot(111)
        ax9.plot(BVolCV[20:-2], DCapf2[20:],'o')
        ax9.set_xlabel("Bias Voltage (V)")
        ax9.set_ylabel("d$^2$C/dV$^2$ (arbitrary)")
        ax9.set_title("Double Derivative, %(freq2)s(Hz)" %locals())

    if PlotDopingProfile == 1:
        fig8 = plt.figure(8)
        ax10 = fig8.add_subplot(111)
        ax10.plot(Wf1[:-1], Neff1W[:-1], 'o')
        ax10.set_xlabel("Depth (micron)")
        ax10.set_ylabel("Doping Density (10$^{11}$cm$^{-3}$)")
        ax10.set_title("%(SenName)s_Doping Profile\n%(Annealing)s_%(Environment)s, %(freq1)s(Hz)" %locals())
        fig8.savefig('%(PlotSavingLocation)s%(SenName)s_CV_IV_%(Annealing)s_%(Environment)s_%(MeasMode)s_DopingProfile_%(freq1)s(Hz).png' %locals())

        fig9 = plt.figure(9)
        ax11 = fig9.add_subplot(111)
        ax11.plot(Wf2[:-1], Neff2W[:-1], 'o')
        ax11.set_xlabel("Depth (micron)")
        ax11.set_ylabel("Doping Density (10$^{11}$cm$^{-3}$)")
        ax11.set_title("%(SenName)s_Doping Profile\n%(Annealing)s_%(Environment)s, %(freq2)s(Hz)" %locals())
        fig9.savefig('%(PlotSavingLocation)s%(SenName)s_CV_IV_%(Annealing)s_%(Environment)s_%(MeasMode)s_DopingProfile_%(freq2)s(Hz).png' %locals())

    ##fig10 = plt.figure(10)
    ##ax12 = fig10.add_subplot(111)
    ##ax12.plot(BVolIV,TempIV,'k.')
    ##ax12.set_ylabel("Temperature ($^0$C)")
    ##ax12.set_xlim(BVolIV[VminIdx],BVolIV[VmaxIdx])
    ##ax12.set_ylim(min(TempIV[:-1])-0.1,max(TempIV)+0.1)
    ##ax12.set_xlabel("Voltage (V)")
    ##ax12.spines['right'].set_color('blue')
    ##
    ##if len(DewTempCV) > 0:
    ##    ax122 = fig10.add_subplot(111, sharex=ax12, frameon=False)
    ##    ax122.plot(BVolIV,DewTempIV,'b.')
    ##    ax122.yaxis.tick_right()
    ##    ax122.yaxis.set_label_position("right")
    ##    ax122.set_ylabel("Dew ($^0$C)")
    ##    ax122.set_xlim(BVolIV[VminIdx],BVolIV[VmaxIdx])
    ##    ax122.set_ylim(min(DewTempIV[:-1])-0.1,max(DewTempIV)+0.1)
    ##    ax122.set_xlabel("Voltage (V)")
    ##    ax122.spines['right'].set_color('blue')
    ##    ax122.yaxis.label.set_color('blue')
    ##    ax122.tick_params(axis='y', colors='blue')

    show()
