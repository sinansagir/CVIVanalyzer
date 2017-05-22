#!/usr/bin/python

import os
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

CAFile = "CA_Extracted_Data.txt"
CDFile = "CD_Extracted_Data.txt"
Epi50File = "6336_Extracted_Data.txt"
MCzFile = "8556_Extracted_Data_Corrected.txt"
Epi150File = "261636_Extracted_Data.txt"
Sensors = ["FZ300um","DOFZ300um","Epi50um","MCZ300um","Epi150um"]
StrFluences = ['1.8x10$^{13}$ MeV n$_{eq}$/cm$^{2}$','4.3x10$^{13}$ MeV n$_{eq}$/cm$^{2}$','1.3x10$^{14}$ MeV n$_{eq}$/cm$^{2}$','4.6x10$^{14}$ MeV n$_{eq}$/cm$^{2}$']

def FileRead(File):
    f = open(File, 'r')
    lines = f.readlines()
    f.close()

    AnnTime = []
    DPVols = []
    DPVolErrs = []
    DPCurs = []
    DPCurErrs = []
    
    count = 0
    for line in lines:
        count +=1
        if line.startswith("Sample:"):
            data = line.strip().split()
            SensorType = data[1]
        if line.startswith("Temperature(C):"):
            data = line.strip().split()
            Temp = data[1]
        if line.startswith("Fluences"):
            data = line.strip().split()
            Fluences = data[1:]
            Fluences = map(float, Fluences)
        if line.startswith("Neff0"):
            data = line.strip().split()
            Neff0s = data[1:]
            Neff0s = map(float, Neff0s)
        if line.startswith("MEASUREMENT : DEPLETION VOLTAGE"):
            L1 = count
        if line.startswith("Depletion Voltage Errors:"):
            L2 = count
        if line.startswith("MEASUREMENT : DEPLETION CURRENT"):
            L3 = count
        if line.startswith("Depletion Current Errors:"):
            L4 = count
        if line.startswith("End"):
            L5 = count

    count = 0
    for line in lines:
        count +=1
        data = line.strip().split()
        if count >= L1+5 and count <= L2-1:
            AnnTime.append(float(data[0]))
            datatemp = map(float, data[1:])
            DPVols.append(datatemp)
        if count >= L2+1 and count <= L3-3:
            datatemp = map(float, data[1:])
            numcoulums = len(datatemp)
            DPVolErrs.append(datatemp)
        if count >= L3+5 and count <= L4-1:
            datatemp = map(float, data[1:])
            DPCurs.append(datatemp)
        if count >= L4+1 and count <= L5-3:
            datatemp = map(float, data[1:])
            DPCurErrs.append(datatemp)
    DPVolsArray = np.array(DPVols) # converting list to numpy array
    DPVolErrsArray = np.array(DPVolErrs)
    DPCursArray = np.array(DPCurs)
    DPCurErrsArray = np.array(DPCurErrs)
    return Fluences, Neff0s, AnnTime, DPVolsArray, DPVolErrsArray, DPCursArray, DPCurErrsArray, numcoulums

MCzFluences, MCzNeff0s, MCzAnnTime, MCzDPVolsArray, MCzDPVolErrsArray, MCzDPCursArray, MCzDPCurErrsArray, MCznumcoulums = FileRead(MCzFile)
CAFluences, CANeff0s, CAAnnTime, CADPVolsArray, CADPVolErrsArray, CADPCursArray, CADPCurErrsArray, CAnumcoulums = FileRead(CAFile)
CDFluences, CDNeff0s, CDAnnTime, CDDPVolsArray, CDDPVolErrsArray, CDDPCursArray, CDDPCurErrsArray, CDnumcoulums = FileRead(CDFile)
Epi50Fluences, Epi50Neff0s, Epi50AnnTime, Epi50DPVolsArray, Epi50DPVolErrsArray, Epi50DPCursArray, Epi50DPCurErrsArray, Epi50numcoulums = FileRead(Epi50File)
Epi150Fluences, Epi150Neff0s, Epi150AnnTime, Epi150DPVolsArray, Epi150DPVolErrsArray, Epi150DPCursArray, Epi150DPCurErrsArray, Epi150numcoulums = FileRead(Epi150File)

epsi0 = 8.85418 * pow(10,-12)
epsi = 11.9 # for Silicon
qe = 1.60218 * 1e-19 # C

CADPVolF1 = CADPVolsArray[:,CAnumcoulums/2].tolist()
CADPVolF2 = CADPVolsArray[:,CAnumcoulums/2+1].tolist()
CADPVolF3 = CADPVolsArray[:,CAnumcoulums/2+2].tolist()
CADPVol8min = CADPVolsArray[3,CAnumcoulums/2:CAnumcoulums/2+3].tolist()
CADPVolF1Err = CADPVolErrsArray[:,CAnumcoulums/2].tolist()
CADPVolF2Err = CADPVolErrsArray[:,CAnumcoulums/2+1].tolist()
CADPVolF3Err = CADPVolErrsArray[:,CAnumcoulums/2+2].tolist()
CADPVol8minErr = CADPVolErrsArray[3,CAnumcoulums/2:CAnumcoulums/2+3].tolist()

dCA = 300 * 1e-6 # m
CANeffF1 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF1]
CANeffF2 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF2]
CANeffF3 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF3]
CANeffF1Err = [0.045 * item for item in CANeffF1] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF1Err]
CANeffF2Err = [0.045 * item for item in CANeffF2] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF2Err]
CANeffF3Err = [0.045 * item for item in CANeffF3] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCA**2) for x in CADPVolF3Err]
CANeffF1Err[-4] = 2.0 * epsi * epsi0 * CADPVolF1Err[-4] * 1e-17 / (qe * dCA**2)
CANeffF1Err[-3] = 2.0 * epsi * epsi0 * CADPVolF1Err[-3] * 1e-17 / (qe * dCA**2)
CANeffF1Err[-2] = 2.0 * epsi * epsi0 * CADPVolF1Err[-2] * 1e-17 / (qe * dCA**2)
CANeffF2Err[-4] = 2.0 * epsi * epsi0 * CADPVolF2Err[-4] * 1e-17 / (qe * dCA**2)
CANeffF2Err[-3] = 2.0 * epsi * epsi0 * CADPVolF2Err[-3] * 1e-17 / (qe * dCA**2)
CANeffF2Err[-2] = 2.0 * epsi * epsi0 * CADPVolF2Err[-2] * 1e-17 / (qe * dCA**2)
CANeffF3Err[-4] = 2.0 * epsi * epsi0 * CADPVolF3Err[-4] * 1e-17 / (qe * dCA**2)
CANeffF3Err[-3] = 2.0 * epsi * epsi0 * CADPVolF3Err[-3] * 1e-17 / (qe * dCA**2)
CANeffF3Err[-2] = 2.0 * epsi * epsi0 * CADPVolF3Err[-2] * 1e-17 / (qe * dCA**2)

CADeltaNeffF1 = []
CADeltaNeffF2 = []
CADeltaNeffF3 = []
for item in CANeffF1:
    if item in CANeffF1[:-7]:
        CADeltaNeffF1.append(CANeff0s[0] - item)
    else:
        CADeltaNeffF1.append(CANeff0s[0] + item)
for item in CANeffF2:
    CADeltaNeffF2.append(CANeff0s[1] + item)
for item in CANeffF3:
    CADeltaNeffF3.append(CANeff0s[2] + item)

CDDPVolF1 = CDDPVolsArray[:,CDnumcoulums/2].tolist()
CDDPVolF2 = CDDPVolsArray[:,CDnumcoulums/2+1].tolist()
CDDPVolF3 = CDDPVolsArray[:,CDnumcoulums/2+2].tolist()
CDDPVol8min = CDDPVolsArray[3,CDnumcoulums/2:CDnumcoulums/2+3].tolist()
CDDPVolF1Err = CDDPVolErrsArray[:,CDnumcoulums/2].tolist()
CDDPVolF2Err = CDDPVolErrsArray[:,CDnumcoulums/2+1].tolist()
CDDPVolF3Err = CDDPVolErrsArray[:,CDnumcoulums/2+2].tolist()
CDDPVol8minErr = CDDPVolErrsArray[3,CDnumcoulums/2:CDnumcoulums/2+3].tolist()

dCD = 300 * 1e-6 # m
CDNeffF1 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF1]
CDNeffF2 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF2]
CDNeffF3 = [2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF3]
CDNeffF1Err = [0.045 * item for item in CDNeffF1] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF1Err]
CDNeffF2Err = [0.045 * item for item in CDNeffF2] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF2Err]
CDNeffF3Err = [0.045 * item for item in CDNeffF3] #[2.0 * epsi * epsi0 * x * 1e-17 / (qe * dCD**2) for x in CDDPVolF3Err]
CDNeffF1Err[-4] = 2.0 * epsi * epsi0 * CDDPVolF1Err[-4] * 1e-17 / (qe * dCD**2)
CDNeffF1Err[-3] = 2.0 * epsi * epsi0 * CDDPVolF1Err[-3] * 1e-17 / (qe * dCD**2)
CDNeffF1Err[-2] = 2.0 * epsi * epsi0 * CDDPVolF1Err[-2] * 1e-17 / (qe * dCD**2)
CDNeffF2Err[-4] = 2.0 * epsi * epsi0 * CDDPVolF2Err[-4] * 1e-17 / (qe * dCD**2)
CDNeffF2Err[-3] = 2.0 * epsi * epsi0 * CDDPVolF2Err[-3] * 1e-17 / (qe * dCD**2)
CDNeffF2Err[-2] = 2.0 * epsi * epsi0 * CDDPVolF2Err[-2] * 1e-17 / (qe * dCD**2)
CDNeffF3Err[-4] = 2.0 * epsi * epsi0 * CDDPVolF3Err[-4] * 1e-17 / (qe * dCD**2)
CDNeffF3Err[-3] = 2.0 * epsi * epsi0 * CDDPVolF3Err[-3] * 1e-17 / (qe * dCD**2)
CDNeffF3Err[-2] = 2.0 * epsi * epsi0 * CDDPVolF3Err[-2] * 1e-17 / (qe * dCD**2)

CDDeltaNeffF1 = []
CDDeltaNeffF2 = []
CDDeltaNeffF3 = []
for item in CDNeffF1:
    if item in CDNeffF1[:-7]:
        CDDeltaNeffF1.append(CDNeff0s[0] - item)
    else:
        CDDeltaNeffF1.append(CDNeff0s[0] + item)
for item in CDNeffF2:
    CDDeltaNeffF2.append(CDNeff0s[1] + item)
for item in CDNeffF3:
    CDDeltaNeffF3.append(CDNeff0s[2] + item)

Fluence1 = StrFluences[0]
Fluence2 = StrFluences[1]
Fluence3 = StrFluences[2]
Fluence4 = StrFluences[3]

def general_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
    
    popt, pcov = curve_fit(f, xdata, ydata, p0, sigma, maxfev=100000)
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

def Hamburg(x, ga, NC, gy, ta, ty):
    fluence = CurrentFluence
    return ga*exp(-x/ta)*fluence + gy*(1.-1./(1.+x/ty))*fluence + NC

def Hamburg2(x, ga, NC, gy, ta, ty, gy2, ty2):
    fluence = CurrentFluence
    return ga*exp(-x/ta)*fluence + gy*(1.-1./(1.+x/ty))*fluence + NC + gy2*log(fluence*x/ty2)

CurrentFluence = CAFluences[0]*100.0
p0CAF1 = [ 0.001, 11., 0.04, 15., 500.]
poptCAF1, puncCAF1, rcCAF1, dCAF1 = general_fit(Hamburg, np.array(CAAnnTime[:-1]), np.array(CADeltaNeffF1[:-1]), p0CAF1, np.array(CANeffF1Err[:-1]))
CAAnnTimeFit = np.linspace(min(np.array(CAAnnTime[:-1])),max(np.array(CAAnnTime[:-1])),40000)
CADeltaNeffF1Fit = Hamburg(CAAnnTimeFit, poptCAF1[0], poptCAF1[1], poptCAF1[2], poptCAF1[3], poptCAF1[4])
print "************************FZ-1.8E13MeV/cm-2**************************"
print "ga = ", poptCAF1[0], "+\-", puncCAF1[0]
print "NC = ", poptCAF1[1], "+\-", puncCAF1[1]
print "gy = ", poptCAF1[2], "+\-", puncCAF1[2]
print "ta = ", poptCAF1[3], "+\-", puncCAF1[3]
print "ty = ", poptCAF1[4], "+\-", puncCAF1[4]
print "chi2 = ", rcCAF1
print "dof = ", dCAF1
print "chi2_red = ", rcCAF1/dCAF1
print "******************************************************************"
CurrentFluence = CAFluences[1]*100.0
p0CAF2 = [ 0.1, 10., 0.06, 2., 800.]
poptCAF2, puncCAF2, rcCAF2, dCAF2 = general_fit(Hamburg, np.array(CAAnnTime[:-1]), np.array(CADeltaNeffF2[:-1]), p0CAF2, np.array(CANeffF2Err[:-1]))
CAAnnTimeFit = np.linspace(min(np.array(CAAnnTime[:-1])),max(np.array(CAAnnTime[:-1])),40000)
CADeltaNeffF2Fit = Hamburg(CAAnnTimeFit, poptCAF2[0], poptCAF2[1], poptCAF2[2], poptCAF2[3], poptCAF2[4])
print "************************FZ-4.3E13MeV/cm-2**************************"
print "ga = ", poptCAF2[0], "+\-", puncCAF2[0]
print "NC = ", poptCAF2[1], "+\-", puncCAF2[1]
print "gy = ", poptCAF2[2], "+\-", puncCAF2[2]
print "ta = ", poptCAF2[3], "+\-", puncCAF2[3]
print "ty = ", poptCAF2[4], "+\-", puncCAF2[4]
print "chi2 = ", rcCAF2
print "dof = ", dCAF2
print "chi2_red = ", rcCAF2/dCAF2
print "******************************************************************"
CurrentFluence = CAFluences[2]*100.0
p0CAF3 = [ 0.01, 11., 0.04, 15., 100.]
poptCAF3, puncCAF3, rcCAF3, dCAF3 = general_fit(Hamburg, np.array(CAAnnTime[:-1]), np.array(CADeltaNeffF3[:-1]), p0CAF3, np.array(CANeffF3Err[:-1]))
CAAnnTimeFit = np.linspace(min(np.array(CAAnnTime[:-1])),max(np.array(CAAnnTime[:-1])),40000)
CADeltaNeffF3Fit = Hamburg(CAAnnTimeFit, poptCAF3[0], poptCAF3[1], poptCAF3[2], poptCAF3[3], poptCAF3[4])
print "************************FZ-1.3E14MeV/cm-2**************************"
print "ga = ", poptCAF3[0], "+\-", puncCAF3[0]
print "NC = ", poptCAF3[1], "+\-", puncCAF3[1]
print "gy = ", poptCAF3[2], "+\-", puncCAF3[2]
print "ta = ", poptCAF3[3], "+\-", puncCAF3[3]
print "ty = ", poptCAF3[4], "+\-", puncCAF3[4]
print "chi2 = ", rcCAF3
print "dof = ", dCAF3
print "chi2_red = ", rcCAF3/dCAF3
print "******************************************************************"
CurrentFluence = CDFluences[0]*100.0
p0CDF1 = [ 0.01, 11., 0.04, 3., 500.]
poptCDF1, puncCDF1, rcCDF1, dCDF1 = general_fit(Hamburg, np.array(CDAnnTime[:-1]), np.array(CDDeltaNeffF1[:-1]), p0CDF1, np.array(CDNeffF1Err[:-1]))
CDAnnTimeFit = np.linspace(min(np.array(CDAnnTime[:-1])),max(np.array(CDAnnTime[:-1])),40000)
CDDeltaNeffF1Fit = Hamburg(CDAnnTimeFit, poptCDF1[0], poptCDF1[1], poptCDF1[2], poptCDF1[3], poptCDF1[4])
print "***********************DOFZ-1.8E13MeV/cm-2*************************"
print "ga = ", poptCDF1[0], "+\-", puncCDF1[0]
print "NC = ", poptCDF1[1], "+\-", puncCDF1[1]
print "gy = ", poptCDF1[2], "+\-", puncCDF1[2]
print "ta = ", poptCDF1[3], "+\-", puncCDF1[3]
print "ty = ", poptCDF1[4], "+\-", puncCDF1[4]
print "chi2 = ", rcCDF1
print "dof = ", dCDF1
print "chi2_red = ", rcCDF1/dCDF1
print "******************************************************************"
CurrentFluence = CDFluences[1]*100.0
p0CDF2 = [ 0.01, 8., 0.04, 15., 100.]
poptCDF2, puncCDF2, rcCDF2, dCDF2 = general_fit(Hamburg, np.array(CDAnnTime[:-1]), np.array(CDDeltaNeffF2[:-1]), p0CDF2, np.array(CDNeffF2Err[:-1]))
CDAnnTimeFit = np.linspace(min(np.array(CDAnnTime[:-1])),max(np.array(CDAnnTime[:-1])),40000)
CDDeltaNeffF2Fit = Hamburg(CDAnnTimeFit, poptCDF2[0], poptCDF2[1], poptCDF2[2], poptCDF2[3], poptCDF2[4])
print "***********************DOFZ-4.3E13MeV/cm-2*************************"
print "ga = ", poptCDF2[0], "+\-", puncCDF2[0]
print "NC = ", poptCDF2[1], "+\-", puncCDF2[1]
print "gy = ", poptCDF2[2], "+\-", puncCDF2[2]
print "ta = ", poptCDF2[3], "+\-", puncCDF2[3]
print "ty = ", poptCDF2[4], "+\-", puncCDF2[4]
print "chi2 = ", rcCDF2
print "dof = ", dCDF2
print "chi2_red = ", rcCDF2/dCDF2
print "******************************************************************"
CurrentFluence = CDFluences[2]*100.0
p0CDF3 = [ 0.01, 11., 0.04, 15., 400.]
poptCDF3, puncCDF3, rcCDF3, dCDF3 = general_fit(Hamburg, np.array(CDAnnTime[:-1]), np.array(CDDeltaNeffF3[:-1]), p0CDF3, np.array(CDNeffF3Err[:-1]))
CDAnnTimeFit = np.linspace(min(np.array(CDAnnTime[:-1])),max(np.array(CDAnnTime[:-1])),40000)
CDDeltaNeffF3Fit = Hamburg(CDAnnTimeFit, poptCDF3[0], poptCDF3[1], poptCDF3[2], poptCDF3[3], poptCDF3[4])
print "***********************DOFZ-1.3E14MeV/cm-2*************************"
print "ga = ", poptCDF3[0], "+\-", puncCDF3[0]
print "NC = ", poptCDF3[1], "+\-", puncCDF3[1]
print "gy = ", poptCDF3[2], "+\-", puncCDF3[2]
print "ta = ", poptCDF3[3], "+\-", puncCDF3[3]
print "ty = ", poptCDF3[4], "+\-", puncCDF3[4]
print "chi2 = ", rcCDF3
print "dof = ", dCDF3
print "chi2_red = ", rcCDF3/dCDF3
print "******************************************************************"

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
#ax1.set_title("Depletion Voltage vs. Annealing Time (FZ300um, 10kHz, 80$^0$C)")
ax1.errorbar(CAAnnTime[:-1],CADPVolF1[:-1], yerr=CADPVolF1Err[:-1], fmt='ko--', label='%(Fluence1)s' %locals())
ax1.errorbar(CAAnnTime[:-1],CADPVolF2[:-1], yerr=CADPVolF2Err[:-1], fmt='go--', label='%(Fluence2)s' %locals())
ax1.errorbar(CAAnnTime[:-1],CADPVolF3[:-1], yerr=CADPVolF3Err[:-1], fmt='bo--', label='%(Fluence3)s' %locals())
ax1.set_xscale("log", nonposx='clip')
ax1.set_xlabel("Annealing Time [min]")
ax1.set_ylabel("$V_{dep}$ [V]")
ax1.legend(loc=2,prop={'size':12})
ax1.set_xlim(1,10000)
ax12 = ax1.twinx()
ax12.set_ylabel("$N_{eff}$ [$x10^{11}\,cm^{-3}$]")
ax12.set_ylim(ax1.get_ylim())
ax12.set_yticks(ax1.get_yticks())
ax12.set_yticklabels([int(2.0 * epsi * epsi0 * item * 1e-17 / (qe * dCA**2)) for item in ax1.get_yticks()])

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
#ax2.set_title("Depletion Voltage vs. Annealing Time (DOFZ300um, 10kHz, 80$^0$C)")
ax2.errorbar(CDAnnTime[:-1],CDDPVolF1[:-1], yerr=CDDPVolF1Err[:-1], fmt='ko--', label='%(Fluence1)s' %locals())
ax2.errorbar(CDAnnTime[:-1],CDDPVolF2[:-1], yerr=CDDPVolF2Err[:-1], fmt='go--', label='%(Fluence2)s' %locals())
ax2.errorbar(CDAnnTime[:-1],CDDPVolF3[:-1], yerr=CDDPVolF3Err[:-1], fmt='bo--', label='%(Fluence3)s' %locals())
ax2.set_xscale("log", nonposx='clip')
ax2.set_xlabel("Annealing Time [min]")
ax2.set_ylabel("$V_{dep}$ [V]")
ax2.legend(loc=2,prop={'size':12})
ax2.set_xlim(1,10000)
ax2.set_ylim(ax2.get_ylim())
ax22 = ax2.twinx()
ax22.set_ylabel("$N_{eff}$ [$x10^{11}\,cm^{-3}$]")
ax22.set_ylim(ax2.get_ylim())
ax22.set_yticks(ax2.get_yticks())
ax22.set_yticklabels([int(2.0 * epsi * epsi0 * item * 1e-17 / (qe * dCD**2)) for item in ax2.get_yticks()])


fig7 = plt.figure(7)
ax7 = fig7.add_subplot(111)
#ax7.set_title("$\\Delta$$N_{eff}$ vs. Annealing Time (FZ300um&DOFZ300um, 10kHz, 80$^0$C)")
ax7.errorbar(CAAnnTime[:-1],CADeltaNeffF1[:-1], yerr=CANeffF1Err[:-1], fmt='ko', label='FZ-%(Fluence1)s' %locals())
ax7.semilogx(CAAnnTimeFit,CADeltaNeffF1Fit,'k')
ax7.errorbar(CAAnnTime[:-1],CADeltaNeffF2[:-1], yerr=CANeffF2Err[:-1], fmt='go', label='FZ-%(Fluence2)s' %locals())
ax7.semilogx(CAAnnTimeFit,CADeltaNeffF2Fit,'g')
ax7.errorbar(CAAnnTime[:-1],CADeltaNeffF3[:-1], yerr=CANeffF3Err[:-1], fmt='bo', label='FZ-%(Fluence3)s' %locals())
ax7.semilogx(CAAnnTimeFit,CADeltaNeffF3Fit,'b')
ax7.errorbar(CDAnnTime[:-1],CDDeltaNeffF1[:-1], yerr=CDNeffF1Err[:-1], fmt='ro', label='DOFZ-%(Fluence1)s' %locals())
ax7.semilogx(CDAnnTimeFit,CDDeltaNeffF1Fit,'r')
ax7.errorbar(CDAnnTime[:-1],CDDeltaNeffF2[:-1], yerr=CDNeffF2Err[:-1], fmt='co', label='DOFZ-%(Fluence2)s' %locals())
ax7.semilogx(CDAnnTimeFit,CDDeltaNeffF2Fit,'c')
ax7.errorbar(CDAnnTime[:-1],CDDeltaNeffF3[:-1], yerr=CDNeffF3Err[:-1], fmt='yo', label='DOFZ-%(Fluence3)s' %locals())
ax7.semilogx(CDAnnTimeFit,CDDeltaNeffF3Fit,'y')
ax7.set_xlim(1,10000)
ax7.set_xscale("log", nonposx='clip')
ax7.set_xlabel("Annealing Time [min]")
ax7.set_ylabel("$\\Delta$$N_{eff}$ [$x10^{11}\,cm^{-3}$]")
ax7.legend(loc=2,prop={'size':12})

show()
