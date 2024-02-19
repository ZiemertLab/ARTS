#!/usr/bin/env python
# Copyright (C) 2015,2016 Mohammad Alanjary
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
#
# This file is part of ARTS
# ARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
#
# License: You should have received a copy of the GNU General Public License v3 with ARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import argparse, glob, tempfile, numpy as np

def runall(finput, fout=None, lt=1.0, gt=1.0, method="NG86"):
    Allrecs={}
    if type(finput)==list:
        flist=finput
    else:
        flist=glob.glob(finput)
    for fname in flist:
        recs={}
        recs['NG86']=[] #(A) Nei-Gojobori (1986) method
        recs['YN2K']=[] #(B) Yang & Nielsen (2000) method
        # recs['LWL85']=[] #(C) LWL85, LPB93 & LWLm methods
        recs['LWL85M']=[] #(C) LWL85, LPB93 & LWLm methods
        recs['LPB93']=[] #(C) LWL85, LPB93 & LWLm methods

        #yang_kappa=None
        mthd=""
        seqlabels=[]
        L=[0,0]
        with open(fname,"r") as fil:
            for line in fil:
                try:
                    if line.startswith("(A)"):
                        mthd="A"
                    elif line.startswith("(B)"):
                        mthd="B"
                    elif line.startswith("(C)"):
                        mthd="C"

                    #Nei-Gojobori
                    if mthd=="A" and "      " in line:
                        x = line.replace("("," ").replace(")"," ").split()
                        seqlabels.append(x[0])
                        if len(x)>1:
                            w = x[1::3]
                            dN = x[2::3]
                            dS = x[3::3]
                            for i in range(len(w)):  #for i in xrange(len(w)):
                                # Store record: [seq#1, seq#2, omega, dN, dS]
                                recs["NG86"].append([len(seqlabels),i,float(w[i]),float(dN[i]),float(dS[i])])

                    #Yang and Nielsen
                    if mthd=="B" and "+-" in line and not "seq." in line:
                        x=line.split()
                        #Store record: [seq#1, seq#2, omega, dN, dS, dN std error, dS std error,]
                        omega=float(x[6])
                        if float(x[10])==99: #Flag extreme dS error
                            omega=-1.0
                        recs["YN2K"].append([int(x[0]),int(x[1]),float(x[6]),float(x[7]),omega,float(x[9]),float(x[-1])])
                        #if not yang_kappa:
                            #yang_kappa=float(x[5])

                    #LWL & LPB
                    if mthd=="C":
                        if "vs." in line:
                            L=line.split("vs.")
                            L[0]=int(L[0][0:L[0].index("(")])
                            L[1]=int(L[1][0:L[1].index("(")])
                        #if line.startswith("LWL85") or line.startswith("LWL85M") or line.startswith("LPB93"):
                        if line.startswith("LWL85m") or line.startswith("LPB93"):
                            x = line.split()
                            recs[x[0][0:x[0].index(":")].upper()].append(L+[float(v) for v in x[3:10:3][::-1]])
                except Exception as e:
                    print("error reading values %s"%e)
        #Merge and summarize
        def isval(v):
            return (v>=0 and v!=(v+1))
        klist = recs.keys()
        for k in klist:
            try:
                omegaAV = np.median([x[2] for x in recs[k] if isval(x[2])])
                omegaSD = np.std([x[2] for x in recs[k] if isval(x[2])])
                dNAV = np.median([x[3] for x in recs[k] if isval(x[3])])
                dNSD = np.median([x[3] for x in recs[k] if isval(x[3])])
                dSAV = np.median([x[4] for x in recs[k] if isval(x[4])])
                dSSD = np.median([x[4] for x in recs[k] if isval(x[4])])
                recs[k+"_AVG"]=[omegaAV,dNAV,dSAV]
                recs[k+"_STD"]=[omegaSD,dNSD,dSSD]
            except Exception as e:
                print("Error calculating median %s"%e)
        Allrecs[fname[0:fname.index(".")]]=recs

    #Get gene list under threshold and calculate method agreement:
    LTlist = []
    GTlist = []
    for k,recs in Allrecs.items():
        x = recs[method+"_AVG"][0]
        if np.median(x) <= lt:
           LTlist.append(k)
        elif np.median(x) > gt:
           GTlist.append(k)

    #Output result
    if fout:
        with open("median_"+fout,"w") as afil:
            with open("std_"+fout,"w") as sfil:
                #head="#GENE\t#Nei-Gojobori_omega_dN_dS\t#Yang-Nielsen_omega_dN_dS\t#Li-Wu-Luo85_omega_dN_dS\t#Li-Wu-Luo93_omega_dN_dS\t#Pamilo-Bianchi_omega_dN_dS\t#Method_median_RSD"
                head="#GENE\t#Nei-Gojobori_omega_dN_dS\t#Yang-Nielsen_omega_dN_dS\t#Li-Wu-Luo93_omega_dN_dS\t#Pamilo-Bianchi_omega_dN_dS\t#Method_median_RSD"
                afil.write(head+"\n")
                sfil.write(head+"\n")
                for k,recs in Allrecs.items():
                    mOmega = [recs[m][0] for m in recs if m.endswith("_AVG")]
                    mOAR = [np.median(mOmega),100*np.std(mOmega)/np.median(mOmega)] # median and RSD omega of all methods
                    afil.write("%s\t%.4f %.4f %.4f\t%.4f %.4f %.4f\t%.4f %.4f %.4f\t%.4f %.4f %.4f\t%.4f %.2f%%\n"%tuple([k]+list(recs.get("NG86_AVG",[-1,-1,-1]))+list(recs.get("YN2K_AVG",[-1,-1,-1]))+list(recs.get("LWL85M_AVG",[-1,-1,-1]))+list(recs.get("LPB93_AVG",[-1,-1,-1]))+mOAR))
                    sfil.write("%s\t%.4f %.4f %.4f\t%.4f %.4f %.4f\t%.4f %.4f %.4f\t%.4f %.4f %.4f\n"%tuple([k]+list(recs.get("NG86_STD",[-1,-1,-1]))+list(recs.get("YN2K_STD",[-1,-1,-1]))+list(recs.get("LWL85M_STD",[-1,-1,-1]))+list(recs.get("LPB93_STD",[-1,-1,-1]))))
                afil.write("##Filtered list using "+method+"\n")
                afil.write("#dN/dS>"+str(gt)+":\t"+",".join(GTlist)+"\n")
                afil.write("#dN/dS<="+str(lt)+":\t"+",".join(LTlist)+"\n")
    else:
        return LTlist,GTlist,Allrecs
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PAML yn00 output and return dictionary or summary table")
    parser.add_argument("input",nargs="?",help="List expression of files for conversion ex: somedirectory/*.txt (default: currentdir/*.kaks)",default="*.kaks")
    parser.add_argument("-o","--out",help="Store results in OUT (default: summarykaks.txt)",default="summarykaks.txt")
    parser.add_argument("-m","--method",help="Use this method to threshold list ['NG86'=Nei-Gojobori(1986),'YN2K'=Yang-Nielsen(2000),'LWL85M'=Li-Wu(1993),'LPB93'=Pamilo-Bianchi(1993)] (default='NG86')",default="NG86")
    parser.add_argument("-lt","--lthrsh",help="Threshold value to print all genes under:  avg(omega)+std(omega) <= thrsh (default=1.0)",default=1.0)
    parser.add_argument("-gt","--gthrsh",help="Threshold value to print all genes over:  avg(omega)+std(omega) > thrsh (default=1.0)",default=1.0)
    args = parser.parse_args()
    runall(args.input,args.out,args.lthrsh,args.gthrsh,args.method)