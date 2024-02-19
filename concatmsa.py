#!/usr/bin/python
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

import glob, argparse

def runall(finput,fout,frp,fsplit,aam):
    if type(finput) is list:
        flist=finput
    else:
        flist=glob.glob(finput)
    if flist:
        aseqs={}
        plens=[0]
        ptype=[]

        for fname in flist:
            with open(fname,'r') as fil:
                l=[]
                t=0
                for line in fil:
                    if line.startswith(">"):
                        if t:
                            l.append(t)
                            t=0
                        if fsplit:
                            org=line.strip().split(fsplit)[0][1:]
                        else:
                            org=line.strip()[1:]
                        if org not in aseqs.keys():
                            aseqs[org]="-"*plens[-1]
                    else:
                        aseqs[org]+=line.strip().upper()
                        t+=len(line.strip())
                #FillGaps
                maxL=max([len(x) for x in aseqs.values()])
                for x in aseqs:
                    if len(aseqs[x]) < maxL:
                        aseqs[x]+="-"*(maxL-len(aseqs[x]))
                #crude check for type
                ptype.append(all((c=="A" or c=="T" or c=="G" or c=="C" or c=="-") for c in "".join([r[plens[-1]:] for r in aseqs.values()])))
                #Check all same length
                #l=[len(x) for x in aseqs.values()]
                if max(l)==min(l):
                    plens.append(len(aseqs.values()[0]))
                else:
                    print ("Sequences in alignment file: "+fname+" are not same length. exiting...")
                    return False
        #Write out supermatrix
        with open(fout,"w") as fout:
            for org,seq in aseqs.items():
                fout.write(">%s\n%s\n"%(org,seq))
        with open(frp,"w") as fout:
            for i,isdna in enumerate(ptype):
                part="p%s=%s-%s\n"%(i+1,plens[i]+1,plens[i+1])
                if isdna:
                    fout.write("DNA, "+part)
                else:
                    fout.write(aam.strip().upper()+", "+part)
    return (fout,frp)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Concat all Fasta format MSA files into supermatrix. Titles must match for concatenation to take place and only single copy genes allowed.")
    parser.add_argument("input",nargs="?",help="List expression of files for conversion ex: somedirectory/*.fa (default: currentdir/*.faa)",default="*.faa")
    parser.add_argument("-o","--out",help="Store results in OUT (default: supermatrix.fa)",default="supermatrix.fa")
    parser.add_argument("-p","--part",help="Store raxml partition file (default: raxmlpart.txt)",default="raxmlpart.txt")
    parser.add_argument("-s","--split",help="Character to split title ex: '|' for >Nocardia|1123 will be read as >Nocardia (default: None)",default=None)
    parser.add_argument("-aa","--aamatrix",help="AA subsitution maxtrix for protein data to use in raxml partion file (default: WAG)",default="WAG")
    args = parser.parse_args()
    runall(args.input,args.out,args.part,args.split,args.aamatrix)