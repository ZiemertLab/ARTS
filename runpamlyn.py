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

import argparse, tempfile, os, glob, subprocess, sys, multiprocessing as mp

def convertandrun(finput,outdir):
    finput = os.path.realpath(finput)
    fpath, fname = os.path.split(finput)
    outfname = os.path.join(os.path.realpath(outdir),fname)
    cmd = ['seqret', '-sequence', str(finput), '-outseq', '%s.phy'%outfname, '-sformat1', 'fasta', '-osf', 'phylipnon']
    #write ctrl file
    with open(os.path.join(outdir,fname+".ctrl"),"w") as fil:
        fil.write("seqfile = %s.phy\noutfile = %s.kaks\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n"%(outfname,outfname))
    with open(os.devnull,"w") as devnull:
        subprocess.call(cmd,stdout=devnull,stderr=sys.stderr)
        subprocess.call(["yn00",outfname+".ctrl"],stdout=devnull,stderr=devnull)
    return outfname+".kaks"

def runall(finput, outdir=None, cpu=mp.cpu_count()):
    if type(finput)==list:
        flist=finput
    else:
        flist=glob.glob(finput)

    #Output result
    if not outdir:
        outdir = tempfile.mkdtemp()
    elif not os.path.isdir(outdir):
        os.makedirs(os.path.relpath(outdir))

    if cpu > 1:
        pool = mp.Pool()
        for filename in flist:
            pool.apply_async(convertandrun, args=(filename, outdir))
        pool.close()
        pool.join()
    print(outdir)
    return outdir

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert Fasta codnon-alignment to Phylip and run yn00 for Ka Ks calculation in PAML")
    parser.add_argument("input",nargs="?",help="List expression of files for conversion ex: somedirectory/*.fna (default: currentdir/*.fna)",default="*.fna")
    parser.add_argument("-o","--out",help="Store results in OUT (default: temporary folder)",default=None)
    parser.add_argument("-c","--cpu",help="Number of parallel cpus to use (default = all)",default=mp.cpu_count())
    args = parser.parse_args()
    runall(args.input,args.out,args.cpu)