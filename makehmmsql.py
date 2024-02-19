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
import argparse, ast, os, setlog, sqlite3 as sql
import time, numpy as np

global log
log = setlog.init(toconsole=True)

def gethmmcov(clist):
    for cl in clist:
        spots = set()
        spots2 = set()
        for x in cl[2].split('&'):
            xt=ast.literal_eval(x)
            spots|=set(range(xt[0],xt[1]))     #spots|=set(xrange(xt[0],xt[1]))
        for x in cl[3].split('&'):
            xt=ast.literal_eval(x)
            spots2|=set(range(xt[0],xt[1]))    #spots2|=set(xrange(xt[0],xt[1]))
        yield int(cl[0]),cl[1],len(spots)+1,len(spots2)+1

def run(fname,ofil,ev=1e-4,bs=0.0,rna=False,filt2=None):
    conn = sql.connect(ofil)
    csr = conn.cursor()
    tabletitle="HMMhits"
    if rna:
        tabletitle="RNAhits"
    try:
        sqlcmd="""
            CREATE TABLE """+tabletitle+""" (hmmhit text, orgname text, seqid int, hmmstart int, hmmend int, hmmlen int,
            geneStart int, geneEnd int, genelen int, evalue real, score real, bias real, iscore real, flags int, hmmcov real, genecov real)"""
        csr.execute(sqlcmd)
        log.info("Creating db...")
    except sql.OperationalError as ex:
        log.info("Table exists. Adding to database...")
    #set filter dict form file:
    if filt2 and os.path.isfile(filt2):
        temp={}
        with open(filt2,"r") as fil:
            for line in fil:
                x=line.strip().split()
                temp[x[0]]=float(x[1])
        filt2=temp
    #Add hmmhits:
    with open(fname,"r") as ifil:
        recs=[]
        for line in ifil:
            if line[0]!="#":
                x = line.split()
                org,seqid = x[0].split("|")
                if filt2 and x[3] in filt2.keys() and filt2[x[3]]>float(x[7]):
                    continue
                if float(x[6]) < ev and float(x[7]) > bs:  #global quality filter
                    recs.append([x[3],org,int(seqid),int(x[15]),int(x[16]),int(x[5]),int(x[19]),int(x[20]),int(x[2]),float(x[6]),float(x[7]),float(x[8]),float(x[13]),0,0,0])
        else:
            csr.executemany("INSERT INTO "+tabletitle+" VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",recs)
            log.info("Added %d records"%len(recs))
            del recs
    # Remove duplicate rows
    try:
        log.info("Removing duplicates...")
        csr.execute("DELETE FROM "+tabletitle+" WHERE rowid NOT IN (SELECT min(t.rowid) FROM "+tabletitle+" t GROUP BY hmmhit,seqid,hmmstart,hmmend,geneStart,geneEnd)")
        log.info("Marking best hits...")
        csr.execute("UPDATE "+tabletitle+" SET flags=0") # clear all previous
        csr.execute("UPDATE "+tabletitle+" SET flags=10 WHERE seqid||hmmhit IN (SELECT seqid||hmmhit FROM "+tabletitle+" GROUP BY seqid HAVING evalue==MIN(evalue))")
    except sql.OperationalError as ex:
        log.error("exception:"%ex)

    temptable="temp_%d%d"%(os.getpid(),time.time()) # make temporary table name 
    try:
        log.info("Calculating total coverage...")
        
        #get hmmalign coords:
        csr.execute("SELECT seqid,hmmhit,Group_concat('('||hmmstart||','||hmmend||')','&'),Group_concat('('||geneStart||','||geneEnd||')','&') FROM "+tabletitle+" GROUP BY seqid,hmmhit")
        seqcovs=gethmmcov(csr.fetchall())
        csr.execute("CREATE TEMP TABLE %s_cov (seqid int, hmmhit text, hmmcov real, genecov real)"%temptable)
        csr.executemany("INSERT INTO %s_cov VALUES (?,?,?,?)"%temptable,seqcovs)

        sqlcmd="""
            CREATE TABLE %s AS SELECT A.hmmhit,A.orgname,A.seqid,A.hmmstart,A.hmmend,A.hmmlen,A.geneStart,A.geneEnd,A.genelen,A.evalue,A.score,A.bias,A.iscore,A.flags,1.*B.hmmcov/A.hmmlen as hmmcov,1.*B.genecov/A.genelen as genecov
            FROM %s A INNER JOIN %s_cov B ON A.seqid=B.seqid AND A.hmmhit=B.hmmhit""" % (temptable,tabletitle,temptable)
        csr.execute(sqlcmd)
        csr.execute("DROP TABLE "+tabletitle)
        # csr.execute("DROP TABLE %s_cov"%temptable)
        csr.execute("ALTER TABLE "+temptable+" RENAME TO "+tabletitle)
        # print "Compacting db size..."
        # csr.execute("VACUUM")
    except sql.OperationalError as ex:
        log.error("Aborting coverage calculation, Error:%s"%ex)

    conn.commit()    
    conn.close()
            
# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Combine hmmresults into single sql db""")
    parser.add_argument("input", help="Input Hmm result file")
    parser.add_argument("out", help="Output database")
    parser.add_argument("-r", "--rna", help="Add to RNA table instead of HK gene table", action='store_true')
    parser.add_argument("-e", "--evalue", help="Remove hits > evalue default(1e-4)", type=float, default=1e-4)
    parser.add_argument("-b", "--bitscore", help="Remove hits < bitscore default(0)", type=float, default=0)
    parser.add_argument("-f2", "--filter", help="Use per-model bitscore cutoffs dictated in file (format: modelname value)", default=None)
    args = parser.parse_args()
    run(args.input,args.out,args.evalue,args.bitscore,args.rna,args.filter)