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

import argparse, sys, os, time, setlog, sqlite3 as sql
import io

global log
log = setlog.init(toconsole=True)

def runlist(finput,ofil,transonly=False,orgname=False):
    if type(finput) is list:
        flist=[x.strip() for x in finput if os.path.exists(x.strip())]
    elif isinstance(finput, io.IOBase):       #elif type(finput) is file:
        flist=[x.strip() for x in finput if os.path.exists(x.strip())]
    elif os.path.exists(finput):
        with open(finput,"r") as ifil:
            flist=[x.strip() for x in ifil if os.path.exists(x.strip())]
    elif "," in finput:
        flist=[x for x in finput.split(",") if os.path.exists(x)]
    else:
        log.error("No File list found")
        return False
    numrecs=len(flist)
    if numrecs:
        conn = sql.connect(ofil)
        csr = conn.cursor()
        #Make table if doesnt exist:
        try:    
            sqlcmd="CREATE TABLE Seqs (seqid INTEGER PRIMARY KEY, orgname text, gene text, description text, source text, loc_start int, loc_end int, loc_strand int, lastscan int, naseq text, aaseq text)"
            csr.execute(sqlcmd)
            log.info("Creating db...")
        except sql.OperationalError as ex:
            log.info("Table exists. Adding to database...")
            
        # ex4istingrecs=set(str(x[0])+"@"+str(x[1]) for x in csr.execute("SELECT orgname,gene FROM Seqs"))
        #read each file
        for nr,fname in enumerate(flist):
            recs=[]
            temp=None
            aminod=None
            fp,fn = os.path.split(fname)
            org,ext = os.path.splitext(fn)
            if orgname:
                org = orgname

            #skip .faa files, test type for gbk or fasta
            if ext.lower()==".faa":
                continue
            elif ext.lower()==".fna":
                #test if .faa exists for .fna pair and get seqs in dictionary
                if len(fp):
                    fp+="/"
                if os.path.isfile(fp+org+".faa"):
                    with open(fp+org+".faa","r") as faafil:
                        aminod=[]
                        reci=-1
                        for line in faafil:
                            if line[0]==">":
                                idx=line.index(" ")
                                gn = line[1:idx]
                                # if gn not in aminod:
                                #     aminod[gn]=""
                                reci+=1
                                aminod.append([gn,""])
                            else:
                                aminod[reci][1]+=line.strip()

                def addamino(xrow,i,transtab=1):
                    if not transonly and aminod and aminod[i][0] == xrow[1]:
                        xrow.append(aminod[i][1])
                    else:
                        from Bio.Seq import Seq
                        from Bio import Data
                        try:
                            aseq=Seq(xrow[-1]).translate(to_stop=True,table=transtab)
                            xrow.append(str(aseq))
                        except Data.CodonTable.TranslationError as e:
                            log.error("Unable to get amino acid translation for reccord: %s,%s"%e)
                    return xrow
                #Read each file and collect sequences
                with open(fname,"r") as ifil:
                    reci=-1
                    transtab=1
                    for line in ifil:
                        if line[0]==">":
                            if temp: #add last temp
                                recs.append(addamino(temp,reci))
                                transtab=1
                            reci+=1
                            line=line.replace('"','').replace("'","")
                            idx=line.index(" ")
                            gn = line[1:idx]
                            source = line[1:idx].split("|")[-1]
                            ds = line[idx+1:].strip()
                            if "_PLASMID_" in ds:
                                source+="_PLASMID"
                            if "transl_table=" in ds:
                                try:
                                    transtab=int(ds.split("transl_table=")[-1])
                                except ValueError:
                                    log.debug("found translation table but is invalid")
                            if "|loc|" in ds:
                                idx = ds.index("|loc|")+5
                                locs = [int(x) for x in ds[idx:].split()[:3]]
                            else:
                                locs = (0,0,0)
                            temp = [org,gn,ds,source,locs[0],locs[1],locs[2],int(time.time()),""]
                        elif temp:
                            temp[-1]+=line.strip().upper()
                    if temp:
                        recs.append(addamino(temp,reci))

            #Insert into database
            if len(recs)>1:
                csr.executemany("INSERT INTO Seqs VALUES (NULL,?,?,?,?,?,?,?,?,?,?)",recs)
            elif len(recs):
                csr.execute("INSERT INTO Seqs VALUES (NULL,?,?,?,?,?,?,?,?,?,?)",recs[0])
            log.info("added %d of %d files (%s  - %d reccords)" % (nr+1,numrecs,fname,len(recs)))
        try:
            log.info("Removing duplicates...")
            csr.execute("DELETE FROM Seqs WHERE seqid NOT IN (SELECT min(t.seqid) FROM Seqs t GROUP BY orgname,gene,description,loc_start,loc_end,loc_strand)")
        except sql.OperationalError as ex:
            log.error("failed sql operation: %s"%ex)
        conn.commit()
        conn.close()
        return True
    else:
        log.error("No reccords found")
        return False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Combine fasta files into single sql db (Sequences from .fna extention use translation for protein seqs if .faa is not present)""")
    parser.add_argument("input", nargs="?", help="Input filelist", default=sys.stdin)
    parser.add_argument("-o", "--out", help="Output sqlite file (default: seqsql.db)", default="seqsql.db")
    parser.add_argument("-org", "--orgname", help="Organism name (default: filename)", default="")
    parser.add_argument("-t", "--trans", help="Only store translation of DNA for protein seqs (default: False)", action='store_true')
    args = parser.parse_args()
    runlist(args.input,args.out,args.trans,args.orgname)