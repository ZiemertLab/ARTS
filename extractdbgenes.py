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

import os, argparse, ast, setlog, tempfile, json, sqlite3 as sql, numpy as np

global log
log = setlog.init(toconsole=True)

def writeall(outdir, db, pct=0.5, pct2=0.5, ev=0.1, bs=0, bh=True, filt2=None, topseq=False, idprfx=""):
    if type(filt2) is str and os.path.isfile(filt2):
        with open(filt2,"r") as filtfil:
            temp={}
            for line in filtfil:
                x=line.strip().split()
                temp[x[0]]=float(x[1])
            filt2=temp
    if not outdir:
        outdir = tempfile.mkdtemp(prefix="results_", dir="./")
        log.warning("NO RESULT FOLDER SPECIFIED. Results will be stored in temporary folder: " + outdir)
    elif not os.path.exists(outdir):
        os.mkdir(outdir)
    if not outdir.endswith("/"):
        outdir += "/"

    conn = sql.connect(db)
    cur = conn.cursor()
    hkcore=[str(x[0]) for x in cur.execute("SELECT DISTINCT hmmhit FROM HMMhits")]
    orgs=[str(x[0]) for x in cur.execute("SELECT DISTINCT orgname FROM HMMhits")]
    orgs = list(orgs)
    oarry = np.zeros((len(hkcore), len(orgs)))
    coredict = {"orgs":orgs,"singles":[],"core":{},"seqs":{}}
    # oarry2 = np.zeros((len(hkcore), len(orgs)))
    addwhere = " AND A.orgname IN " + str(orgs).replace("[", "(").replace("]", ")")
    if bh:
        # best hits only
        addwhere += " AND flags>0"
    if bs:
        addwhere += " AND score>%s"%bs
    # Combine tables into view and write full hits, fragment hits, genematrix
    cur.execute("CREATE TEMP VIEW allgenes AS SELECT A.*,B.* FROM HMMhits A INNER JOIN Seqs B ON A.seqid=B.seqid WHERE evalue<=%s"%ev + addwhere)
    for i, hkgn in enumerate(hkcore):
        results = cur.execute(
            "SELECT orgname,seqid,genecov,hmmcov,MAX(iscore),genelen,source,loc_start,loc_end,loc_strand,description,score,naseq,aaseq FROM allgenes WHERE hmmhit=? GROUP BY seqid",
            (hkgn,)).fetchall()
        with open(outdir + hkgn + ".faa", "w") as aafil, open(outdir + hkgn + ".fna", "w") as nafil:
            orgseqs = {} # store here before writing to check for top hit
            for x in results:
                # add to genematrix
                j = orgs.index(x[0])
                seqs = (x[-2], x[-1])
                # separate fragments and full length
                if filt2 and hkgn in filt2.keys() and float(filt2[hkgn]) > float(x[4]):
                    continue
                if x[2] >= pct or x[3] >= pct2:
                    oarry[i, j] += 1
                    seqid = idprfx+str(x[1])
                    # aafil.write(">%s|%s\n%s\n" % (x[0], seqid, seqs[1]))
                    # nafil.write(">%s|%s\n%s\n" % (x[0], seqid, seqs[0]))
                    if str(x[0]) not in orgseqs:
                        orgseqs[str(x[0])] = {}
                        orgseqs[str(x[0])]["all"] = []
                        orgseqs[str(x[0])]["max"] = [-1]
                    orgseqs[str(x[0])]["all"].append([float(x[-3]),seqid,seqs[1],seqs[0]])
                    if float(x[-3]) > orgseqs[str(x[0])]["max"][0]:
                        orgseqs[str(x[0])]["max"] = [float(x[-3]),seqid,seqs[1],seqs[0]]
                    #Add to core list (seqid, hkgene, source, location start, location end, strand, description, bitscore)
                    if hkgn not in coredict["core"].keys():
                        coredict["core"][hkgn] = {"seqs":[]}
                    if seqid not in coredict["seqs"].keys():
                        coredict["seqs"][seqid] = [seqid,hkgn,str(x[6]),x[7],x[8],x[9],x[10],x[11]]
                    coredict["core"][hkgn]["seqs"].append(seqid)
            #write out stored seqs
            for org in orgseqs.keys():
                if topseq:
                    seqrow = orgseqs[org]["max"]
                    aafil.write(">%s|%s\n%s\n" % (org, seqrow[1], seqrow[2]))
                    nafil.write(">%s|%s\n%s\n" % (org, seqrow[1], seqrow[3]))
                else:
                    seqrow = orgseqs[org]["all"]
                    for row in seqrow:
                        aafil.write(">%s|%s\n%s\n" % (org, row[1], row[2]))
                        nafil.write(">%s|%s\n%s\n" % (org, row[1], row[3]))

        log.debug("Wrote (%d of %d): %s" % (i + 1, len(hkcore), hkgn))
    with open(outdir + "genematrix.txt", "w") as gmfil:
        gmfil.write("#Gene\tMedianCount\tStDev\tSinglesRatio\tUbiquity\t" + "\t".join(orgs) + "\n")
        singleslist = []
        for i, k in enumerate(hkcore):
            gr = [x for x in oarry[i] if x > 0]
            if len(gr):
                grmedian = float(np.median(gr))
                grstd = float(np.std(gr))
                line = "%s\t%.4f\t%.4f\t%.4f\t%.4f\t" % (k, grmedian, grstd, float(list(oarry[i]).count(1)) / len(oarry[i]), float(len(gr))/len(oarry[i]))
                gmfil.write(line + "\t".join([str(x) for x in oarry[i]]) + "\n")

                coredict["core"][k]["medcount"] = grmedian
                coredict["core"][k]["stdev"] = grstd
                coredict["core"][k]["countv"] = [int(x) for x in oarry[i]]
                coredict["core"][k]["ubiquity"] = float(len(gr))/len(oarry[i])
            else:
                os.remove(outdir+k+".faa")
                os.remove(outdir+k+".fna")
                log.info("None found passing coverage thresholds in potential core: %s"%k)
            if len([x for x in oarry[i] if round(x)==1])==len(orgs):
                singleslist.append(k)
        coredict["singles"] = singleslist
        gmfil.write("#Singles List:\t" + ",".join(singleslist) + "\n")
        log.info("Single genes found: %s"%singleslist)
    with open(outdir + "coreresults.json","w") as cfil:
        json.dump(coredict,cfil,indent=2)
    conn.close()
    return coredict


# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Extract core genes from sql database""")
    parser.add_argument("input", help="Hmm hits sql database (READS TABLES: HMMhits, Seqs)")
    parser.add_argument("-p", "--pct", help="percent match of gene length to consider as core hit (default: 0.9)",
                        type=float, default=0.7)
    parser.add_argument("-p2", "--pct2",
                        help="percent match of hmm model length to consider as core hit (default: 0.5)",
                        type=float, default=0.7)
    parser.add_argument("-e", "--evalue", help="Remove hits > evalue default(1e-4)", type=float, default=0.1)
    parser.add_argument("-b", "--bitscore", help="Remove hits < bitscore default(0)", type=float, default=0)
    parser.add_argument("-bh", "--besthit", help="Only write best hit hmm for each gene", action='store_true')
    parser.add_argument("-ts", "--topseq", help="Only write best hit hmm for each gene", action='store_true')
    parser.add_argument("-f2", "--filter", help="Use per-model bitscore cutoffs dictated in file (format: modelname value)", default=None)
    parser.add_argument("-od", "--outdir", help="Store results in OUTDIR (default: currentdir/results_xxx)", default=None)
    args = parser.parse_args()
    writeall(args.outdir, args.input, args.pct, args.pct2, args.evalue, args.bitscore, args.besthit, args.filter, args.topseq)