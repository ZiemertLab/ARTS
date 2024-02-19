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

def getcore(db, ev, bs, thrsh=0.95, pct=0.9, pct2=0.5, maxcopy=10.0, bh=False, limitorgs=[]):
    conn = sql.connect(db)
    cur = conn.cursor()
    orgs = set(str(x[0]) for x in cur.execute("SELECT DISTINCT orgname FROM HMMhits"))
    # get core genes based on orgcount and hmmcoverage
    addwhere = ""
    if len(limitorgs):
        addwhere += " AND orgname IN " + str(limitorgs).replace("[", "(").replace("]", ")")
        orgs = set(limitorgs)
    if bh:
        addwhere += " AND flags>0"
    cur.execute(
        "CREATE TEMP VIEW coreset AS SELECT hmmhit,orgname,seqid from HMMhits WHERE genecov >= %.4f AND hmmcov >= %.4f AND evalue<=%s AND score>=%s" % (pct,pct2, ev, bs) + addwhere)
    sqlcmd = """
		SELECT t.hmmhit FROM
		(SELECT hmmhit,COUNT (DISTINCT orgname) AS orgcount,COUNT (DISTINCT seqid) AS genecount FROM coreset
		GROUP BY hmmhit) t WHERE t.orgcount >= %d AND genecount*1.0/orgcount <= %.4f""" % (round(thrsh * len(orgs)), maxcopy)
    hkcore = list(str(x[0]) for x in cur.execute(sqlcmd))
    conn.close()
    return hkcore, orgs

def getseqs(row):
    spots = set()
    spots2 = set()
    for x in row[-3].split('&'):
        xt = ast.literal_eval(x)
        spots |= set(range(xt[0] - 1, xt[1]))                  #spots |= set(xrange(xt[0] - 1, xt[1]))
        spots2 |= set(range((xt[0] - 1) * 3, xt[1] * 3))       #spots2 |= set(xrange((xt[0] - 1) * 3, xt[1] * 3))
    naary = np.fromstring(row[-2], dtype="|S1")
    aaary = np.fromstring(row[-1], dtype="|S1")
    # spots2=np.array(list(spots|set(np.array(list(spots))*2)|set(np.array(list(spots))*3)))
    spots = np.array(list(spots))
    spots2 = np.array(list(spots2))
    na = "".join(naary[spots2])
    aa = "".join(aaary[spots])
    return na, aa

def writeall(outdir, db, hkcore=None, pct=0.7, pct2=0.7, ev=0.1, bs=0, bh=True, wg=True, orgs=None, filt2=None):
    conn = sql.connect(db)
    cur = conn.cursor()
    if not hkcore:
        hkcore=[str(x[0]) for x in cur.execute("SELECT DISTINCT hmmhit FROM HMMhits")]
    if not orgs:
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
        # results = cur.execute(
        #     "SELECT orgname,seqid,genecov,hmmcov,MAX(iscore),genelen,GROUP_CONCAT('('||geneStart||','||geneEnd||')','&'),naseq,aaseq FROM allgenes WHERE hmmhit=? GROUP BY seqid",
        #     (hkgn,)).fetchall()
        results = cur.execute(
            "SELECT orgname,seqid,genecov,hmmcov,MAX(iscore),genelen,source,loc_start,loc_end,loc_strand,description,score,naseq,aaseq FROM allgenes WHERE hmmhit=? GROUP BY seqid",
            (hkgn,)).fetchall()
        with open(outdir + hkgn + ".faa", "w") as aafil, open(outdir + hkgn + ".fna", "w") as nafil:
            for x in results:
                # add to genematrix
                j = orgs.index(x[0])
                if wg:
                    seqs = (x[-2], x[-1])
                else:
                    seqs = getseqs(x)
                # separate fragments and full length
                if filt2 and hkgn in filt2.keys() and float(filt2[hkgn]) > float(x[4]):
                    continue
                if x[2] >= pct or x[3] >= pct2:
                    oarry[i, j] += 1
                    seqid = str(x[1])
                    aafil.write(">%s|%s\n%s\n" % (x[0], seqid, seqs[1]))
                    nafil.write(">%s|%s\n%s\n" % (x[0], seqid, seqs[0]))
                    #Add to core list (seqid, hkgene, source, location start, location end, strand, description, bitscore)
                    if hkgn not in coredict["core"].keys():
                        coredict["core"][hkgn] = {"seqs":[]}
                    if seqid not in coredict["seqs"].keys():
                        coredict["seqs"][seqid] = [seqid,hkgn,str(x[6]),x[7],x[8],x[9],x[10],x[11]]
                    coredict["core"][hkgn]["seqs"].append(seqid)
                # elif x[2] >= pct2:
                #     with open(outdir + hkgn + ".faa.frag", "w") as aafragfil, open(outdir + hkgn + ".fna.frag","w") as nafragfil:
                #         oarry2[i, j] += float(x[3])
                #         aafragfil.write(">%s|%s\n%s\n" % (x[0], str(x[1]), seqs[1]))
                #         nafragfil.write(">%s|%s\n%s\n" % (x[0], str(x[1]), seqs[0]))
        log.info("Wrote (%d of %d): %s" % (i + 1, len(hkcore), hkgn))
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


def run(db, bh=True, wg=True, co=False, pct=0.7, pct2=0.7, maxcopy=10.0, thrsh=0.95, ev=0.1, bs=0, lo=None, outdir=None, filt2=None):
    if db and os.path.exists(db):
        # correct bad inputs
        if thrsh > 1 or thrsh <= 0:
            thrsh = 0.95
            log.warning("Threshold out of range setting to default 0.95")
        if pct2 > 1 or pct2 < 0:
            pct2 = 0.2
            log.warning("Fragment coverage threshold out of range setting to default 0.2")
        if pct > 1 or pct < 0:
            pct = 0.9
            log.warning("Percent coverage out of range setting to default 0.9")
        if filt2 and os.path.isfile(filt2):
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
        limitlist = []
        if lo and os.path.exists(lo):
            with open(lo, "r") as fil:
                for line in fil:
                    limitlist.append(line.strip())
        hkcore, orgs = getcore(db, ev, bs, thrsh, pct, pct2, maxcopy, bh, limitlist)
        log.info("Found %d core genes: %s" % (len(hkcore), hkcore))
        if not co:
            log.info("Extracting sequences...")
            writeall(outdir, db, hkcore, pct, pct2, ev, bs, bh, wg, orgs, filt2)
    else:
        log.error("No sequence database found see --help for usage")


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
    parser.add_argument("-t", "--thrsh",
                        help="Threshold to calculate core gene set (percent of taxa with gene) default: 0.95",
                        type=float, default=0.95)
    parser.add_argument("-mc", "--maxcopy",
                        help="Threshold exclude core genes that have > maxcopy number of genes on average (default: 10.0)",
                        type=float, default=10.0)
    parser.add_argument("-e", "--evalue", help="Remove hits > evalue default(1e-4)", type=float, default=0.1)
    parser.add_argument("-b", "--bitscore", help="Remove hits < bitscore default(0)", type=float, default=0)
    parser.add_argument("-bh", "--besthit", help="Only write best hit hmm for each gene", action='store_true')
    parser.add_argument("-wg", "--wholegene", help="Write entire gene instead of match location", action='store_true')
    parser.add_argument("-co", "--coreonly", help="Only show core list, do not extract and write sequences",
                        action='store_true')
    parser.add_argument("-lo", "--limitorgs", help="Limit pan-genome to orgs in file", default=[])
    parser.add_argument("-f2", "--filter", help="Use per-model bitscore cutoffs dictated in file (format: modelname value)", default=None)
    parser.add_argument("-od", "--outdir", help="Store results in OUTDIR (default: currentdir/results_xxx)",
                        default=None)
    args = parser.parse_args()
    run(args.input, args.besthit, args.wholegene, args.coreonly, args.pct, args.pct2, args.maxcopy, args.thrsh,
            args.evalue, args.bitscore, args.limitorgs, args.outdir, args.filter)
