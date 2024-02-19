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

import os, argparse, setlog, tempfile, sqlite3 as sql

global log
log = setlog.init(toconsole=True)

def writeall(db,outdir,allcopy):
	conn=sql.connect(db)
	cur=conn.cursor()

	rnagenes=[x[0] for x in cur.execute("SELECT DISTINCT hmmhit FROM RNAhits")]

	#Combine tables into view and write full hits, fragment hits, genematrix
	cur.execute("CREATE TEMP VIEW allgenes AS SELECT A.*,B.* FROM RNAhits A INNER JOIN Seqs B ON A.seqid=B.seqid")
	for i,rna in enumerate(rnagenes):
		cur.execute("SELECT orgname,GROUP_CONCAT(seqid),GROUP_CONCAT(naseq) FROM allgenes WHERE hmmhit=? GROUP BY orgname",(rna,))
		with open(outdir+"RNA_"+rna+".fna","w") as nafil:
			for x in cur:
				seqlist=x[2].split(",")
				idlist=x[1].split(",")
				if allcopy:
					for j,s in enumerate(seqlist):
						nafil.write(">%s|%s\n%s\n"%(x[0],idlist[j],s))
				else:
					# SORT BY Longest then Most Common
					seqlist=sorted([(y,seqlist.count(y),idlist[seqlist.index(y)]) for y in set(seqlist)],key=lambda k: (len(k[0]),k[1]),reverse=True)
					nafil.write(">%s|%s\n%s\n"%(x[0],str(seqlist[0][2]),str(seqlist[0][0])))
		log.info("Finished RNA_"+rna+".fna")
	conn.close()

def runfile(db, allcopy=False, outdir=None):
	if db and os.path.exists(db):
		#correct bad inputs
		if not outdir:
			outdir = tempfile.mkdtemp(prefix="results_",dir="./")
			log.warning("NO RESULT FOLDER SPECIFIED. Results will be stored in temporary folder: "+outdir)
		elif not os.path.exists(outdir):
			os.mkdir(outdir)
		if not outdir.endswith("/"):
			outdir += "/"
		log.info("Extracting sequences...")
		writeall(db,outdir,allcopy)
	else:
		log.error("No database found see --help for usage")

# Commandline Execution
if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="""Extract RNA genes from sql database""")
	parser.add_argument("input", help="Sql database (READS TABLES: Seqs, RNAhits)")
	parser.add_argument("-all", "--allcopy", help="Write all copies from same org (default:False, writes most common copy or first copy)", action='store_true')
	parser.add_argument("-od", "--outdir", help="Store results in OUTDIR (default: currentdir/results_xxx)",default=None)
	args = parser.parse_args()
	runfile(args.input,  args.allcopy, args.outdir)