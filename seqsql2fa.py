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

import argparse, setlog, datetime, os, sqlite3 as sql

global log
log = setlog.init(toconsole=True)

def writefasta(db,outfile,nuc=False,idprfx=""):
    conn = sql.connect(db)
    cur = conn.cursor()
    fn,ext = os.path.splitext(outfile)
    fname = fn+ext
    cur.execute("SELECT * FROM Seqs")
    rec = cur.fetchone()
    if rec:
        log.info("Writing sequences to disk...")
        with open(fname,"w") as ofil:
            r = rec
            ofil.write(">%s|%s%s %s|source|%s|loc|%s_%s_%s\n%s\n" % (r[1], idprfx, r[0],r[2],r[4],r[5],r[6],r[7],r[-1-int(nuc)]))
            while rec:
                for r in cur:
                    ofil.write(">%s|%s%s %s|source|%s|loc|%s_%s_%s\n%s\n" % (r[1], idprfx, r[0],r[2],r[4],r[5],r[6],r[7],r[-1-int(nuc)]))
                else:
                    rec=False
    else:
        log.warning("No sequences to export for database %s"%db)
    #Update lastscan of exported seqs
    # try:
    #     cur.execute("UPDATE Seqs SET lastscan=? WHERE lastscan <= ?",(timestamp,fus))
    # except sql.OperationalError as ex:
    #     print "Error:",ex
    conn.commit()
    conn.close()


# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write fasta file of sequence DB with format >ORG|dbkey")
    parser.add_argument("input", help="sequence database")
    parser.add_argument("outfile", help="Filename of for saved files. Number extension will be added for split files.")
    parser.add_argument("-n","--nuc", help="Write DNA sequences instead of protein",action='store_true')
    args = parser.parse_args()
    writefasta(args.input,args.outfile,args.nuc)