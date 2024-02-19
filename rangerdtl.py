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

import subprocess, setlog, os, multiprocessing as mp
from numpy import median

global log
log = setlog.init(toconsole=True)

def findnodebyname(orgname,intree):
    found = False
    for x in intree.iter_descendants():
        if x.name == orgname:
            found = x
            break
    return found

def gettransfers(rdtlrslt,sptree):
    orgrecs={}
    for k,rslt in rdtlrslt.items():
        hkgene,gid=k.split("_@_")
        for line in rslt.split("\n"):
            if "Transfer," in line:
                x = line.split(" --> ")
                #get receipient and doner
                dx = x[-2][:x[-2].index(",")]
                rx = x[-1]
                dxnode = findnodebyname(dx,sptree)
                rxnode = findnodebyname(rx,sptree)
                distdx = dxnode.get_distance(dxnode.get_common_ancestor(rxnode))
                distrx = rxnode.get_distance(rxnode.get_common_ancestor(dxnode))
                if dx not in orgrecs:
                    orgrecs[dx] = []
                if rx not in orgrecs:
                    orgrecs[rx] = []
                #Add Least common ancestor line
                lca = line.split("]")[0].split("[")[1]
                orgrecs[dx].append([hkgene,gid,dx,"D",rx,distdx,lca])
                orgrecs[rx].append([hkgene,gid,rx,"R",dx,distrx,lca])
    return orgrecs

def callrdtl(instr,k,dbg=False):
    #p = subprocess.Popen(("ranger-dtl-U",), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.Popen(("ranger-dtl-U",), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    out, err = p.communicate(instr)
    if dbg:
        with open(os.path.join(dbg,k+".rdtl"),"w") as fil:
            fil.write(out)
    return (out, err)

def rdtlresults(treedict,mcpu=1):
    results={}
    #if set to directory output raw ranger-dtl results
    dbg = os.getenv("DEBUG_RDTL",False)
    if dbg and os.path.isdir(dbg):
        log.debug("Found DEBUG_RDTL saving to: %s"%dbg)
    elif dbg:
        dbg = False
        log.debug("Found DEBUG_RDTL but directory is invalid. Set to a valid directory")
    if mcpu > 1:
        pool = mp.Pool(mcpu)
        for k,x in treedict.items():
            if type(x) is str:
                results[k]=pool.apply_async(callrdtl, args=(x,k,dbg))
        pool.close()
        pool.join()
        for x in results.keys():
            err = results[x].get()[1]
            results[x] = results[x].get()[0]
            if err:
                log.error("Tree %s: %s"%(x,err))
    else:
        for k,x in treedict.items():
            if type(x) is str:
                # results[k]=pool.apply_async(callrdtl, args=(x,k,dbg))
                results[k]=callrdtl(x,k,dbg)[0]
    return results

def filtertransfers(orgrecs,sptree,orgname,medfilt=False):
    orgleaf = sptree.get_leaves_by_name(orgname)[0]
    orglist = {x.name:"@%s"%orgleaf.get_distance(x) for x in orgleaf.get_ancestors()}
    orglist[orgleaf.name] = ""
    orgrecs2 = {}
    # orgrecs2={k:v for k,v in orgrecs.items() if any(k in x for x in orglist)}
    for k,v in orgrecs.items():
        if k in orglist.keys():
            md = median([val[-2] for val in v])
            if medfilt:
                temp = [val for val in v if val[-2] >= md]
                orgrecs2["%s%s"%(k,orglist[k])] = temp
            else:
                orgrecs2["%s%s"%(k,orglist[k])] = v
    return orgrecs2

# Get only transfers with query org in LCA, rx or dx
def getorgtransfers(rdtlrslt,sptree,qorg):
    orgrecs = {}
    orgleaf = sptree.get_leaves_by_name(qorg)[0]
    #Get list of ancestors and only use monophyletic
    temp = {anc.name:set(x.genus for x in anc.get_leaves() if x.genus) for anc in orgleaf.get_ancestors()}
    orglist = [anc for anc,v in temp.items() if len(v)<=1]
    orglist.append(orgleaf.name)
    for k,rslt in rdtlrslt.items():
        hkgene,gid=k.split("_@_")
        for line in rslt.split("\n"):
            if "Transfer," in line and qorg in line:
                x = line.split(" --> ")
                #get receipient and doner
                dx = x[-2][:x[-2].index(",")]
                rx = x[-1]
                #Filter only transfers within monphyletic clade based on genus value
                if dx in orglist or rx in orglist:
                    if dx not in orgrecs:
                        orgrecs[dx] = []
                    if rx not in orgrecs:
                        orgrecs[rx] = []
                    #Add Least common ancestor line
                    lca = line.split("]")[0].split("[")[1]
                    #Get LCA distance from root
                    x = lca.split(", ")
                    org1 = sptree.get_leaves_by_name(x[0])[0]
                    org2 = sptree.get_leaves_by_name(x[1])[0]
                    anc = org1.get_common_ancestor(org2)
                    lcadist = sptree.get_distance(anc)
                    orgrecs[dx].append([hkgene,gid,dx,"D",rx,lca,lcadist])
                    orgrecs[rx].append([hkgene,gid,rx,"R",dx,lca,lcadist])
    return orgrecs

def runallrdtl(sptree,orgname,treedict,mcpu=1):
    orgname = orgname.replace("-","_")
    results = rdtlresults(treedict,mcpu)
    return getorgtransfers(results,sptree,orgname)
    # orgrecs = gettransfers(results,sptree)
    # return filtertransfers(orgrecs,sptree,orgname)