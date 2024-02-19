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

import re, os
from ete3 import Tree, EvolTree

def maketree(x):
    with open(x,"r") as fil:
        x = fil.readline()      #x = fil.next()
    #Replace bad chars from newick
    x = re.sub("\[I\d+?\]|\"|'|QUERY___","",x)
    T = EvolTree(x)
    T.set_species_naming_function(lambda spname: spname.split("|")[0])
    T.resolve_polytomy(recursive=True) ##Make bifurcated
    return T

def relabeltree(T,prfx,outgroup=None,qorg=None):
    i=0
    og=0
    T=T.copy("deepcopy")
    T.hasroot=False
    if not len(T.name):
        T.name=str(prfx)+str(i)
        i+=1
    for node in T.iter_descendants("postorder"):
        if not len(node.name):
            node.name=str(prfx)+str(i)
            i+=1
        else:
            y=node.name.split("|")
            node.name=y[0].replace("-","_")
            if len(y)>1:
                node.gid=y[1]
            #set genus of all reference leafs except unknown query if qorg is assigned
            if qorg and qorg not in node.name:
                #simply take the genus name separated by underscore
                node.genus = node.name.split("_")[1]
            else:
                node.genus = None
        if "OUTGROUP" in node.name or (outgroup and outgroup in node.name):
            outgroup = node.name
            og+=1
    if og==1:
        T.set_outgroup(outgroup)
        T.hasroot=True
    return T

def splitdups(T,org,spset,outgroup=None):
    trees={}
    #T.collapse_lineage_specific_expansions() #Condense expansions
    for node in T:
        if org in node.name:
            gid=node.name.split("|")[-1]
            trees[gid]=relabeltree(T,"G",outgroup)
    for k,tree in trees.items():
        for node in tree:
            if (org in node.name and k not in node.gid) or (node.name not in spset):
                node.delete()
    return trees

def writetree(T):
    # To ensure readable branch lengths and internal nodes for ranger-dtl
    # Also latterize to ensure constant ordering
    T.ladderize()
    return Tree(T.write(format=3),format=3).write(format=3,dist_formatter="%0.16f")

# def renderpng(infile, width, spname):
#     T = Tree(infile)
#     nstyle = NodeStyle()
#     nstyle["bgcolor"] = "LightSteelBlue"
#     nstyle["size"] = 10
#
#     for node in T:
#         if spname in node.name:
#             node.set_style(nstyle)
#
#     T.render(infile+".png", w=width, units="px")

def mergetrees(sptree,tlist,org,outgroup=None):
    sptree=maketree(sptree)
    sptree=relabeltree(sptree,"SP",outgroup,qorg=org)
    treedict = {}
    # treedict2 = {}
    sptline = writetree(sptree)

    for HKtree in tlist:
        hk = os.path.splitext(os.path.split(HKtree)[-1])[0] #get from filename
        temp = maketree(HKtree)
        temp = splitdups(temp,org,sptree.get_species(),outgroup)
        for gid,x in temp.items():
            pfx=""
            if not x.hasroot:
                pfx="[&U]"
            treedict[hk+"_@_"+gid] = "%s\n%s%s\n" % (sptline,pfx,writetree(x))
            # treedict2[hk+"_@_"+gid] = x
    return (treedict,sptree)

def splitdups2(T,org,sptree,outgroup=None):
    trees={}
    sptrees={}
    #T.collapse_lineage_specific_expansions() #Condense expansions
    for node in T:
        if org in node.name:
            gid=node.name.split("|")[-1]
            trees[gid]=relabeltree(T,"G",outgroup)
            sptrees[gid]=sptree.copy(method="deepcopy")
    for k,tree in trees.items():
        for node in tree:
            if (org in node.name and k not in node.gid) or (node.name not in sptree.get_species()):
                node.delete()
    for k,tree in sptrees.items():
        for node in tree:
            if node.name not in trees[k].get_species():
                node.delete()
    return trees,sptrees