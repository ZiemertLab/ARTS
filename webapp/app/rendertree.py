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

from ete3 import Tree
from app import app
from PIL import Image, ImageDraw, ImageFont
import argparse, os

#FAILOVER when PyQT cannot import
def getsize(txt,fnt,padding=0):
    im = Image.new("RGB", (1,1), "white")
    d = ImageDraw.Draw(im)
    return tuple([x+padding for x in list(d.textsize(txt,fnt))])

def getcoords(txt,fnt,findtxt,offsety=0):
    global test
    coordlist = []
    height = 0
    if findtxt in txt:
        txtlist = txt.split(findtxt)
        test = txtlist
        for p in txtlist[:-1]:
            p1,p2 = tuple(p.rsplit("\n",1))
            p2size = getsize(p2,fnt)
            p1size = getsize(p1,fnt)
            height = p2size[1]
            coordlist.append((p2size[0],p1size[1]+offsety))
            offsety = offsety + p1size[1]
    return coordlist,height

def printpng(fname,txt,spname,fontfil,fontsize,offset=0):
    fnt = ImageFont.truetype(fontfil, fontsize)
    imsize = getsize(txt,fnt,padding=offset)
    im = Image.new("RGB", imsize, "white")
    d = ImageDraw.Draw(im,"RGBA")
    coords,height = getcoords(txt,fnt,spname,offset)
    for c in coords:
        d.rectangle([c,(imsize[0],c[1]+height)],fill=(0,200,0,50))
    d.text((offset,offset),txt,font=fnt,fill="black")
    with open(fname,"w") as fil:
        im.save(fname,"png")

def rendertree(infile, width=800, spname=False):
    fntfile = os.path.join(os.path.split(os.path.realpath(__file__))[0],"fonts","roboto.ttf")
    if infile and os.path.exists(infile):
        T = Tree(infile)
        if app.config.get("DISABLE_PYQT",False):
            printpng(infile+".png",str(T),spname,fntfile,16,offset=10)
            return True
        else:
            try:
                if spname:
                    for node in T:
                        if "OUTGROUP" in node.name:
                            T.set_outgroup(node)
                            T.ladderize()
                        if spname in node.name:
                                #Use styles if availible (PyQT) otherwise use TEXT
                                from ete3 import NodeStyle
                                nstyle = NodeStyle()
                                nstyle["bgcolor"] = "#00CC00"
                                nstyle["size"] = 10
                                node.set_style(nstyle)
                    T.render(infile+".png", w=width, units="px")
                    return True
            except ImportError:
                printpng(infile+".png",str(T),spname,fntfile,16,offset=10)
                return True
    else:
        return False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Convert Newick files to PNG""")
    parser.add_argument("-in", "--infile", help="Newick file to use", default=None)
    parser.add_argument("-w", "--width", help="Width in pixels to output", type=int, default=800)
    parser.add_argument("-spname", "--speciesname", help="Highlight nodes matching species name", default=False)
    args = parser.parse_args()
    rendertree(args.infile,args.width,args.speciesname)