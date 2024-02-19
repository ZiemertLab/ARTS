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

import uuid

class ArtsJob(object):
    def __init__(self,**kwargs):
        self.infile = kwargs.get("infile",'')
        self.asjob = kwargs.get("asjob",'')
        if self.asjob:
            self.id = self.asjob
        else:
            self.id = kwargs.get("jobid",unicode(uuid.uuid4()))
        self.asrun = kwargs.get("asrun","")
        self.ref = kwargs.get("ref",'')
        self.cut = kwargs.get("cut","")
        self.custmdl = kwargs.get("custmdl",'')
        self.custcoremdl = kwargs.get("custcoremdl",'')
        self.options = kwargs.get("options",'')
        self.started = ''
        self.finished = ''
        self.worker = ''

    def getdict(self):
        return self.__dict__