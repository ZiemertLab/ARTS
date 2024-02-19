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
import logging

def init(logfile=None,toconsole=False,level="debug",logname=""):
    if logname:
        log=logging.getLogger(logname)
    else:
        log=logging.getLogger("root")
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
    if logfile:
        log.handlers = [] #clear all handlers
        handler=logging.FileHandler(logfile)
        handler.setFormatter(formatter)
        log.addHandler(handler)
    if toconsole:
        handler=logging.StreamHandler()
        handler.setFormatter(formatter)
        log.addHandler(handler)
    if level.upper()=="CRITICAL":
        log.setLevel(logging.CRITICAL)
    elif level.upper()=="ERROR":
        log.setLevel(logging.ERROR)
    elif level.upper()=="WARNING":
        log.setLevel(logging.WARNING)
    elif level.upper()=="INFO":
        log.setLevel(logging.INFO)
    elif level.upper()=="DEBUG":
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.NOTSET)
    return log