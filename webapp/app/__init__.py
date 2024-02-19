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

import os
from flask import Flask
from flask_mail import Mail

app = Flask(__name__)
configfile = os.environ.get('ARTS_SETTINGS',False)
if configfile and os.path.exists(configfile):
    app.config.from_envvar('ARTS_SETTINGS')
else:
    try:
        parentdir = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-2])
        if os.path.exists(os.path.join(parentdir,"config","activeconfig.conf")):
            app.config.from_pyfile(os.path.join(parentdir,"config","activeconfig.conf"), silent=True)
        elif os.path.exists("/etc/artsapp.conf"):
            app.config.from_pyfile("/etc/artsapp.conf", silent=True)
        elif os.path.exists(os.path.join(parentdir,"config","artsapp_default.conf")):
            app.config.from_pyfile(os.path.join(parentdir,"config","artsapp_default.conf"), silent=True)
        else:
            raise ImportError
    except ImportError:
        raise ImportError("Error: Could not find settings config file. Set ARTS_SETTINGS environment variable or add activeconfig.py to config folder")

#Check folder configs, defaults to /tmp if default dirs are not found
parentdir = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-3])
temp = app.config.get("RESULTS_FOLDER",os.path.join(parentdir,"results"))
if os.path.exists(temp):
    app.config["RESULTS_FOLDER"] = temp
else:
    app.config["RESULTS_FOLDER"] = "/tmp"
temp = app.config.get("UPLOAD_FOLDER",os.path.join(parentdir,"uploads"))
if os.path.exists(temp):
    app.config["UPLOAD_FOLDER"] = temp
else:
    app.config["UPLOAD_FOLDER"] = "/tmp"


mail = Mail(app)

from app import views
from app import models
from app import api