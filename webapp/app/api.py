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

from app import app
from app import routines
from flask import request
from flask_restful import Resource, Api, reqparse
from werkzeug import datastructures
import os

api = Api(app)
parser = reqparse.RequestParser()
parser.add_argument('seqfile',type=datastructures.FileStorage,help='Input genbank file without antismash annotation',location='files')
parser.add_argument('asseqfile',type=datastructures.FileStorage,help='Input genbank file with existing antismash annotation',location='files')
parser.add_argument('asjobid',type=str,help='Existing antismash jobid',location="form")
parser.add_argument('reference',type=str,help='Reference genomes set',location="form",default="actinobacteria")

testdict = {}

class api_v1p0(Resource):
    def get(self):
        return {"server_status":routines.getservstats()}

#API endpoints for job tables and reports
class api_v1p0_jobreports(Resource):
    def getresultfile(self,fname):
        rf = os.path.join(self.rd,"tables",fname)
        if self.jobid and routines.checkresult(self.jobid) and os.path.exists(rf):
            return rf
        else:
            return False

    def act_jobstatus(self):
        if self.jobid and routines.checkresult(self.jobid):
            return routines.getjobstatus(self.jobid)
        return {"data":[]}

    def act_krtab(self):
        rf = self.getresultfile("knownhits.json")
        if rf:
            return routines.json2dict(rf)
        else:
            return {"data":[]}

    def act_duftab(self):
        rf = self.getresultfile("dufhits.json")
        if rf:
            return routines.json2dict(rf)
        else:
            return {"data":[]}

    def act_dupmatrix(self):
        rf = self.getresultfile("duptable.json")
        if rf:
            return routines.json2dict(rf)
        else:
            return {"data":[]}

    def act_bgctable(self):
        rf = self.getresultfile("bgctable.json")
        if rf:
            return routines.json2dict(rf)
        else:
            return {"data":[]}

    def act_funcstats(self):
        rf = self.getresultfile("coretable.json")
        if rf:
            return routines.json2dict(rf)["funcstats"]
        else:
            return {"data":[]}

    def act_summarytab(self):
        rf = self.getresultfile("coretable.json")
        if rf:
            return routines.json2dict(rf)
        else:
            return {"data":[]}

    def get(self,jobid,action):
        self.jobid = jobid
        self.action = action
        self.rd = os.path.join(app.config['RESULTS_FOLDER'],self.jobid)

        #alias endpoints (act_"endpoint"):
        self.act_sumtab = self.act_sumtable = self.act_summarytable = self.act_summarytab

        try:
            return getattr(self,"act_%s"%action.lower().replace("__",""))()
        except AttributeError:
            return {"error":"Invalid parameters. Usage api/<jobid>/<action>","availible_actions":[x.replace("act_","") for x in dir(self) if x.startswith("act_")]}

class api_v1p0_startjob(Resource):
    def post(self):
        version = "1.0"
        args = parser.parse_args()
        args.asseqfile = str(args.asseqfile)
        args.seqfile = str(args.seqfile)
        filename,asrun = routines.getinfile()
        if filename:
            artsjob = routines.addjob(infile=filename,ref="actinobacteria",asrun=asrun)
            return {"parsed_request":args,"jobid":artsjob.id,"status":"success","version":version,"report_url":"http://"+request.url.split("/")[2]+"/results/"+artsjob.id}, 201
        elif not asrun == "invalid":
            return {"parsed_request":args,"error":"Antismash run id is invalid","status":"error","version":version}, 400
        elif not asrun == "pending":
            return {"parsed_request":args,"error":"Antismash job is not finished","status":"error","version":version}, 400
        else:
            return {"parsed_request":args,"status":"failed","version":version}, 400

api.add_resource(api_v1p0,"/api/1.0","/api/1.0/","/api","/api/","/api/1.0/serverstatus","/api/serverstatus")
api.add_resource(api_v1p0_startjob,"/api/1.0/startjob","/api/startjob")
api.add_resource(api_v1p0_jobreports,"/api/1.0/job/<string:jobid>/<string:action>","/api/job/<string:jobid>/<string:action>","/results/<string:jobid>/<string:action>")
