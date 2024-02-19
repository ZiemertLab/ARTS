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

import argparse, sys, time, os, logging, shutil
from Daemon import Daemon
from logging.handlers import RotatingFileHandler
from redis import Redis
import artspipeline1

class dictobj(object):
    """Convert dictionary to object"""
    def __init__(self,dict):
        for k,v in dict.items():
            setattr(self,k,v)

class ArtsDaemon(Daemon):
    def __init__(self, pidfile, stdin='/dev/null', stdout='/dev/null', stderr='/dev/null', redis=None):
        cfgfile = os.environ.get('ARTS_SETTINGS',False)
        self.pdir = os.path.dirname(os.path.realpath(__file__))
        if not cfgfile and os.path.exists(os.path.join(self.pdir,"webapp","config","activeconfig.conf")):
            cfgfile = os.path.join(self.pdir,"webapp","config","activeconfig.conf")
        elif not cfgfile and os.path.exists("/etc/artsapp.conf"):
            cfgfile = "/etc/artsapp.conf"
        elif not cfgfile and os.path.exists(os.path.join(self.pdir,"webapp","config","artsapp_default.conf")):
            cfgfile = os.path.join(self.pdir,"webapp","config","artsapp_default.conf")
        else:
            print("Could not find config. Set 'ARTS_SETTINGS' env var to config location. ex:\nexport ARTS_SETTINGS='/home/user/artsapp.cfg'")
            exit(1)
        with open(cfgfile,'r') as fil:
            self.config = {x.split("=")[0].strip().upper():x.split("=",1)[1].strip().strip('"').strip("'") for x in fil if '=' in x and not x.startswith('#')}
        #check refdir and results folder default to local
        temp = self.config.get("REF_FOLDER",os.path.join(self.pdir,"reference"))
        if os.path.exists(temp):
            self.config["REF_FOLDER"] = temp
        else:
            print("Could not find reference folder. Check config")
            exit(1)
        temp = self.config.get("RESULTS_FOLDER",os.path.join(self.pdir,"results"))
        if os.path.exists(temp):
            self.config["RESULTS_FOLDER"] = temp
        else:
            print("Could not find results folder. Storing in /tmp")
            self.config["RESULTS_FOLDER"] = "/tmp"
        Daemon.__init__(self,pidfile)
        if not redis:
            redis = self.config.get("REDISURL","redis://localhost:6379/0")
        self.redis = Redis.from_url(redis)
        self.runningjob = None

        #Set logging
        self.log = logging.getLogger("artsdaemon")
        formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
        handler = RotatingFileHandler("%s.log"%pidfile, mode='a', maxBytes=5*1024*1024, backupCount=2, encoding=None, delay=0)
        handler.setFormatter(formatter)
        self.log.addHandler(handler)
        self.log.setLevel(logging.DEBUG)

    def getjobs(self,Q):
        if self.redis.llen(Q):
            jobid = self.redis.rpop(Q)
            jobargs = self.redis.hgetall("artsjob:%s"%jobid)
            return jobid,jobargs
        return False,False

    def clean(self):
        RQ = self.redis.lrange("RQ",0,-1)
        SQ = self.redis.lrange("SQ",0,-1)
        PQ = self.redis.lrange("PQ",0,-1)
        for job in self.redis.keys():
            if "artsjob:" in job:
                jobid = job.split(":")[-1]
                if jobid not in (RQ+SQ+PQ):
                    self.redis.delete(job)
                    # self.redis.lrem("DQ",jobid)
                    self.redis.lrem("EQ",jobid)

        now = time.time()
        def removeall(maxage,results,dir=True):
            if f and len(results):
                for result in results:
                    #skip example
                    if "example" in os.path.split(result)[-1]:
                        continue
                    if (now - os.path.getmtime(result))/(60*60*24) > maxage:
                        if dir:
                            shutil.rmtree(result)
                        else:
                            os.remove(result)

        #Check and clean old results
        f = self.config.get("RESULTS_FOLDER",False)
        results = [os.path.join(f,x) for x in os.listdir(f) if os.path.isdir(os.path.join(f,x))]
        maxage = self.config.get("RESULT_AGE",30)
        removeall(maxage,results)

        #Check and clean old uploads
        f = self.config.get("UPLOAD_FOLDER",False)
        results = [os.path.join(f,x) for x in os.listdir(f) if os.path.isdir(os.path.join(f,x))]
        maxage = self.config.get("RESULT_AGE",30)
        removeall(maxage,results)

        #Check and clean old archived results
        f = self.config.get("ARCHIVE_FOLDER",False)
        results = [os.path.join(f,x) for x in os.listdir(f) if x.endswith(".zip")]
        maxage = self.config.get("ARCHIVE_AGE",100)
        removeall(maxage,results,dir=False)

    def pause(self):
        """Move unstarted to temporary pause Queue"""
        for x in self.redis.lrange("SQ",0,-1):
            strt = self.redis.hget("arstjob:%s"%x,"started")
            if not strt:
                self.redis.lpush("PQ",x)
                self.redis.lrem("SQ",x)

    def resume(self):
        """Move unstarted to temporary pause Queue"""
        for x in self.redis.lrange("PQ",0,-1):
            self.redis.rpush("SQ",x)
            self.redis.lrem("PQ",x)

    def info(self):
        report = "All keys: %s\nStart Queue: %s\nPause Queue: %s\nError Queue: %s\nDone: %s\n"%\
                 (self.redis.keys(),self.redis.llen("SQ"),self.redis.llen("PQ"),self.redis.llen("EQ"),self.redis.llen("DQ"))
        return report

    def run(self):
        jobid = ""
        #Stop looping if pid file is removed
        while os.path.exists(self.pidfile):
            try:
                if self.redis.llen("SQ"):
                    jobid, jobargs = self.getjobs("SQ")
                    if jobid:
                        self.log.info("Started %s"%jobid)
                        self.redis.hset("artsjob:%s"%jobid,"started",int(time.time())) # Mark start time (epoch time)
                        self.redis.lpush("RQ",jobid)
                        self.runningjob = jobid
                        #### Options from job arguments ####
                        options = jobargs.get("options",False)
                        asrun = jobargs.get("asrun",False)
                        aspath = self.config.get("ANTISMASH_PATH",os.path.join(self.pdir,"antismash"))
                        if not os.path.exists(aspath):
                            aspath = False
                        ## Run bigscape if this is a multianalysis job
                        run_bsc = False
                        if "," in jobargs["infile"]:
                            run_bsc = True
                        ## This path requires more attention at the result part
                        bcp = self.config.get("BIG_SCAPE_PATH",os.path.join(self.pdir,"BiG-SCAPE-master"))
                        if not os.path.exists(bcp):
                            bcp = False
                        knownhmms = os.path.join(self.config["REF_FOLDER"],"knownresistance.hmm")
                        if not os.path.exists(knownhmms):
                            knownhmms=False
                        dufhmms = os.path.join(self.config["REF_FOLDER"],"dufmodels.hmm")
                        if not os.path.exists(dufhmms):
                            dufhmms=False

                        custhmms = jobargs.get("custmdl",False)
                        if not os.path.exists(custhmms):
                            custhmms = False
                        custcorehmms = jobargs.get("custcoremdl",False)
                        if custcorehmms and not os.path.exists(custcorehmms):
                            custcorehmms = None
                        rnahmm = os.path.join(self.config["REF_FOLDER"],"barnap_bact_rRna.hmm")
                        if not os.path.exists(rnahmm):
                            rnahmm = None
                        ## Do the job
                        argdict = {"input":jobargs["infile"], "refdir":os.path.join(self.config["REF_FOLDER"],jobargs["ref"]), "hmmdblist":None, "orgname":None,
                                   "tempdir":self.config.get("TEMPDIR",None), "resultdir":os.path.join(self.config["RESULTS_FOLDER"],jobargs["id"]), "prebuilttrees":False,
                                   "rnahmmdb":rnahmm, "thresh":jobargs.get("cut","TC"), "astral":self.config.get("ASTJAR",False), "toconsole":False,
                                   "multicpu":self.config.get("MCPU",1), "runantismash":asrun, "knownhmms":knownhmms, "dufhmms":dufhmms, "custcorehmms":custcorehmms,
                                   "customhmms":custhmms, "antismashpath":aspath, "options":options, "bigscapepath":bcp, "runbigscape":run_bsc}
                        argobj = dictobj(argdict)
                        artspipeline1.call_startquery(argobj)

                        #artspipeline1.startquery(infile=jobargs["infile"],refdir=os.path.join(self.config["REF_FOLDER"],jobargs["ref"]),hmmdbs=None,td=self.config.get("TEMPDIR",None),
                        #                         rd=os.path.join(self.config["RESULTS_FOLDER"],jobargs["id"]),rnahmm=rnahmm,cut=jobargs.get("cut","TC"),
                        #                         astjar=self.config.get("ASTJAR",False),toconsole=False,mcpu=self.config.get("MCPU",1),asrun=asrun,knownhmms=knownhmms,dufhmms=dufhmms,
                        #                         custcorehmms=custcorehmms,custhmms=custhmms,aspath=aspath,options=options, bcp=bcp, run_bsc = run_bsc)

                        self.runningjob = False
                        self.redis.lrem("RQ",jobid)
                        self.redis.lpush("DQ",jobid)
                        self.redis.hset("artsjob:%s"%jobid,"finished",int(time.time())) # Mark finish time (epoch time)
                        self.log.info("Finished %s"%jobid)
                        self.log.info("Compressing job to Archive %s"%jobid)
                        try:
                            ad = self.config.get("ARCHIVE_FOLDER","/tmp")
                            shutil.make_archive(os.path.join(ad,str(jobid)),"zip",os.path.join(self.config["RESULTS_FOLDER"],str(jobid)))
                            self.log.info("Archived %s at %s"%(jobid,ad))
                        except Exception as e:
                            self.log.error("Failed to make Archive %s"%jobid)
                            self.log.exception("exception")
            except Exception as e:
                self.log.error("Unexpected error: %s"%e)
                self.log.exception("exception")
                #move to error queue
                if self.runningjob:
                    self.redis.hset("artsjob:%s"%self.runningjob,"error",int(time.time())) # Mark error
                    self.redis.lrem("RQ",self.runningjob)
                    self.redis.lpush("EQ",self.runningjob)
            time.sleep(3)
        #If finished send exit
        self.log.info("Pidfile not found, finishing last job and exiting")
        exit(0)

def rundaemon(action, redis, pidfile, cpu=None, resultage=30, archiveage=100):
    artsdmn = ArtsDaemon(pidfile,redis=redis)
    if resultage:
        artsdmn.config["RESULT_AGE"] = resultage
    if archiveage:
        artsdmn.config["ARCHIVE_AGE"] = archiveage
    if cpu:
        artsdmn.config["MCPU"] = cpu
    elif os.environ.get('ARTS_CPU',False):
        artsdmn.config["MCPU"] = os.environ.get('ARTS_CPU',1)
    if action == "run":
        # Run without forking - use this with systemd / daemon manager
        with open(pidfile,"w") as fil:
            fil.write("%s\n"%os.getpid())
        artsdmn.run()
    if action == "start":
        artsdmn.start()
    elif action == "stop":
        artsdmn.stop()
    elif action == "restart":
        artsdmn.restart()
    elif action == "clean":
        artsdmn.clean()
    elif action == "pause":
        artsdmn.pause()
    elif action == "resume":
        artsdmn.resume()
    elif action == "info":
        print(artsdmn.info())
    sys.exit(0)

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""ARTS workflow daemon. Consumes jobs from redis store""")
    parser.add_argument("input", help="Action for daemon (start|stop|restart|info|clean|pause|resume)", choices=["start","stop","restart","info","clean","pause","resume","run"])
    parser.add_argument("-rd", "--redis", help="Redis store url ( default: from config file )", default=None)
    parser.add_argument("-ra", "--resultage", help="Age to keep oldest result jobs in days ( default: 30 )",type=int, default=30)
    parser.add_argument("-aa", "--archiveage", help="Age to keep oldest archived jobs in days ( default: 100 )",type=int, default=100)
    parser.add_argument("-pid", "--pidfile", help="Process id file (default: /tmp/artsdaemon-1.pid)", default="/tmp/artsdaemon-1.pid")
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing (default: from config file)", default=None)
    args = parser.parse_args()
    rundaemon(args.input,args.redis,args.pidfile,args.multicpu,args.resultage,args.archiveage)