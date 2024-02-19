
import os, threading, shutil, json, re
#from rendertree import rendertree
from app import rendertree
from flask import render_template, jsonify, Response, url_for, request, redirect, abort, make_response, send_from_directory, flash, copy_current_request_context
from app import app
from app import routines

@app.route('/')
@app.route('/index')
def index():
    return render_template("front.html")

@app.route('/upload', methods=['POST'])
def upload():
    #print request.files.getlist("seqfile1")
    filename,asrun = routines.getinfile()
    name = [os.path.split(x)[-1] for x in filename]
    #    filename = routines.getNCBIgbk(request.form["ncbiacc1"])
    return json.dumps({"filename": filename,"name": name})

@app.route('/results/<jobid>/trees/<tree>')
def gettreepng(jobid,tree,qorg=None):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
    treedir = os.path.join(rd,"trees")
    fname = [x for x in os.listdir(treedir) if (tree.lower()+".tree")==x.lower()]
    if len(fname):
        fname = fname[0]
    else:
        fname = False
    qorg = routines.getqorg(rd)

    if fname and os.path.exists(os.path.join(treedir,fname)):
        if os.path.exists(os.path.join(treedir,fname+".png")):
            return send_from_directory(treedir,fname+".png",mimetype='image/png')
        elif rendertree.rendertree(os.path.join(treedir,fname),800,qorg):
            return send_from_directory(treedir,fname+".png",mimetype='image/png')
    return app.send_static_file("images/blank.png")

@app.route('/static/exampleinputs/<fname>')
def sendexampleinputs(fname):
    return app.send_static_file(os.path.join("exampleinputs",fname))

@app.route('/analyze/')
@app.route('/analyze')
def analyze():
    return render_template("analyze.html")

@app.route('/results/')
@app.route('/results')
def results():
    getall = app.config.get("SHOW_ALL_RESULTS",False)
    if type(getall) == str:
        getall = True
    if getall:
        lastresult = routines.getallresults()
    else:
        lastresult = routines.getlastresults()
    return render_template("results.html",recentresults=lastresult)

@app.route('/serverstatus')
def serverstatus():
    status = routines.getservstats()
    return jsonify(status)

@app.route('/results/<jobid>')
@app.route('/results/<jobid>/')
def gotoreport(jobid):
    return redirect("/results/%s/report"%jobid)

@app.route('/results/<jobid>/report')
def showresult(jobid):
    if jobid and routines.checkresult(jobid):
        # resp = make_response(render_template("showresult.html",resid=jobid))
        # saverun = request.form.get('saverun',False)
        # if saverun and saverun.lower() != "false":
        #     lastresult = request.cookies.get('arts.lastresult')
        #     if lastresult and jobid not in lastresult:
        #         lastresult += ";%s"%jobid
        #     elif not lastresult:
        #         lastresult = jobid
        #     resp.set_cookie('arts.lastresult',lastresult)
        # return resp

        #return multi-job report
        if "multi" in routines.getjobtype(jobid):
            return render_template("showmultiresult.html",jobid=jobid)
        return render_template("showresult.html",jobid=jobid)
    return abort(404)

@app.route('/results/<jobid>/tables/<path:path>')
def gettables(jobid,path):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid,"tables")
    if os.path.exists(os.path.join(rd,path)):
        return send_from_directory(rd,path)
    return abort(404)

@app.route('/results/<jobid>/xlfile')
def getxlfile(jobid):
    if jobid and routines.checkresult(jobid):
        tpath = os.path.join(app.config["RESULTS_FOLDER"],jobid,"tables")
        xlfile = routines.makexltable(tpath)
        if xlfile:
            return send_from_directory(tpath,xlfile)
        return abort(404)

@app.route('/results/<jobid>/status')
def jobstatus(jobid):
    if jobid and routines.checkresult(jobid):
        #rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
        # statuscache = os.path.join(rd,"jobstatus.json")
        # if not os.path.exists(statuscache) or os.path.getmtime(statuscache) < os.path.getmtime(os.path.join(rd,"arts-query.log")):
        status = routines.getjobstatus(jobid)
        return jsonify(status)
        # else:
        #     return send_from_directory(rd,"jobstatus.json")
    return jsonify({"data":[]})

@app.route('/results/<jobid>/multistatus')
def multijobstatus(jobid):
    if jobid and routines.checkresult(jobid):
        #rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
        # statuscache = os.path.join(rd,"jobstatus.json")
        # if not os.path.exists(statuscache) or os.path.getmtime(statuscache) < os.path.getmtime(os.path.join(rd,"arts-query.log")):
        status = routines.getmultistatus(jobid)
        return jsonify(status)
        # else:
        #     return send_from_directory(rd,"jobstatus.json")
    return jsonify({"data":[]})

@app.route('/results/<jobid>/antismash')
@app.route('/results/<jobid>/antismash/')
def gotoantismashfiles(jobid):
    return redirect("/results/%s/antismash/index.html"%jobid)

@app.route('/results/<jobid>/antismash/<path:path>')
def getantismash(jobid,path):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid,"antismash")
    if os.path.exists(os.path.join(rd,path)):
        return send_from_directory(rd,path)
    return abort(404)

@app.route('/results/<jobid>/bigscape')
@app.route('/results/<jobid>/bigscape/')
def gotobigscapefiles(jobid):
    return redirect("/results/%s/bigscape/index.html"%jobid)

@app.route('/results/<jobid>/bigscape/<path:path>')
def getbigscape(jobid,path):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid,"antismash_bigscape_result")
    if os.path.exists(os.path.join(rd,path)):
        return send_from_directory(rd,path)
    return abort(404)

@app.route('/results/<jobid>/log')
def joblog(jobid):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"arts-query.log")):
        return send_from_directory(rd,"arts-query.log",mimetype="text/plain")
    return "No log found"

@app.route('/results/<jobid>/export/<expfil>')
def exportfile(jobid,expfil):
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,expfil)):
        return send_from_directory(rd,expfil)
    return abort(404)

@app.route('/archive/<jobfile>')
def getarchive(jobfile):
    jobid = os.path.splitext(jobfile)[0]
    rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
    ad = os.path.join(app.config.get("ARCHIVE_FOLDER","/tmp"))
    if os.path.exists(os.path.join(ad,jobfile)):
        return send_from_directory(ad,jobfile)
    elif jobid and routines.checkresult(jobid):
        shutil.make_archive(os.path.join(ad,str(jobid)),"zip",rd)
        return send_from_directory(ad,jobfile)
    return abort(404)

# @app.route('/results/<jobid>/<action>')
# def api_action(jobid,action):
#     return redirect("/api/job/%s/%s"%(jobid,action))

@app.route('/results/<jobid>/genometable')
def genometab(jobid):
    rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"tables","summary_table.tsv")):
        rows = []
        with open(os.path.join(rd,"tables","summary_table.tsv"),"r") as fil:
            for line in fil:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                rows.append(cols)
        return jsonify({"data":rows})
    elif jobid and routines.checkresult(jobid):
        #otherwise list the finished job folders without stats
        startedjobs = [[x]+[x]+["-"]*9 for x in os.listdir(rd) if jobid in x]
        return jsonify({"data":startedjobs})
    else:
        return jsonify({"data":[]})

@app.route('/results/<jobid>/multicoretab')
def multicoretab(jobid):
    rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"tables","combined_core_table.tsv")):
        rows = []
        with open(os.path.join(rd,"tables","combined_core_table.tsv"),"r") as fil:
            for line in fil:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                IDs = {}
                for ind,desc in enumerate(cols[7].split(";")):
                    f = desc.split("|")
                    title = f[0]
                    if title not in IDs:
                        IDs[title] = ([],ind)
                    IDs[title][0].append(f[-3])
                # #Replace description with gene ID in fasta header
                # regx = 'ID=[0-9]*_[0-9]*'
                # IDs = re.findall(regx,cols[7])
                #Remove full paths
                cols = cols[:8] + [[str(os.path.split(x)[-1]) for x in json.loads(col.replace("'",'"'))] for col in cols[8:]]
                k,l = zip(*[(x[1:],y[0]) for x,y in sorted(IDs.items(), key=lambda t: t[1][1])])
                cols[7] = list(k)
                cols.append([["Loc: "+y.replace(" # ID="," @Contig-seq: ") for y in x] for x in list(l)])
                rows.append(cols)
        return jsonify({"data":rows})
    else:
        return jsonify({"data":[]})

@app.route('/results/<jobid>/multikrtab')
def multikrtab(jobid):
    rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"tables","combined_known_table.tsv")):
        rows = []
        with open(os.path.join(rd,"tables","combined_known_table.tsv"),"r") as fil:
            for line in fil:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")[:3]
                rows.append(cols)
        return jsonify({"data":rows})
    else:
        return jsonify({"data":[]})

@app.route('/results/<jobid>/multiduptab')
def multiduptab(jobid):
    rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"tables","combined_dup_table.tsv")):
        rows = []
        with open(os.path.join(rd,"tables","combined_dup_table.tsv"),"r") as fil:
            for line in fil:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")[:3]
                rows.append(cols)
        return jsonify({"data":rows})
    else:
        return jsonify({"data":[]})

@app.route('/results/<jobid>/multibgctab')
def multibgctab(jobid):
    rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    bgctable = False
    if os.path.exists(os.path.join(rd,"tables","new_combined_bgc_table.tsv")):
        bgctable = os.path.join(rd,"tables","new_combined_bgc_table.tsv")
    elif os.path.exists(os.path.join(rd,"tables","combined_bgc_table.tsv")):
        bgctable = os.path.join(rd,"tables","combined_bgc_table.tsv")
    if jobid and routines.checkresult(jobid) and bgctable:
        rows = []
        # hitdict = {}
        with open(bgctable,"r") as fil:
            for line in fil:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                # Fix shared counts, (remove single family member shared count)
                if "1" == cols[1]:
                    cols[7] = "0"
                    cols[8] = "0"
                #Fix string for proper JSON
                cols = [json.loads(x.replace("u'",'"').replace("'",'"')) for x in cols]
                rows.append(cols)
        # return jsonify({"data":rows,"hits":hitdict})
        return jsonify({"data":rows})
    else:
        return jsonify({"data":[]})

#
# @app.route('/results/<jobid>/bgctable')
# def bgctable(jobid):
#     rd = os.path.join(app.config['RESULTS_FOLDER'],jobid)
#     if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"resultsummary.json")):
#         return jsonify(routines.getbgctable(rd))
#     else:
#         return jsonify({"data":[]})
#
# @app.route('/results/<jobid>/summarytab')
# def summarytab(jobid):
#     rd = os.path.join(app.config["RESULTS_FOLDER"],jobid)
#     if jobid and routines.checkresult(jobid) and os.path.exists(os.path.join(rd,"arts-query.log")):
#         st = routines.formatsumtab(jobid)
#         return jsonify(st)
#     return jsonify({"data":[]})

@app.route('/download')
def download():
    return render_template("download.html")

@app.route('/example')
def example():
    return redirect("/results/example/report")

@app.route('/about')
def about():
    return render_template("about.html")

@app.route('/help')
def help():
    return render_template("help.html")

@app.route('/startjob', methods=['GET', 'POST'])
def startjob():
    if request.method == 'POST':
        multifiles = request.form.getlist("upfiles")
        if multifiles:
            filename = ",".join(multifiles)
            asrun = True
        else:
            filename,asrun = routines.getinfile()
        custmdl = routines.getcustmdl("custhmm")
        custcoremdl = routines.getcustmdl("custcorehmm")
        options = routines.getoptions()
        searchmode = request.form.get("searchmode","")
        if filename:
            if custmdl and custmdl=="Failed":
                redirect("/analyze")
            if custcoremdl and custcoremdl=="Failed":
                redirect("/analyze")
            artsjob = routines.addjob(infile=filename,ref=request.form.get("refset","auto"),asjob=request.form.get('asjobid',""),asrun=asrun,custmdl=custmdl,custcoremdl=custcoremdl,options=options,cut=searchmode)

            resp = make_response(redirect("/results/%s/report"%artsjob.id))
            saverun = request.form.get('saverun',False)
            if saverun and saverun.lower() != "false":
                lastresult = request.cookies.get('arts.lastresult')
                if lastresult and artsjob.id not in lastresult:
                    lastresult += ";%s"%artsjob.id
                else:
                    lastresult = artsjob.id
                resp.set_cookie('arts.lastresult',lastresult)
            email = request.form.get('email',False)
            if email:
                # Run threaded so mail server resp does not block process
                @copy_current_request_context
                def sendmail(x,y,z):
                    routines.sendnotifymail(x,y,z)

                mailer = threading.Thread(name='mail_sender', target=sendmail, args=("",artsjob.id,email))
                mailer.start()
            validchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 _-()"
            jobtitle = request.form.get('jobtitle',False)
            jobtitle = ''.join([c for c in str(jobtitle[:30]) if c in validchars])
            if jobtitle and len(jobtitle):
                with open(os.path.join(app.config['RESULTS_FOLDER'],artsjob.id,"jobtitle.txt"),"w") as fil:
                    fil.write(jobtitle+"\n")
            return resp

            # return redirect("/results/%s"%artsjob.id)
        elif asrun == "invalid":
            flash("No such Antismash job found, please check Antismash id and try again")
            return redirect("/analyze")
        elif asrun == "pending":
            flash("Antismash run: %s is pending, please try again when results are available"%request.form.get("asjobid",""))
            return redirect("/analyze")
        else:
            flash("Error getting sequence file, please ensure file is in correct format and try again")
            return redirect("/analyze")
    else:
        flash("Please fill in required fields to start a job")
    return redirect("/analyze")

@app.errorhandler(404)
@app.errorhandler(401)
@app.errorhandler(500)
def page_not_found(e):
    return render_template('error.html',title="",errormsg=e)

