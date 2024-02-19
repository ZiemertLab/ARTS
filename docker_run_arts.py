#!/usr/bin/env python
# Copyright (C) 2023 Turgut Mesut YILMAZ
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

import argparse, os

def call_startquery(args):
    reference= "/arts/reference/actinobacteria/"
    container_volume_list = []
    container_argument_list = []
    docker_run_cmd = "docker run -it --rm "
###Input volume - path in container

    input_paths = []
    if "," in args.input:
        input_list = args.input.split(",")
        for i in input_list:
            path,file = os.path.split(i)
            input_paths.append("/arts/uploads/%s" % file)
            container_volume_list.append("-v %s:/arts/uploads/" % path)
    else:
        path, file = os.path.split(args.input)
        input_paths.append("/arts/uploads/%s" % file)
        container_volume_list.append("-v %s:/arts/uploads/" % path)


 ##Reference volume - path
    if args.refdir != reference:
        path, file = os.path.split(args.refdir)
        container_ref_arg = "/custom_reference/%s" % file
        container_volume_list.append("-v %s:/custom_reference/" % path)

    else:
        container_ref_arg = reference


### Directory to store results
    path, file = os.path.split(args.resultdir)
    container_volume_list.append("-v %s:/arts/results/" % path)
    container_argument_list.append("-rd /arts/results/%s" % file)

### Antismash path
    if args.runantismash is True:
        container_argument_list.append("-ras -asp /usr/local/bin/antismash")
    else:
        container_argument_list.append("-asp /usr/local/bin/antismash")

### Custom hmm models list
    if args.hmmdblist is not None:
        path, file = os.path.split(args.hmmdblist)
        container_argument_list.append("-hmms /arts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/arts/custom_files/" % path)

### Known Resistance models hmms
    if "kres" in args.options:
        container_argument_list.append("-khmms /arts/reference/knownresistance.hmm")

### Domains of unknown function hmm file
    if "duf" in args.options:
        container_argument_list.append("-duf /arts/reference/dufmodels.hmm")

### User supplied core models
    if args.custcorehmms is True:
        path, file = os.path.split(args.custcorehmms)
        container_argument_list.append("-cchmms /arts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/arts/custom_files/" % path)

### User supplied resistance models
    if args.customhmms is True:
        path, file = os.path.split(args.customhmms)
        container_argument_list.append("-chmms /arts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/arts/custom_files/" % path)

### RNA hmm models to run
    if args.rnahmmdb is True:
        path, file = os.path.split(args.rnahmmdb)
        container_argument_list.append("-rhmm /arts/custom_files/%s"%file)
        container_volume_list.append("-v %s:/arts/custom_files/" % path)

### Hmm reporting threshold
    if args.thresh is not None:
        container_argument_list.append("-t %s"%args.thresh)

### Multi processing set - CPU
    if args.multicpu > 1:
        container_argument_list.append("-cpu %s" % args.multicpu)
        docker_run_cmd += "--cpus=%s"%args.multicpu + " "

### Opions for Analysis to run
    container_argument_list.append("-opt %s"%args.options)

### Astral jar path
    container_argument_list.append("-ast /arts/astral/astral.5.7.7.jar")

### Explicitly specify organism name
    if args.orgname is not None:
        container_argument_list.append("-org %s"%args.orgname)

### Directory of prebuilt trees
    if args.prebuilttrees is not None:
        path, file = os.path.split(args.prebuilttrees)
        container_argument_list.append("-pbt /arts/custom_files/%s" % file)
        container_volume_list.append("-v %s:/arts/custom_files/" % path)

### Run BiG-Scape Path:
    if args.runbigscape is True:
        container_argument_list.append("-rbsc -bcp /arts/BiG-SCAPE/bigscape.py")


## Docker Run ARTS
    docker_run_cmd += " ".join(list(set(container_volume_list)))
    docker_run_cmd += " ziemertlab/arts-beta:3.0.b2 "
    docker_run_cmd += ",".join(input_paths) + " " + container_ref_arg + " "
    docker_run_cmd += " ".join(container_argument_list)
    print(docker_run_cmd)
    os.system(docker_run_cmd)


# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Start from genbank file and compare with pre-computed reference for Duplication and Transfers""")
    parser.add_argument("input", help="gbk file to start query. For multi-genome input, put commas without any space between the paths.")
    parser.add_argument("resultdir", help="Directory to store results")
    parser.add_argument("-refdir", help="Directory of precomputed reference files (default='~/actinobacteria/')",default="/arts/reference/actinobacteria/")
    parser.add_argument("-hmms","--hmmdblist", help="hmm file, directory, or list of hmm models for core gene id",default=None)
    #parser.add_argument("-khmms","--knownhmms", help="Resistance models hmm file",default="/arts/reference/knownresistance.hmm")
    #parser.add_argument("-duf","--dufhmms", help="Domains of unknown function hmm file",default="/arts/reference/dufmodels.hmm")
    parser.add_argument("-cchmms","--custcorehmms", help="User supplied core models. hmm file",default=False)
    parser.add_argument("-chmms","--customhmms", help="User supplied resistance models. hmm file",default=False)
    parser.add_argument("-rhmm","--rnahmmdb", help="RNA hmm models to run (default: None)",default=None)
    parser.add_argument("-t","--thresh", help="Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted= TC, noise= NC(default: none)",default=None)
    # parser.add_argument("-td", "--tempdir", help="Directory to create unique results folder", default=None)
    #parser.add_argument("-ast", "--astral", help="Location of Astral jar executable default: Value of environment var 'ASTRALJAR' ", default=None)
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing set # Cpus (default: Off, 1)", type=int, default=1)
    parser.add_argument("-opt", "--options", help="Analysis to run. phyl=phylogeny, kres=known resistance, duf=Domain of unknown function, expert=Exploration mode, (default: phyl,kres,duf)", default="phyl,kres,duf")
    parser.add_argument("-org", "--orgname", help="Explicitly specify organism name", default=None)
    parser.add_argument("-pbt", "--prebuilttrees", help="Directory of prebuilt trees", default=None)
    parser.add_argument("-ras", "--runantismash", help="Run input file through antismash first", action='store_true', default=False)
    #parser.add_argument("-asp", "--antismashpath", help="location of antismash 'run_antismash.py' script", default=False)
    #parser.add_argument("-bcp", "--bigscapepath", help="location of bigscape 'bigscape.py' script", default=False)
    parser.add_argument("-rbsc", "--runbigscape",help="Run antismash results through bigscape", action='store_true', default=False )
    args = parser.parse_args()
    call_startquery(args)