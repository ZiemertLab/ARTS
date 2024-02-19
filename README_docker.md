# Antibiotic Resistant Target Seeker (ARTS) Overview

ARTS is a webserver and analysis pipeline for screening for known and
putative antibiotic resistance markers in order to identify and prioritize
their corresponding biosynthetic gene clusters. ARTS can be installed locally
or you can use the free public webserver located at https://arts.ziemertlab.com

See https://bitbucket.org/ziemertlab/arts for a guide on installing the main analaysis pipeline independantly.

See https://bitbucket.org/ziemertlab/artswebapp for a guide on installing the webserver independantly.

## Quick start with docker

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- To run ARTS Image, you should download the "docker_run_arts.py" file from the command line or from the repository using a web browser.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://arts.ziemertlab.com/static/docker_run/docker_run_arts.py
```
**Note:** Python 3.x is needed to run "docker_run_arts.py".

- ARTS Image include only actinobacteria reference set. If you need other reference sets, please download and unzip all of them.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Enter the required arguments and run the script
```bash
    python docker_run_arts.py [-h] [input] [resultdir] [-optional_arguments]
```

## Arguments of docker_run_arts.py
````
usage: docker_run_arts.py [-h] [-refdir REFDIR] [-hmms HMMDBLIST] [-cchmms CUSTCOREHMMS] [-chmms CUSTOMHMMS]
                          [-rhmm RNAHMMDB] [-t THRESH] [-cpu MULTICPU] [-opt OPTIONS] [-org ORGNAME]
                          [-pbt PREBUILTTREES] [-ras] [-rbsc]
                          input resultdir

Start from genbank file and compare with pre-computed reference for Duplication and Transfers

positional arguments:
  input                 gbk file to start query. For multi-genome input, put commas without any space between the
                        paths.
  resultdir             Directory to store results

optional arguments:
  -h, --help            show this help message and exit
  -refdir REFDIR        Directory of precomputed reference files (default="~/actinobacteria/")
  -hmms HMMDBLIST, --hmmdblist HMMDBLIST
                        hmm file, directory, or list of hmm models for core gene id
  -cchmms CUSTCOREHMMS, --custcorehmms CUSTCOREHMMS
                        User supplied core models. hmm file
  -chmms CUSTOMHMMS, --customhmms CUSTOMHMMS
                        User supplied resistance models. hmm file
  -rhmm RNAHMMDB, --rnahmmdb RNAHMMDB
                        RNA hmm models to run (default: None)
  -t THRESH, --thresh THRESH
                        Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA,
                        trusted= TC, noise= NC(default: none)
  -cpu MULTICPU, --multicpu MULTICPU
                        Turn on Multi processing set # Cpus (default: Off, 1)
  -opt OPTIONS, --options OPTIONS
                        Analysis to run. phyl=phylogeny, kres=known resistance, duf=Domain of unknown function,
                        expert=Exploration mode, (default: phyl,kres,duf)
  -org ORGNAME, --orgname ORGNAME
                        Explicitly specify organism name
  -pbt PREBUILTTREES, --prebuilttrees PREBUILTTREES
                        Directory of prebuilt trees
  -ras, --runantismash  Run input file through antismash first
  -rbsc, --runbigscape  Run antismash results through bigscape
````
## Usage 

- For basic run with positional arguments;
````
    python docker_run_arts.py /PATH/input_genome.gbk /PATH/result_folder
````

- To use antiSMASH, please use: `-ras`, `--runantismash`;
````
    python docker_run_arts.py /PATH/input_genome.gbk /PATH/result_folder -ras
````
- To use other reference sets (default: actinobacteria);
````
    python docker_run_arts.py /PATH/input_genome.gbk /PATH/result_folder -refdir /PATH/reference_set_file
````
- For multi-genome input, it is enough to put commas without any space between the paths of genome files;
````
    python docker_run_arts.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk,/PATH/input_genome3.gbk /PATH/result_folder
````

- It is possible to run ARTS if there is a .json files of antiSMASH 6 results;
````
    python docker_run_arts.py /PATH/antismash_result.json /PATH/result_folder
````

- To run ARTS analysis options, please use `-opt`, `--options` parameter (default: phyl,kres,duf);
````
    python docker_run_arts.py /PATH/input_genome.gbk /PATH/result_folder -ras -opt kres,duf,expert 
````

- To run the BiG-SCAPE algorithms, please use `-rbsc`, `--runbigscape`;
````
    python docker_run_arts.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk /PATH/result_folder -ras -rbsc
````
# Support
If you have any issues please feel free to contact us at arts-support@ziemertlab.com

# Licence
This software is licenced under the GPLv3. See LICENCE.txt for details.

# Publication
If you found ARTS to be helpful, please [cite us](https://doi.org/10.1093/nar/gkaa374):

Mungan,M.D., Alanjary,M., Blin,K., Weber,T., Medema,M.H. and Ziemert,N. (2020) ARTS 2.0: feature updates and expansion 
of the Antibiotic Resistant Target Seeker for comparative genome mining. 
Nucleic Acids Res.,[10.1093/nar/gkaa374](https://doi.org/10.1093/nar/gkaa374)

Alanjary,M., Kronmiller,B., Adamek,M., Blin,K., Weber,T., Huson,D., Philmus,B. and Ziemert,N. (2017) The Antibiotic 
Resistant Target Seeker (ARTS), an exploration engine for antibiotic cluster prioritization and novel drug target discovery. 
Nucleic Acids Res.,[10.1093/nar/gkx360](https://doi.org/10.1093/nar/gkx360)