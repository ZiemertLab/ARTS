# Antibiotic Resistant Target Seeker (ARTS) Overview

ARTS is a webserver and analysis pipeline for screening for known and putative antibiotic resistance markers in order to identify and prioritize their corresponding biosynthetic gene clusters. 
ARTS allows for specific and efficient genome mining for antibiotics with interesting and novel targets by rapidly linking housekeeping and known resistance genes to BGC proximity, duplication and horizontal gene transfer (HGT) events.

ARTS can be installed locally, or you can use the free public webserver located at https://arts.ziemertlab.com

See https://github.com/ziemertlab/artswebapp for a guide on installing the webserver independently.


# Installation of ARTS

There are three options for installing ARTS:

- Using Docker Images 
- Using Anaconda/Miniconda
- Manual Installation for Linux/Ubuntu

## 1- Using Docker Image:

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- To run ARTS Image, you should download the "docker_run_arts.py" file from the command line or from the repository using a web browser.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://github.com/ziemertlab/arts/raw/master/docker_run_arts.py
```
**Note:** Python 3.x is needed to run "docker_run_arts.py".

- ARTS Image include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    mkdir ARTSdocker && cd ARTSdocker
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Enter the required arguments and run the script
```bash
    python docker_run_arts.py [-h] [input] [resultdir] [-optional_arguments]
```
- You can see the other details [on Docker Hub](https://hub.docker.com/r/ziemertlab/arts)

## 2- Using Anaconda/Miniconda:
We recommend [Anaconda3/Miniconda3](https://docs.anaconda.com/free/anaconda/install/index.html) (with python >=3.8) and 
it is necessery for the [conda](https://docs.conda.io/en/latest/index.html) package manager.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/arts
```
- Enter the arts folder:
```bash
    cd arts
```
- ARTS GitHub repository include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Create a new environment and install all the packages using the environment.yml file with conda:
```bash
    conda env create -f environment.yml
```
- Activate arts environment and run ARTS (See [Usage](https://github.com/ZiemertLab/ARTS/tree/master#usage) for more):
```bash
    conda activate arts
```
```bash
    python artspipeline1.py [-h] [input] [refdir] [-optional_arguments]
```

## 3- Manual Installation for Linux/Ubuntu:
The analysis server will start a local antiSMASH job if cluster annotation is not already provided as input. We recommend antiSMASH version >= 6.0.1.
See [antiSMASH](https://docs.antismash.secondarymetabolites.org/install/) for installation instructions.

**Note:** Python version 3.8 or higher is recommended.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/arts
```
- Enter the arts folder:
```bash
    cd arts
```
- Prepare Actinobacteria reference set and related HMMs:
```bash
    unzip reference/'*.zip' -d reference/ 
```
- ARTS GitHub repository include only Actinobacteria reference set. If you need other reference sets, please download (~2.2GB) and unzip all of them.
```bash
    wget https://arts.ziemertlab.com/static/zip_refsets/all_references.zip
    unzip all_references.zip 
```
- Install required libraries and applications (root / sudo required):
```bash
    apt-get update
    apt-get install -y python3-dev liblzma-dev default-jdk hmmer2 hmmer diamond-aligner fasttree prodigal ncbi-blast+ muscle mafft
    pip install -r requirements.txt
```
- Install required binaries or use pre-compiled linux64bit bins (root / sudo required):
  - Dependencies:
    - TrimAl : trimal => https://github.com/inab/trimal
    - RaxML : raxmlHPC-SSE3 => https://github.com/stamatak/standard-RAxML
    - Ranger-DTL : ranger-dtl-U => http://compbio.mit.edu/ranger-dtl/
    - Glimmer : glimmer3 => https://ccb.jhu.edu/software/glimmer/index.shtml
    - GlimmerHMM : glimmerhmm => https://ccb.jhu.edu/software/glimmerhmm/
  - Pre-compiled bins:
    ```bash
        tar -zxvf linux_64bins.tar.gz -C /usr/local/bin/
    ```
- 
- Run ARTS (See [Usage](https://github.com/ZiemertLab/ARTS/tree/master#usage) for more):
```bash
    python artspipeline1.py [-h] [input] [refdir] [-optional_arguments]
```


## Optional: For comparing the results of multi-genome analysis:

The BiG-SCAPE algorithm is used to compare the results of multi-genome analysis. 
All clustered BGCs from antiSMASH results are analyzed to determine BGC similarity. 
The BiG-SCAPE algorithm generates sequence similarity networks of BGCs and classifies them into gene cluster families (GCFs).

To install [the BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org/index.html), please see https://github.com/medema-group/BiG-SCAPE/wiki/installation

**Note:** Make sure that the Pfam database is in the same folder as bigscape.py

# Running ARTS
ARTS uses a webserver to queue jobs to the analysis pipeline. Details on webserver usage can be found at: https://arts.ziemertlab.com/help 

Alternatively jobs can be run directly using the artspipeline1.py script (see -h for options).

````
usage: artspipeline1.py [-h] [-hmms HMMDBLIST] [-khmms KNOWNHMMS] [-duf DUFHMMS] [-cchmms CUSTCOREHMMS] [-chmms CUSTOMHMMS] [-rhmm RNAHMMDB] [-t THRESH]
                        [-td TEMPDIR] [-rd RESULTDIR] [-ast ASTRAL] [-cpu MULTICPU] [-opt OPTIONS] [-org ORGNAME] [-pbt PREBUILTTREES] [-ras]
                        [-asp ANTISMASHPATH] [-bcp BIGSCAPEPATH] [-rbsc]
                        input refdir

Start from genbank file and compare with pre-computed reference for Duplication and Transfers

positional arguments:
  input                 gbk file to start query
  refdir                Directory of precomputed reference files

optional arguments:
  -h, --help            show this help message and exit
  -hmms HMMDBLIST, --hmmdblist HMMDBLIST
                        hmm file, directory, or list of hmm models for core gene id
  -khmms KNOWNHMMS, --knownhmms KNOWNHMMS
                        Resistance models hmm file
  -duf DUFHMMS, --dufhmms DUFHMMS
                        Domains of unknown function hmm file
  -cchmms CUSTCOREHMMS, --custcorehmms CUSTCOREHMMS
                        User supplied core models. hmm file
  -chmms CUSTOMHMMS, --customhmms CUSTOMHMMS
                        User supplied resistance models. hmm file
  -rhmm RNAHMMDB, --rnahmmdb RNAHMMDB
                        RNA hmm models to run (default: None)
  -t THRESH, --thresh THRESH
                        Hmm reporting threshold. Use global bitscore value or Model specific options: gathering= GA, trusted= TC, noise= NC(default: none)
  -td TEMPDIR, --tempdir TEMPDIR
                        Directory to create unique results folder
  -rd RESULTDIR, --resultdir RESULTDIR
                        Directory to store results
  -ast ASTRAL, --astral ASTRAL
                        Location of Astral jar executable default: Value of environment var 'ASTRALJAR'
  -cpu MULTICPU, --multicpu MULTICPU
                        Turn on Multi processing set # Cpus (default: Off, 1)
  -opt OPTIONS, --options OPTIONS
                        Analysis to run. phyl=phylogeny, kres=known resistance, duf=Domain of unknown function, expert=Exploration mode (default: phyl,kres,duf)
  -org ORGNAME, --orgname ORGNAME
                        Explicitly specify organism name
  -pbt PREBUILTTREES, --prebuilttrees PREBUILTTREES
                        Directory of prebuilt trees
  -ras, --runantismash  Run input file through antismash first
  -asp ANTISMASHPATH, --antismashpath ANTISMASHPATH
                        Location of the executable file of antismash or location of antismash 'run_antismash.py' script
  -bcp BIGSCAPEPATH, --bigscapepath BIGSCAPEPATH
                        location of bigscape 'bigscape.py' script
  -rbsc, --runbigscape  Run antismash results through bigscape
````

# Usage 

- For basic run with positional arguments;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria
````

- To save all output data files: `-rd`, `--resultdir`
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -rd /PATH/result_folder
````

- To use antiSMASH: `-asp`, `--antismashpath` and to run antiSMASH: `-ras`, `--runantismash`
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -rd /PATH/result_folder
````

- If there is an exsiting antiSMASH job, .json files of antiSMASH results are available fo ARTS: `-asp`, `--antismashpath`
````
    python artspipeline1.py /PATH/antismash_result.json /PATH/arts/reference/actinobacteria -asp /PATH/antismash -rd /PATH/result_folder
````

- To run ARTS with exploration mode, please use `-opt`, `--options` parameter;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -opt 'expert' 
````

- To identify known resistance, please use `-khmms`, `--knownhmms` and `-opt`, `--options` parameters;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -khmms /PATH/arts/reference/knownresistance.hmm -opt 'kres'
````

- To identify domain of unknown function(DUF), please use `-duf`, `--dufhmms` and `-opt`, `--options` parameters;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -khmms /PATH/arts/reference/dufmodels.hmm -opt 'duf'
````

- To run ARTS with phylogeny screening, please use `-ast`, `--astral` and `-opt`, `--options` parameter ;
````
    python artspipeline1.py /PATH/input_genome.gbk /PATH/arts/reference/actinobacteria -asp /PATH/antismash -ras -ast /PATH/arts/astral/astral.5.7.7.jar -opt 'phly' 
````

- For multi-genome input, it is enough to put commas without any space between the paths of genome files;
````
    python artspipeline1.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk,/PATH/input_genome3.gbk /PATH/arts/reference/actinobacteria -rd /PATH/result_folder
````

- To run the BiG-SCAPE algorithms, please use `-bcp`, `--bigscapepath` and `-rbsc`, `--runbigscape`
````
    python artspipeline1.py /PATH/input_genome1.gbk,/PATH/input_genome2.gbk /PATH/arts/reference/actinobacteria -bcp /PATH/BiG-SCAPE_1.1.5/bigscape.py -rbsc -rd /PATH/result_folder
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