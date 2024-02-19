# ARTS Web Server Overview

This is a sub repository for the Antibiotic Resistant Target Seeker (ARTS).

This can be used to view results generated from the public server at https://arts.ziemertlab.com
or using output from the main analysis pipeline at https://github.com/ziemertlab/arts

For usage of the web server see https://arts.ziemertlab.com/help


# Installation of ARTS Web Server

There are three options for installing ARTS Web Server:

- Using Docker Images 
- Using Anaconda/Miniconda
- Manual Installation for Linux/Ubuntu

## 1- Using Docker Image:

- Firstly, if you don't have Docker, you should install the Docker engine on your computer. Please check out the latest version of Docker 
[on the official website.](https://docs.docker.com/get-docker/)

- Edit desired paths to run Docker Images:
````
I- Desired Result Directory: /my/path/to/result_folders:/results/
II- Desired Log Directory: /my/path/to/log_file:/run/
````
- Enter the required paths and run the command:

````bash
docker run -it -v /my/path/to/results_folders:/results/ -v /my/path/to/log_file:/run/ -p 5000:5000 --rm ziemertlab/artswebapp-beta:latest
````

- You can see the other details [on Docker Hub](https://hub.docker.com/r/ziemertlab/artswebapp-beta)

## 2- Using Anaconda/Miniconda:
We recommend [Anaconda3/Miniconda3](https://docs.anaconda.com/free/anaconda/install/index.html) (with python >=3.8) and 
it is necessery for the [conda](https://docs.conda.io/en/latest/index.html) package manager.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/artswebapp
```
- Enter the artswebapp folder:
```bash
    cd artswebapp
```
- Create a new environment and install all the packages using the environment.yml file with conda:
```bash
    conda env create -f environment.yml
```
- Activate artswebapp environment:
```bash
    conda activate artswebapp
```
- Edit desired folders in configs (config/artsapp_default.conf and config/uwsgi.conf)
  (See [Confugiration](https://github.com/ZiemertLab/ARTSwebapp/tree/master#Configuration-of-ARTS-Web-Server) for more):
- Run server (from artswebapp folder)
  (See [Usage](https://github.com/ziemertlab/arts#Usage_of_ARTS_Web_Server) for more):
```bash
    uwsgi --ini config/uwsgi.conf
```


## 3- Manual Installation for Linux/Ubuntu:

**Note:** Python version 3.8 or higher is recommended.

- Clone/Download the repository (root / sudo required):
```bash
    git clone https://github.com/ziemertlab/artswebapp
```
- Enter the artswebapp folder:
```bash
    cd artswebapp
```
- Install required libraries and applications (root / sudo required):
```bash
    pip install -r requirements.txt
```
- Edit desired folders in configs (config/artsapp_default.conf and config/uwsgi.conf)
  (See [Confugiration](https://github.com/ZiemertLab/ARTSwebapp/tree/master#Configuration-of-ARTS-Web-Server) for more):
- Run server (from artswebapp folder)
  (See [Usage](https://github.com/ZiemertLab/ARTSwebapp/tree/master#Usage-of-ARTS-Web-Server) for more):
```bash
    uwsgi --ini config/uwsgi.conf
```


# Configuration of ARTS Web Server
- Edit desired folders in configs (config/artsapp_default.conf and config/uwsgi.conf) and write your working directories instead of "~PATH":
````
EXAMPLE:
config/artsapp_default.conf:
        ...
        UPLOAD_FOLDER = "~PATH/uploads"
        RESULTS_FOLDER = "~PATH/results"
        ARCHIVE_FOLDER = "~PATH/archive"
        ...
config/uwsgi.conf:
        ...
        logto = ~PATH/artswebapp.log
        stats = ~PATH/uwsgi.stats.sock
        touch-reload = ~PATH/uwsgi.reload
        pidfile = ~PATH/uwsgi.pid
        ...
````

# Usage of ARTS Web Server

- Run server (from artswebapp folder):
```bash
    uwsgi --ini config/uwsgi.conf
```
**Note:** It may need to run "redis-server" on the terminal.

- Click your local server on your browser:

```
    http://127.0.0.1:5000/
```

**Note:** The link may differ according to your configuration.

- Click result page to view ARTS results:
```
    http://127.0.0.1:5000/results
```
- Enter the result file name and click "View Report":

**Note:** Make sure your result folders are in the specified results folder path of artswebapp.conf.
````
config/artsapp_default.conf:
        ...
        RESULTS_FOLDER = "~PATH/results"
        ...
````
- You can view your ARTS results now!

### Optional:
If you do not have a ARTS result file yet, you can download the sample result from the link below (~55MB). 
After extracting the zip file, follow the relevant steps.

````bash
  wget https://arts.ziemertlab.com/archive/GCF_000015125.1.zip
  unzip ~PATH/GCF_000015125.1.zip -d ~PATH/results/GCF_000015125.1
````

# Optional - Submission a job using local webserver:
To start the analaysis on local webserver, please see https://github.com/ziemertlab/arts and install ARTS.

- Edit desired folders in configs (config/artsapp_default.conf and config/uwsgi.conf) and write all your working directories instead of "~PATH"

- Then, run server (from artswebapp folder)::
```bash
    uwsgi --ini config/uwsgi.conf
```
**Note:** It may need to run "redis-server" on the terminal.

- To submit an input file, run "runjobs.py" on the terminal:
```bash
    cd arts
    python runjobs.py run -pid /tmp/runjobs.pid
```
- Local webserver is ready to analyse!!



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