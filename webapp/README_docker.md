# ARTS Web Server Overview

This is a sub repository for the Antibiotic Resistant Target Seeker (ARTS).

This can be used to view results generated from the public server at https://arts.ziemertlab.com
or using output from the main analysis pipeline at https://bitbucket.org/ziemertlab/arts

For usage of the web server see https://arts.ziemertlab.com/help


# Quickstart with Docker:

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

- Click your local server on your browser:

```
    http://127.0.0.1:5000/
```

- Click result page to view ARTS results:
```
    http://127.0.0.1:5000/results
```
- Enter the result file name and click "View Report":

**Note:** Make sure your result folders are in the specified results folder path in the command.

- You can view your ARTS results now!

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