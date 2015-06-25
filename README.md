# Transcriptome compression

This repository is composed of several directories. Each is listed below with a brief description of its purpose and contents. Each directory has its own README that explains its contents in greater detail.

# Repository directories
1. data_download: Contains the custom scripts necessary to download and quantify raw RNA-Seq reads from NCBI's Short Read Archive (SRA).
2. util: Contains accessory code.
  
# Aspera installation
1. Select the appropriate download of aspera connect from http://downloads.asperasoft.com/en/downloads/8?list
2. For linux and mac, this should be a shell script. Make it executable `chmod +x aspera-connect-3.5.1.92523-linux-64.sh`
3. Run it. It should install in your home directory. `./aspera-connect-3.5.1.92523-linux-64.sh`

# Set up tradict  
## Add tradict directory to your $PATH environment variable  
1. Navigate to the tradict copy on your machine with `cd path/to/tradict`   
2. Display the absolute path to tradict with `pwd`. Copy this path into your memory buffer.
3. Open you bashrc profile. `nano ~/.bashrc` for linux users. `sudo nano /etc/profile` for mac users.
4. Add the line `export PATH=/path/to/tradict:$PATH` to your .bashrc profile, where `/path/to/tradict` is the path that was copied into the memory buffer in step 2.
5. Save and exit with Ctrl+x.
6. Type `source ~/.bashrc` as linux user or `source /etc/profile` as mac user to load the modified profile.

# Misc. notes
1. If getting the error: "name not found while resolving tree within virtual file system module" Try updating sra-toolkit to latest release.

