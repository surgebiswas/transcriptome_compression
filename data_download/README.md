This document describes how to use `srafish.pl` to download and quantify raw RNA-Seq data from NCBI's SRA. We assume you are working in a Linux/Unix or Mac environment. 

We'll need to set up a few things before we can begin.

## Aspera ascp installation
Aspera's `ascp` is a command-line fasp transfer program that's quite a bit faster than FTP. This is the utility `srafish.pl` uses to download the raw data from the SRA. To install:

1. Select the appropriate download of aspera connect from http://downloads.asperasoft.com/en/downloads/8?list
2. For linux and mac, this should be a shell script. Make it executable `chmod +x aspera-connect-3.5.1.92523-linux-64.sh`
3. Run it. It should install in your home directory. `./aspera-connect-3.5.1.92523-linux-64.sh`

## Get `srafish.pl`
To get `srafish.pl` you'll need to download our GitHub repository. If you have Git installed, you can grab our code by entering the following in your command line: `git clone https://github.com/surgebiswas/transcriptome_compression.git`. Otherwise, you can grab the repository manually by navigating to `https://github.com/surgebiswas/transcriptome_compression` and clicking "Clone in Desktop" or "Download Zip" (and then unzipping the archive).

Now add the `data_download` folder of the repository to your `$PATH` variable. For example, if the main repository directory is `~/GitHub/transcriptome_compression`, add the line `export PATH=$PATH:~/GitHub/transcriptome_compression/data_download` to your `~/.bash_profile` or your `~/.bashrc`. Then enter `source ~/.bash_profile` or `source ~/.bashrc` as the case may be. 

## Build a query table
The main argument to `srafish.pl` is a query table that contains, among other things, a list of SRA sample IDs that `srafish.pl` will download. To build this table, enter the following in your command line:

```
qt_name=<query_table_file_name>
sra_url=http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=
organism=<organism_name>
wget -O $qt_name `$url($organism[Organism]) AND "strategy rna seq"[Properties]'
```

# Misc. notes
1. If getting the error: "name not found while resolving tree within virtual file system module" Try updating sra-toolkit to latest release.
