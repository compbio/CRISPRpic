# CRISPRpic

A Python script for fast and precise analysis of CRISPR-induced mutations via prefixed index counting (CRISPRpic)

**CRISPRpic runs on Python 2.7 or later. You can simply download CRISPRpic.py and run it without further installation.**
Make sure to put CRISPRpic.py in your path, so it is executable anywhere. Otherwise, place the script in your working directory.

The script has 2 python package dependencies:
* pandas
* matplotlib

For dependency installation instructions, see Dependencies section (below)

## Test

You can test CRISPRpic.py by running it on the example input files (AAVS1.out.extendedFrags.fastq and AAVS1_input.txt), located in the TEST/ directory.

Command for test:
```
python CRISPRpic.py -i AAVS1_input.txt -f AAVS1.out.extendedFrags.fastq -w 3
```
## Command line usage
```
python CRISPRpic.py -i INPUT -f SEQFILE -w WINDOW
```

**-i INPUT file** contains the following information seperated by tab (\t):
* Locus name such as TP53
* Expected amplicon sequence
* guide RNA seq with PAM site (for enzyme_type 1 and 2) or break point from the 5' end amplicon (for 3)
* the type of enzyme - 1:SpCas9, 2:AsCpf1, 3:Custom

You can find an example input file in TEST/AAVS1_input.txt

**-f SEQFILE** is a fastq file of single-end sequencing data. If you have paired-end sequencing data, you can merge to single-end using  a program called Fast Length Adjustment of SHort reads (FLASH: https://ccb.jhu.edu/software/FLASH/)

You can find an example file at TEST/AAVS1.out.extendedFrags.fastq 

**-w WINDOW** is the size of the mutagenic window from the double strand break (DSB). -w 3 means that we only consider a mutation within 3 bp of both directions from the DSB while other mutations outside of this window will be considered as unmodified.

**-d INDEX_SIZE (optional)** is the starting size of index. Default is 8, but larger size such as 12 should be used for when the amplicon contains lots of homologous or low complexity sequences


## How to interprete outputs
*

#### Dependencies

On Ubuntu or Debian Linux:
```
sudo apt-get install python-matplotlib python-pandas
```
On Mac OS X:
```
conda install pandas matplotlib
```
On Windows 10 using Anaconda Prompt:
```
conda install pandas matplotlib
```
Anaconda Prompt can be installed on Window10 by the following instruction:
https://conda.io/docs/user-guide/install/windows.html

