# CRISPRpic

A Pyrhon script for fast and precise analysis of CRISPR-induced mutations via prefixed index counting (CRISPRpic)

CRISPRpic runs on Python 2.7 or later. You can simply download CRISPRpic.py and run it without further installation.
Only dependecies are Pandas and matplotlib, which can be easily installed.

On Ubuntu or Debian Linux:
```
sudo apt-get install python-matplotlib python-pandas
```
On Mac OS X:
```
conda install pandas matplotlib
```
On Window10 using Anaconda Prompt:
```
conda install pandas matplotlib
```
Anaconda Prompt can be installed on Window10 by the following instruction:
https://conda.io/docs/user-guide/install/windows.html

## Command line usage
python CRISPRpic.py -i INPUT -s SEQFILE -w WINDOW

-i INPUT file contains the following information seperated by tab (\t):
Locus name such as PVT1
Expected amplicon sequence
guide RNA seq with PAM site (for 1 and 2) or break point from the 5' end amplicon (for 3)
the type of enyme - 1:SpCas9, 2:AsCpf1, 3:Custom

You can find these example cases in example_input.txt

-s SEQFILE is the either single end sequencing data or merged pair end sequencing data by FLASH (see bleow)

You can find an example file, AAVS1.out.extendedFrags.fastq in the test/ directory.

-w WINDOW is the size of mutagenic window from the double strand break (DSB). -w 3 means that we only consider a mutation within 3 bp for both direction from the DSB while other mutations outside of this window will be considered as unmodified.

-d INDEX_SIZE is the starting size of index. Default is 8, but larger size such as 12 should be used for when the amplicon contains homo

## Test

You can run CRISPRpic.py by running it on the example file in the test/ directory.
Download CRISPRpic.py and AAVS1.out.extendedFrags.fastq into a folder.

python CRISPRpic.py -i cas9_list.txt -w 3 -s AAVS1.out.extendedFrags.fastq

AAVS1.out.extendedFrags.fastq is the output of a program called "FLASH (Fast Length Adjustment of SHort reads)"

More information and installation can be found:
https://ccb.jhu.edu/software/FLASH/



## How to run script


