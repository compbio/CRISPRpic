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

INPUT file contains the following information seperated by tab (\t):
Locus name such as PVT1
Expected amplicon sequence
guide RNA seq with PAM site (for 1 and 2) or break point from the 5' end amplicon (for 3)
the type of enyme - 1:SpCas9, 2:AsCpf1, 3:Custom

For instance, AAVS1_input.txt in the test/ directory has the following cotents:
AAVS1	TTCTGGGAGAGGGTAGCGCAGGGTGGCCACTGAGAACCGGGCAGGTCACGCATCCCCCCCTTCCCTCCCACCCCCTGCCAAGCTCTCCCTCCCAGGATCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAATGGGGGTGTGTCACCAGATAAGGAATCTGCCTAACAGGAGGTGGGGGTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGC	TGCTTACGATGGAGCCAGAGagg	1



## Test

You can run CRISPRpic.py by running it on the example file in the test/ directory.
Download CRISPRpic.py and AAVS1.out.extendedFrags.fastq into a folder.

python CRISPRpic.py -i cas9_list.txt -w 3 -s AAVS1.out.extendedFrags.fastq

AAVS1.out.extendedFrags.fastq is the output of a program called "FLASH (Fast Length Adjustment of SHort reads)"

More information and installation can be found:
https://ccb.jhu.edu/software/FLASH/



## How to run script


