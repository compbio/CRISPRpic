# CRISPRpic

A Pyrhon script for fast and precise analysis of CRISPR-induced mutations via prefixed index counting (CRISPRpic).

CRISPRpic runs on Python 2.7 or later. 
Only dependecies are Pandas and matplotlib, which can be easily installed.

On Ubuntu or Debian Linux:
'''
sudo apt-get install python-matplotlib python-pandas
'''
On Mac OS X:
'''
conda install pandas matplotlib
'''
On Window10 using Anaconda Prompt:
'''
conda install pandas matplotlib
'''
Anaconda Prompt can be installed on Window10 by the following instruction:
https://conda.io/docs/user-guide/install/windows.html

## Test

You can run CRISPRpic.p by running it on the example filein the test/ directory.
Download CRISPRpic.py and AAVS1.out.extendedFrags.fastq into a folder.

python CRISPRpic.py -i cas9_list.txt -w 3 -s AAVS1.out.extendedFrags.fastq

AAVS1.out.extendedFrags.fastq is the output of a program called "FLASH (Fast Length Adjustment of SHort reads)"

More information and installation can be found:
https://ccb.jhu.edu/software/FLASH/



## How to run script


