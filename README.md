# NFA Tandem Repeat Detector v5 + MaxRepv3

This is the fifth version of ppatern matching which use NFA algorythm to find occurences of a tandem repeat pattern in a DNA sequence.
Algorithm come from Baeza-Yates and was improved by Keikki Hyyrö and rewrite in python by Grégoire Prunier to be use in case of Tandem Repeat.
It can be use to find Tandem repeats or CRISPR with addition of gap.

MaxRepv3 is a programme which can be use to find maximal repeat to help te get some motif to test with Tandem Repeat Detector

## Note about the version

This is the fifth version with lot of adds!
v1 - Simply implement the NFA algorithm
v2 - Try to implement a detection of tamdem
v3 - go back whith spacer from 1st version of Baeza-Yates
v4 - Implementation of tandem repeat detection and pattern occurence count + parallelization of process (cut sequence in smaller pieces)
v5 - spacer suppression and upgrade of diferent things

## Getting Started

The program is contain in the python script NFA_TandemRepeat_detector_v5.py 

onelinedfasta.py is used to transforme a fasta file which sequence is write in multiple line. That put all the sequence in one line under the sequence ID.

The run.sh is used to find TR of a pattern gived in argument, in a simple fasta sequence or multifasta.

### Prerequisites

This algorithm need python 3  (tested with python 3.8.5)
for use with pypy like with the run.sh, the package pypy is mandatory, and must be in stable version of python (python 3.7.9)

Packages required :

* multiprocessing
* getopt
* sys
* threading



## Use the algorithm

The program simply need to launch the python script and enter different arguments.

* `--file` : Input sequence file (can contain multiple sequence). It need to be in Fasta format.
* `--out` : The output name where are store results and informations about argument uses.
* `--pattern` : The patern you want to search in sequences.
* `--k` : The maximum number of errors allowed
* `--gap` : size of gap between 2 occurence of motif (case of CRISPR for example)

example of execution command line :
```
python pattern_matching_NFA.py --file sequences.fa --out result.txt --pattern TGT --k 1
```
## Use the algorithm with run.sh
Arguments are simplier to enter :

- `-f = --file`
- `-o = --out`
- `-p = --pattern`
- `-k = --k`
- `-g = --gap`


example of execution command line :
```
./run.sh -f Human_genome.fna -p GTGTCCCCGCGCCAGC -k 5 -o result.txt
```

## Versioning

This is th 5th version of the algorythm.
It take a pattern and search all tandem repeats with this pattern. Each tandem repeat is repetition of the same pattern side to side.
It take a k value of maximal error allowed in 1 pattern.

## Authors

* **Grégoire Prunier** 



## Acknowledgments

* Baeza-Yates, Navarro 1996
* Heikki Hyrrô 2006 and 2008
*
