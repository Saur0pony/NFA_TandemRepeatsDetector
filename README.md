# NFA Tandem Repeat Detector v5 + MaxRepv3

This is the fifth version of pattern matching which use NFA algorithm to find occurrences of a tandem repeat pattern in a DNA sequence.
Algorithm come from Baeza-Yates and was improved by Keikki Hyyrö and rewrite in python by Grégoire Prunier to be use in case of Tandem Repeat.
It can be used to find Tandem repeats or CRISPR with possible addition of a gap.

MaxRepv3 is a programme which can be used to find maximal repeat which to get some pattern to test with Tandem Repeat Detector

## Note about the version

This is the fifth version with a lot of adding!
v1 - Simply implement the NFA algorithm
v2 - Try to implement a detection of tandem
v3 - go back with spacer from the 1st version of Baeza-Yates
v4 - Implementation of tandem repeat detection and pattern occurrence count + parallel processes (cut sequence in smaller pieces)
v5 - spacer suppression and upgrade different things

## Getting Started

The program is containes in the python script NFA_TandemRepeat_detector_v5.py 

"onelinedfasta.py" is used to transform a fasta file in which the sequence is written in multiple line. That put all the sequence lines in one line under the sequence ID.

The "run.sh" is used to find TR of a pattern given in argument, in a simple fasta sequence or multifasta.

### Prerequisites

This algorithm need python 3  (tested with python 3.8.5)
for use with pypy like with the "run.sh", the package pypy is mandatory, and must be in a stable version of python (python 3.7.9)

Required Packages for TR Detector:

* multiprocessing
* getopt
* sys
* threading

Package required for MaxRep:

* Sufarray-kkto
* Numpy

Sequences must be in one line.
To do that, use "onelinedfasta.py" script.
```
python onelinedfasta.py --file downloadedsequence.fasta --out sequence_onelined.fasta
```


## Use the algorithms
### Tandem Repeat Detector
The program simply needs to launch the python script and enter different arguments.

* `--file` : Input sequence file (can contain multiple sequences). It is needed to be in Fasta format.
* `--out` : The output name file in which are stored results and informations about used arguments.
* `--pattern` : The pattern you want to search in sequences.
* `--k` : The maximum number of errors allowed
* `--gap` : size of th gap between 2 occurrences of the pattern (case of CRISPR for example)

example of execution command line :
```
python pattern_matching_NFA.py --file sequences.fa --out result.txt --pattern TGT --k 1
```
### Use the algorithm Tandem repeat detector with run.sh
This method need to have pypy package.

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
### MaxRep

Simply need to launch python script with different arguments :

* `--file` : Input sequence file (can contain multiple sequences). It is needed to be in Fasta format (onelined).
* `--out` : The output name file in which are stored results and informations about argument uses.
* `--ml` : Minimal length of maximal repeats

example of execution command line :
```
python --file Human_genome.fna --out output.out --ml 20
```

## Versioning

This is the 5th version of the algorithm.
It takes a pattern and search all tandem repeats with this pattern. Each tandem repeat is a repetition of the same pattern side to side.
It takes a k value of maximal error allowed in 1 pattern.

## Authors

* **Grégoire Prunier** 
* **Jacques Nicolas**


## Acknowledgments

* Baeza-Yates, Navarro 1996
* Heikki Hyrrô 2006 and 2008
* Becher, Deymonnaz, Ariel Heber 2012
