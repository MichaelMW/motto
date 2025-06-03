# Motto
Represent motifs in consensus sequences

## Requirement
Any Python (developed in python v3.7; backward compatible with python v2.7)

## Getting started
```
git clone https://github.com/MichaelMW/motto.git
cd motto
./test.sh

### examples
## default
cat data/ctcf.meme | ./motto.py
# CTCF_HUMAN.H10MO.A	NNNCC[AG][CG][CTA]AG[GA][GT]GGC[GA][CG][TC][AGC]N

## output style to IUPAC (see http://www.bioinformatics.org/sms/iupac.html)
cat data/ctcf.meme | ./motto.py -s IUPAC
# CTCF_HUMAN.H10MO.A	NNNCCRSHAGRKGGCRSYVN

## trim off flanking Ns.
cat data/ctcf.meme | ./motto.py -s IUPAC --trim
# CTCF_HUMAN.H10MO.A	CCRSHAGRKGGCRSYV

```
## Preview
Run Motto instantly in your browser at
[michaelmw.github.io/motto](https://michaelmw.github.io/motto/).

## Options
```
usage: motto.py [-h] [-i INMEME] [-m METHOD] [-s STYLE] [-d DELIMITER]
                [-p PENALTY] [-M MAXCHARACTER] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -i INMEME, --inMeme INMEME
                        Input meme format motif file
  -m METHOD, --method METHOD
                        Method: Motto, MSE, Cavener, Max
  -s STYLE, --style STYLE
                        Output style: compact, regex, IUPAC (see:
                        http://www.bioinformatics.org/sms/iupac.html)
  -d DELIMITER, --delimiter DELIMITER
                        Optional delimiter. Eg. -d '_' gives C_T_[AG]; use -d
                        $' ' for special chars.
  -p PENALTY, --penalty PENALTY
                        Penalty for ambiguity at each position [0,1]. Only
                        works with -m Motto or MSE. Set to 0 (default) for no
                        penalty. Set to 1 will give maximal penalty, and gives
                        a single nucleotide with max frequency as the
                        consensus.
  -M MAXCHARACTER, --maxCharacter MAXCHARACTER
                        Max number of characters allowed for consensus at each
                        position, only works when -m is set to Motto or MSE.
                        Default: None, no restriction.
  -t, --trim            Trim off ambiguous consensus characters (eg. [ACTG] or
                        N) at both ends

```


## Features
### Overview
![Overview](https://github.com/MichaelMW/motto/blob/master/figures/Fig1.overview.png "Overview")

### Example Usage
![Example usage: CTCF](https://github.com/MichaelMW/motto/blob/master/figures/Fig2.1.ctcf.png "Example usage: CTCF")
![Example usage: lipoprotein](https://github.com/MichaelMW/motto/blob/master/figures/Fig2.2.lipoAA.png "Example usage: lipoprotein")
### Benchmark
![Benchmark](https://github.com/MichaelMW/motto/blob/master/figures/Fig3.benchmark.png "Benchmark")

## Supplementary webpage
[Webpage](http://wanglab.ucsd.edu/star/motto/ "Webpage")

## Citation
Mengchi Wang, David Wang, Kai Zhang, Vu Ngo, Shicai Fan, Wei Wang, 
"Motto: Representing Motifs in Consensus Sequences with Minimum Information Loss." 
Genetics, Volume 216, Issue 2, 1 October 2020, Pages 353–358, 
[Genetics](https://doi.org/10.1534/genetics.120.303597)
