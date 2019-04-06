# Motto
Represent motifs in consensus sequences

## Getting started
```
git clone https://github.com/MichaelMW/motto.git
cd motto
./test.sh
```
## Options
```
usage: motto.py [-h] [-i INMEME] [-m METHOD] [-s STYLE] [-p PENALTY]
                [-M MAXALPHABET] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -i INMEME, --inMeme INMEME
                        Input meme format motif file
  -m METHOD, --method METHOD
                        Method: JSD, MSE, Douglas, Max
  -s STYLE, --style STYLE
                        Output style: compact, regex, IUPAC (see:
                        http://www.bioinformatics.org/sms/iupac.html)
  -p PENALTY, --penalty PENALTY
                        Penalty for ambiguity at each position [0,1]. Only
                        works with -m JSD or MSE. Set to 0 (default) for no
                        penalty. Set to 1 will give maximal penalty, and gives
                        a single nucleotide with max frequency as the
                        consensus.
  -M MAXALPHABET, --maxAlphabet MAXALPHABET
                        Max number of alphabets allowed for consensus at each
                        position, only works when -m is set to JSD or MSE.
                        Default: None, no restriction.
  -t, --trim            Trim off ambiguous consensus alphabets (eg. [ACTG] or
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

## Web service
[Web service](http://wanglab.ucsd.edu/star/motto/ "Web service")

## Citation
Biorxiv: WIP
[Biorxiv]()

