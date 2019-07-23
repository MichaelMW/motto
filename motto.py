#!/usr/bin/env python
from __future__ import division
from sys import stdin
from signal import signal, SIGPIPE, SIG_DFL
from warnings import filterwarnings
import numpy as np
import argparse
import re
signal(SIGPIPE, SIG_DFL)


def main():
    # parsing input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inMeme",
                        help="Input meme format motif file")
    parser.add_argument("-m", "--method", default="JSD",
                        help="Method: JSD, MSE, Douglas, Max")
    parser.add_argument("-s", "--style", default="compact",
                        help="Output style: compact, regex, IUPAC (see: http://www.bioinformatics.org/sms/iupac.html)")
    parser.add_argument("-d", "--delimiter", default="",
                        help="Optional delimiter. Eg. -d '_' gives C_T_[AG]; use -d $'\t' for special chars.")
    parser.add_argument("-p", "--penalty", default=0,
                        help="Penalty for ambiguity at each position [0,1]. Only works with -m JSD or MSE. Set to 0 (default) for no penalty. Set to 1 will give maximal penalty, and gives a single nucleotide with max frequency as the consensus.")
    parser.add_argument("-M", "--maxAlphabet", default=None,
                        help="Max number of alphabets allowed for consensus at each position, only works when -m is set to JSD or MSE. Default: None, no restriction.")
    parser.add_argument("-t", "--trim", action='store_true',
                        help="Trim off ambiguous consensus alphabets (eg. [ACTG] or N) at both ends")
    args = parser.parse_args()

    # checking input arguments
    method = args.method.lower()
    if method not in ["jsd", "mse", "douglas", "max"]:
        print("Unrecognized method: {}".format(method))
        exit()
    style = args.style.lower()
    if style not in ["compact", "regex", "iupac"]:
        print("Unrecognized style: {}".format(style))
        exit()
    penalty = float(args.penalty)**2  # power for smoother transition
    maxAlphabet = int(args.maxAlphabet) if args.maxAlphabet else None
    trim = args.trim
    if args.inMeme:
        inMeme = args.inMeme
        with open(inMeme) as f:
            inData = f.read()
    else:
        inData = stdin.read()
    delimiter = args.delimiter

    # parsing input meme & init parameters
    alphabet, motifPairs = parseMeme(inData)
    if maxAlphabet:
        refDists = [
            [1/m if i < m else 0 for i in range(len(alphabet))] for m in range(1, maxAlphabet+1)]
    else:
        refDists = [
            [1/m if i < m else 0 for i in range(len(alphabet))] for m in range(1, len(alphabet)+1)]
    IUPACdict = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K",
                 "AC": "M", "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}

    # the main loop
    # for each motif
    for motifID, motifPWM in motifPairs:
        motifKmerList = []
        kmerLists = []
        # get kmerList
        for rates in motifPWM:
            if method == "douglas":
                kmerList = Douglas(rates, alphabet)
            elif method == "jsd":
                kmerList = rates2kmer(rates, refDists, alphabet, penalty)
            elif method == "mse":
                kmerList = rates2kmer_mse(rates, refDists, alphabet, penalty)
            elif method == "max":
                kmerList = rates2kmer_max(rates, alphabet)
            else:
                print(("Unrecognized method: {}".format(method)))
                exit()
            kmerLists.append(kmerList)
        # check trim
        if trim:
            while len(kmerLists[0]) == len(alphabet):
                del kmerLists[0]
            while len(kmerLists[-1]) == len(alphabet):
                del kmerLists[-1]
        # convert to style
        for kmerList in kmerLists:
            kmer = kmerStyle(kmerList, style, IUPACdict, alphabet)
            motifKmerList.append(kmer)
        print(("{}\t{}".format(motifID, delimiter.join(motifKmerList))))
    pass


def kmerStyle(kmerList, style, IUPACdict, alphabet):
    kmerSortedString = ''.join(sorted(kmerList))
    kmerString = ''.join(kmerList)
    kmer = "[{}]".format(kmerString)
    if len(kmerList) == 1:
        kmer = kmerString
    elif style == "iupac":
        kmer = IUPACdict.setdefault(kmerSortedString, kmer)
    elif style == "compact":
        if len(kmerList) == len(alphabet):
            kmer = "N"
    return kmer

# parsing meme file


def parseMeme(inData):
    header = inData.split("MOTIF")[0]
    motifPairs = []
    for motifInfo in inData.split("MOTIF")[1:]:
        motifID = motifInfo.split('\n')[0].strip()
        motifPWM = list(list(map(float, line.split())) for line in motifInfo.split(
            '\n')[2:] if re.match("^\s*[.0-9]", line))
        motifPairs.append([motifID, motifPWM])
    alphabet = ''.join([line.lstrip("ALPHABET=").strip()
                        for line in header.split("\n") if re.match("ALPHABET*", line)])
    pattern = '^[{}]'.format(alphabet)
    return alphabet, motifPairs

# convert a list of rates to regex style kmer at one position


def rates2kmer(L, refDists, alphabet, penalty):
    keptInd = indMostInfo(L, refDists, penalty)
    return [alphabet[i] for i in keptInd]


def rates2kmer_mse(L, refDists, alphabet, penalty):
    keptInd = indMSE(L, refDists, penalty)
    return [alphabet[i] for i in keptInd]


def rates2kmer_max(L, alphabet):
    sList, sIndex = list(
        zip(*sorted([(l, i) for i, l in enumerate(L)], reverse=True)))
    return [alphabet[sIndex[0]]]

# Douglas 1976 heuristic method


def Douglas(L, alphabet):
    if len(L) != 4:
        print(("input format error! Has to be nucleotides! len(L)={}".format(len(L))))
        exit()
    else:
        sList, sIndex = list(
            zip(*sorted([(l, i) for i, l in enumerate(L)], reverse=True)))
        if sList[0] > 0.5 and sList[0] > sList[1]*2:
            return [alphabet[sIndex[0]]]
        elif sList[0]+sList[1] > 0.75:
            return [alphabet[sIndex[0]], alphabet[sIndex[1]]]
        else:
            return [alphabet[sIndex[i]] for i in range(len(L))]

# return index of L that should be kept, based on min jsd between reference distribution and L


def indMostInfo(L, refDists, penalty):
    sList, sIndex = list(
        zip(*sorted([(l, i) for i, l in enumerate(L)], reverse=True)))
    jsds = [[jsd(sList, refDist) + penalty*penalty*(i+1), i+1]
            for i, refDist in enumerate(refDists)]
    keepN = min(jsds)[1]
    return sIndex[:keepN]

# return index of L that should be kept, based on MSE between reference distribution and L


def indMSE(L, refDists, penalty):
    sList, sIndex = list(
        zip(*sorted([(l, i) for i, l in enumerate(L)], reverse=True)))
    SEs = [[SE(sList, refDist) + penalty*penalty*(i+1), i+1]
           for i, refDist in enumerate(refDists)]
    keepN = min(SEs)[1]
    return sIndex[:keepN]

# Jensen-shannon divergence


def jsd(x, y):
    filterwarnings("ignore", category=RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    d1 = x*np.log(2*x/(x+y))
    d2 = y*np.log(2*y/(x+y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = 0.5*np.sum(d1+d2)
    return d

# Squared errors:


def SE(x, y):
    assert len(x) == len(y)
    d = sum((x[i]-y[i])**2 for i in range(len(x)))
    return d


if __name__ == '__main__':
    main()
