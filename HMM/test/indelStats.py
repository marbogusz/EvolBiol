#!/usr/bin/env python3

sv1 = "ATC--AA-CC-C-CTG----AAAAT"
sv2 = "ATCC-AA-CC-C-CTG--AAAA--T"

st1 = "AT--A-A-CCC-C-TG----AAAAT"
st2 = "ATC-A-A-CC-C-CTG--AAAA--T"

pv = (sv1, sv2)
pt = (st1, st2)

def extractGapHistogram(seqPair):
    s1 = seqPair[0];
    s2 = seqPair[1];
    hist = dict();
    seg_len = 0;
    in_gap = false;

    if s1[0] == '-':
        in_gap = true;
        seg_len += 1;

    i = 1
    while i < len(s1):
        if(in_gap):
            if(s1[i]) 

        

    
def extractFragmentHistogram(seqPair):
    hist = {};
    return hist;


def getAlignmentDistance(trueAl, otherAl):
    #2 tuples of sequences
    i = 0
    while i < len(str):
            print str[i]
                i += 1
