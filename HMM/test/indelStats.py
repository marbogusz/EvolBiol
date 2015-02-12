#!/usr/bin/env python3

sv1 = "ATC--AA-CC-C-CTG----AAAAT"
sv2 = "ATCC-AA-CC-C-CTG--AAAA--T"

st1 = "AT--A-A-CCC-C-TG----AAAAT"
st2 = "ATC-A-A-CC-C-CTG--AAAA--T"

pv = (sv1, sv2)
pt = (st1, st2)

gap_char = '-'

class Counter(dict):
    def __missing__(self, key):
        return 0

def extractGapHistogram(seqPair):
    s1 = seqPair[0];
    s2 = seqPair[1];
    hist = Counter();
    #segmant length
    seg_len = 0;
    in_gap = false;

    if s1[0] == gap_char:
        in_gap = true;
        seg_len += 1;
    i = 1
    while i < len(s1):
        if(in_gap):
            if(s1[i] == gap_char):
                seg_len +=1
            else:
                in_gap = false
                #add to dict
                hist[seg_len] +=1
                seg_len = 0
        else
            if(s1[i] == gap_char):
                in_gap = true;
                seg_len +=1

    #do the same for s2
    seg_len = 0;
    in_gap = false;

    if s2[0] == gap_char:
        in_gap = true;
        seg_len += 1;
    i = 1
    while i < len(s2):
        if(in_gap):
            if(s2[i] == gap_char):
                seg_len +=1
            else:
                in_gap = false
                #add to dict
                hist[seg_len] +=1
                seg_len = 0
        else
            if(s2[i] == gap_char):
                in_gap = true;
                seg_len +=1
    #do the same for s2

    return hist;   

        

    
def extractFragmentHistogram(seqPair):
    s1 = seqPair[0];
    s2 = seqPair[1];
    hist = Counter();
    #segmant length
    seg_len = 0;
    in_gap = false;

    if s1[0] == gap_char:
        in_gap = true;
    i = 1
    while i < len(s1):
        if(in_gap):
            if(s1[i] != gap_char):
                in_gap = false
                seg_len +=1;
        else
            if(s1[i] == gap_char):
                in_gap = true
            else:



def getAlignmentDistance(trueAl, otherAl):
    #2 tuples of sequences
    i = 0
    while i < len(str):
            print str[i]
                i += 1
