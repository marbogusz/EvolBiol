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

def gapsFromSeq(s1, dct):
    seg_len = 0;
    in_gap = false;

    if s1[0] == gap_char:
        in_gap = true;
        seg_len += 1;
    i = 1
    while i < len(s1):
        if in_gap:
            if s1[i] == gap_char :
                seg_len +=1
            else:
                in_gap = false
                #add to dict
                dct[seg_len] +=1
                seg_len = 0
        else
            if s1[i] == gap_char :
                in_gap = true;
                seg_len +=1

def fragsFromSeq(s1, dct):
    seg_len = 0;
    in_seg = false;

    if s1[0] != gap_char:
        in_gap = true;
        seg_len += 1;
    i = 1
    while i < len(s1):
        if in_seg :
            if s1[i] != gap_char :
                seg_len +=1
            else:
                in_seg = false
                #add to dict
                dct[seg_len] +=1
                seg_len = 0
        else
            if s1[i] != gap_char :
                in_seg = true;
                seg_len +=1

def getScore(ch1,ch2):
    if ch1 == gap_har:
        return -1;
    else if ch2 == gap_char:
        return 1;
    else return 0;

def getAlignmentDistance(trueAl, otherAl):
    #2 tuples of sequences
    i = 0
    while i < len(str):
            print str[i]
                i += 1

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

hist = Counter()        
gapsFromSeq(s1v, hist)

knights = {'gallahad': 'the pure', 'robin': 'the brave'}
>>> for k, v in knights.iteritems():
    ...     print k, v

