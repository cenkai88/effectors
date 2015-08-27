# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 22:21:16 2015

@author: cenkai
"""

import csv

class proteome:
    def __init__(self, spe):
        self.proteins=[]
        ind = []
        signalr = []
        with open('C:\\data\\'+spe.upper()+'.unique.fa') as f:
            seq = f.readlines()
            for i in xrange(len(seq)):
               if seq[i][0]=='>':
                   ind.append(i)
            names = [seq[i][1:].strip() for i in ind]
            seqs = [seq[i+1].strip() for i in ind]
        with open('C:\\data\\'+spe.upper()+'_signal.txt') as f:
            signal = f.readlines()
            for sig in signal:
                signalr.append(sig[74])
        for i in xrange(len(names)):
            self.proteins.append(protein(names[i], seqs[i], signalr[i]))
    def filter_e(self, le, c):
        result = []
        for protein in self.proteins:
            if protein.cys > c and protein.length<le and protein.signal=='Y':
                result.append(protein)
        return result

class protein:
    def __init__(self, name, seq, signal):
        self.name = name
        self.seq = seq
        self.cys = seq.upper().count('C')/float(len(seq))
        self.length = len(seq)
        self.signal = signal

def filter_effectors(spe, le, c):
    prot=proteome(spe.upper())
    r = prot.filter_e(le, c)
    with open(spe.upper() + '_filter.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile,dialect='excel')
        spamwriter.writerow(['name', 'cystine_per', 'length', 'signal'])
        for row in r:
            spamwriter.writerow([row.name, row.cys, row.length, row.signal])