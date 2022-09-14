#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 14:53:05 2021

@author: kirill
"""

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="input mpileup")
parser.add_argument("--output", "-o", help="output picture", default='indels.png')
parser.add_argument("--dpi", "-d", help="dpi", default=150, type=int)
args = parser.parse_args()


insertions = []
deletions = []
coverage = []
position = []

Mpileup = namedtuple('Mpileup', 'chrom pos ref cover read')

def calc(line):
    m = Mpileup(*line.strip().split('\t')[:5])
    position.append(int(m.pos))
    coverage.append(int(m.cover))
    deletions.append(-m.read.count('*') / int(m.cover) if m.cover != '0' else 0)
    insertions.append(m.read.count('+') / int(m.cover) if m.cover != '0' else 0)

if args.input is not None:
    with open(args.input, 'r') as fi:
        for line in fi:
            calc(line)
else:
    for line in sys.stdin:
        calc(line)

insertions = np.array(insertions)
deletions = np.array(deletions)
extra_ticks = []
extra_ticks.extend(np.argwhere(insertions > 0.5).flatten())
extra_ticks.extend(np.argwhere(deletions < -0.5).flatten())
removed = True
while removed:
    for t in reversed(extra_ticks):
        if t - 1 in extra_ticks:
            extra_ticks.remove(t)
            removed = True
            break
    else:
        removed = False

plt.figure(figsize=(20, 15))
ax1 = plt.subplot(211)
plt.plot(position, insertions, label='Insertions', color='blue', alpha=0.8)
plt.plot(position, deletions, label='Deletions', color='red', alpha=0.8)
plt.ylabel('Frequency (#/coverage)')
plt.xlabel('Position')
plt.legend()
plt.xlim(0, 30000)
plt.xticks(list(plt.xticks()[0]) + extra_ticks)
plt.grid(True)
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(position, coverage)
plt.ylabel('Coverage')
plt.xlabel('Position')
plt.grid(True)
plt.savefig(args.output, dpi=args.dpi)
