#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print(f"usage: {sys.argv[0]} <exact> <approx>")
    exit(1)

def load_fingerprints(filename: str) -> set[int]:
    fps = set()
    with open(filename, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        row = line.split(",")
        if len(row) > 1:
            fps.add(int(row[0]))

    return fps

# load exact top fingerprints
exact = list()

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

for line in lines:
    row = line.split(",")
    if len(row) > 1:
        exact.append(int(row[0]))

# load approximation
approx = set()
with open(sys.argv[2], "r") as f:
    lines = f.readlines()

for line in lines:
    row = line.split(",")
    if len(row) > 1:
        approx.add(int(row[0]))

# compute cut
k = len(approx)
if k > len(exact):
    print(f"too few fingerprints in exact top list (have {len(exact)}, need {k})")
    exit(1)

cut = 0
for i in range(0,k):
    if exact[i] in approx:
        cut += 1

score = float(cut) / float(k)
print(f"the approximation contains {cut} out of the top {k} fingerprints ({score * 100}%)")
