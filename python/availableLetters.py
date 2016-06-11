#!/usr/bin/env python
import sys

alphabet = [chr(i) for i in range(97,123)]
used = list(sys.argv[1].lower())
for c in alphabet:
    if c not in used:
        print c+", ",
