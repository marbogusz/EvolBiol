#!/usr/local/bin/python3

import hmmTools
import sys

if len(sys.argv) != 5:
    print('Usage: script.py seqence_length replicate_no indel_rate model ' , file=sys.stderr)
    raise SystemExit(1)

runner = hmmTools.HmmGenerator(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],True)
runner.run()
alsr = hmmTools.HmmAnalyzer(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4],True)
alsr.run()
