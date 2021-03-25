#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Console entry point"""

import sys
from . import core

#! This should be expanded to accept all arguments
def main():
    """Anayze given file"""
    if len(sys.argv) == 2:
        print(' ')
        a=core.PartitionAnalysis(sys.argv[1])
        a.launch()
    else:
        print('Usage: asapy FILE')
        print('Ex:    asapy tests/test.fas')
