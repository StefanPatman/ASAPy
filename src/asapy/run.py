"""Console entry point"""

import sys
from . import core

#! This should be expanded to accept all arguments
def main():
    """Anayze given file"""
    if len(sys.argv) == 2:
        print(' ')
        a=core.PartitionAnalysis(sys.argv[1])
        a.launch())
    else:
        print('Usage: asap FILE')
        print('Ex:    asap tests/test.fas')
