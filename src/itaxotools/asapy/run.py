
"""Console entry point"""

import sys
from pathlib import Path
from . import core

#! This should be expanded to accept all arguments
#! THIS IS CURRENTLY NOT WORKING PROPERLY

def main():
    """Anayze given file"""
    if len(sys.argv) == 2:
        a = core.PartitionAnalysis(sys.argv[1])
        out = Path('out')
        out.mkdir(exist_ok=True)
        a.target = out.as_posix()
        a.run()
    else:
        print('Usage: asapy FILE')
        print('Ex:    asapy tests/test.fas')
