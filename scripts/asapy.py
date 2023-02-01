#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Launch the MolD GUI"""

import multiprocessing

from itaxotools.asapy.gui import run

if __name__ == '__main__':
    multiprocessing.freeze_support()
    run()
