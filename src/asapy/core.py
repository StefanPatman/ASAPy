#-----------------------------------------------------------------------------
# ASAPy - Assemble Species by Automatic Partitioning with ASAP
# Copyright (C) 2021  Patmanidis Stefanos
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------


from multiprocessing import Process

import tempfile
import shutil
import pathlib
import os
import sys
from contextlib import contextmanager
from datetime import datetime

from . import asap
from . import param
from . import params


@contextmanager
def _redirect(stream='stdout', dest=None):
    """Redirect system stream to file stream"""
    # High-level redirection
    original = getattr(sys, stream)
    original.flush()
    setattr(sys, stream, dest)
    # Low-level redirection
    duplicate = os.dup(original.fileno())
    os.dup2(dest.fileno(), original.fileno())
    try:
        yield dest
    finally:
        # Restore stream
        os.dup2(duplicate, original.fileno())
        dest.flush()
        setattr(sys, stream, original)

@contextmanager
def redirect(stream='stdout', dest=None, mode='w'):
    """
    Redirect system stream according to `dest`:
    - If None: Do nothing
    - If String: Open file and redirect
    - Else: Assume IOWrapper, redirect
    """
    if dest is None:
        yield getattr(sys, stream)
    elif isinstance(dest, str):
        with open(dest, mode) as file, _redirect(stream, file) as f:
            yield f
    else:
        with _redirect(src, tar) as f:
            yield f


class PartitionAnalysis():
    """
    Container for input/output of ASAP core.
    """

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state

    def __init__(self, file):
        """
        """
        self.file = file
        self.useLogfile = False
        self.target = None
        self.results = None
        # self.time_format = '%x - %I:%M%p'
        self.time_format = '%FT%T'
        self.param = param.ParamList(params.params)

    def fetch(self, destination):
        """
        Copy results as a new directory.
        """
        if self.results is None:
            raise RuntimeError('No results to fetch.')
        shutil.copytree(self.results, destination)

    def run(self):
        """
        Run the ASAP core with given params,
        save results to a temporary directory.
        """
        kwargs = self.param.as_dictionary()
        kwargs['time'] = datetime.now().strftime(self.time_format)
        if self.target is not None:
            kwargs['out'] = self.target
        # When running as Windows GUI app, stdout and stderr are NullWriters
        # that don't define fileno(). In this case, freopen must be called
        # before redirecting or "Bad file descriptor" error will occur.
        if (not hasattr(sys.stdout, 'fileno') or
            not hasattr(sys.stderr, 'fileno')):
            out, err = asap._freopen()
            sys.stdout = open(os.devnull, 'w')
            sys.stderr = open(os.devnull, 'w')
            sys.stdout.fileno = lambda: out
            sys.stderr.fileno = lambda: err
        with open(str(pathlib.Path(self.target) / "asap.log"), 'w') as logfile, \
                _redirect('stdout', logfile), _redirect('stderr', logfile):
            asap.main(self.file, **kwargs)
        self.results = self.target

    def launch(self):
        """
        Should always use a seperate process to launch the ASAP core,
        since it uses exit(1) and doesn't always free allocated memory.
        Save results on a temporary directory, use fetch() to retrieve them.
        """
        # When the last reference of TemporaryDirectory is gone,
        # the directory is automatically cleaned up, so keep it here.
        self._temp = tempfile.TemporaryDirectory(prefix='asap_')
        self.target = pathlib.Path(self._temp.name).as_posix()
        p = Process(target=self.run)
        p.start()
        p.join()
        if p.exitcode != 0:
            raise RuntimeError('ASAP internal error, please check logs.')
        # Success, update analysis object for parent process
        self.results = self.target
