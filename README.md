# ASAPy

**WIP: being prepared for launching on PyPI**

Assemble species by automatic partitioning.

This is a Python wrapper for [ASAP](https://bioinfo.mnhn.fr/abi/public/asap/).</br>
Currently using source code from 27 October 2022.


### Windows and macOS executables
Download and run the standalone executables without installing Python.</br>
[See the latest release here.](https://github.com/iTaxoTools/ASAPy/releases/latest)


### Installing from source
Clone and install the latest version (requires Python 3.8.6 or later):
```
git clone https://github.com/iTaxoTools/ASAPy.git
cd ASAPy
pip install . -f packages.html
```

*(you will need a C compiler when building from source)*


## Usage
To launch the GUI from an installation, please use:
```
asapy-gui
```

Click "Open" to select a fasta or MEGA input file. Configure the ASAP parameters on the left sidebar. Then click "Run" and wait. A list of files will be reported at the end of the analysis. Double-click a file to preview it. Click "Save" to store all results to disk with the suffix of your choice.


### Python module

You may import and use the ASAP module in your python scripts.

From the root directory, launch the Python interpreter:
```
$ python -i
```

Initialize an analysis on your file:
```
>>> a = asapy.PartitionAnalysis('tests/test.fas')
```

Browse and change parameters:
```
>>> a.param.keys()
>>> a.param.general.keys()
>>> a.param.general.all = True
```

Run the analysis:
```
>>> a.launch()
```

You can find the results inside the folder `a.results`.
Save them in a new directory:
```
>>> print(a.results)
>>> a.fetch('./my_results')
```


### Packaging

It is recommended to use PyInstaller from within a virtual environment:
```
pip install ".[dev]" -f packages.html
pyinstaller scripts/asapy.spec
```


## Acknowledgements

Puillandre N, Brouillet S, Achaz G. ASAP,\
assemble species by automatic partitioning,
Molecular Ecology Resources 2021.

- Original C code by Guillaume Achaz and Sophie Brouillet
- BIONJ by Olivier Gascuel
