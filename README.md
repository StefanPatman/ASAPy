# ASAPy

Assemble species by automatic partitioning.

This is a Python wrapper for ASAP: <https://bioinfo.mnhn.fr/abi/public/asap/>

Currently using source code from 27 October 2022.

*(you will need a C compiler when building from source)*

## Quick start

Install using pip:

```
$ pip install .
```

Run the GUI:

```
$ asapy-qt
```

Simple command line tool:

```
$ asapy tests/test.fas
```

## Launch without installing

Before the first time you use the program, you must install any required modules, build the ASAP core and auto-compile the Qt resource files:
```
$ pip install -r requirements.txt
$ python setup.py build_ext --inplace
$ python setup.py build_qt
```

You can now launch the GUI:
```
$ python launcher.py
```

## Packaging

You must first compile the ASAP C module, auto-compile Qt resources,
then use PyInstaller on the launcher **spec** file:
```
$ pip install pyinstaller
$ python setup.py build_ext --inplace
$ python setup.py build_qt
$ pyinstaller launcher.spec
```

## Module

You may import and use the ASAP module in your python scripts.

To launch the GUI:
```
>>> import asapy.qt
>>> asapy.qt.main.show()
```

More examples to follow soon.

### Python interactive example

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

## Acknowledgements

Puillandre N, Brouillet S, Achaz G. ASAP,\
assemble species by automatic partitioning,
Molecular Ecology Resources 2021.

- Original C code by Guillaume Achaz and Sophie Brouillet
- BIONJ by Olivier Gascuel
