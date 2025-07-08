from setuptools import setup, Extension

# An bug with macro definition forces the usage of setup.py
# https://github.com/pypa/setuptools/issues/4810

asapmodule = Extension('itaxotools.asapy._asap',
        include_dirs = ['src/asap'],
        define_macros = [
            ('ismodule', '1')
            ],
        sources = [
            'src/asap/asapmodule.c',
            'src/asap/asap_common.c',
            'src/asap/asap_core.c',
            'src/asap/asap.c',
            'src/asap/draw.c',
            'src/asap/gdtosvg.c',
            'src/asap/oldfns.c',
            'src/asap/wrapio.c',
            ],
        extra_compile_args = ["-w"],
    )

setup(ext_modules = [asapmodule])
