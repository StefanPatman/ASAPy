"""A setuptools based setup module."""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages, Extension, Command
import pathlib

here = pathlib.Path(__file__).parent.resolve()


asapmodule = Extension('itaxotools.asapy.asap',
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

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='asapy',
    version='0.1.2',
    description='A Python wrapper for ASAP',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Patmanidis Stefanos',
    author_email='stefanpatman91@gmail.com',
    classifiers = [
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    package_dir={'': 'src'},
    packages=find_namespace_packages(
        include=('itaxotools*',),
        where='src',
    ),
    ext_modules = [asapmodule],
    python_requires='>=3.8.6, <4',
    install_requires=[
        'itaxotools-common==0.2.4',
        'PySide6>=6.1.3',
    ],
    extras_require={
        'dev': [
            'pyinstaller'
        ]
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'asapy=itaxotools.asapy.run:main',
            'asapy-gui=itaxotools.asapy.gui:run',
        ],
        'pyinstaller40': [
            'hook-dirs = itaxotools.__pyinstaller:get_hook_dirs',
            'tests = itaxotools.__pyinstaller:get_pyinstaller_tests',
        ],
    },
)
