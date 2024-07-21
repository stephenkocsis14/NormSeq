from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11

ext_modules = [
    Extension(
        'norm_seq.normalization',
        ['norm_seq/normalization.cpp'],
        include_dirs=[
            pybind11.get_include(),
        ],
        language='c++',
        extra_compile_args=['-std=c++11']
    ),
]

setup(
    name='norm_seq',
    version='0.1',
    description='RNA-seq Data Normalization Tool',
    author='Your Name',
    author_email='your.email@example.com',
    packages=['norm_seq'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    entry_points={
        'console_scripts': [
            'norm_seq=cli:main',
        ],
    },
    install_requires=[
        'pandas',
        'pybind11'
    ],
)
