from setuptools import setup, Extension
from pybind11.setup_helpers import build_ext

ext_modules = [
    Extension(
        'norm_seq.normalization',
        ['norm_seq/normalization.cpp'],
        include_dirs=['pybind11/include'],
        language='c++'
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
