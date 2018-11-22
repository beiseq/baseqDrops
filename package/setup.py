# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'readme.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='baseqDrops',
    version='1.5',
    description='Processing Drop-seq, 10X(3prime) and inDrop RNA-seq dataset',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gene.pku.edu.cn',
    author='Xiannian Zhang',
    author_email='friedpine@gmail.com',

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],

    keywords = 'sample setuptools development',

    packages = find_packages(include=['baseq*'], exclude=[]),

    install_requires=[
        'pandas',
        'pysam',
        'numpy',
        'click',
        'configparser',
        'matplotlib'
    ],

    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    entry_points={
        'console_scripts': [
            'baseqDrops=baseqDrops.cmd:cli',
        ]
    }
)