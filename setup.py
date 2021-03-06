
import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


setup(
    name='sling',
    version='2.0.1',
    description='sling: a tool to search for linked gene arrays in bacterial datasets',
    packages = find_packages(),
    package_data={'sling': ['data/*']},
    author='Gal Horesh',
    author_email='gh11@sanger.ac.uk',
    url='https://github.com/ghoresh11/sling/wiki',
    scripts=glob.glob('scripts/*'),
    install_requires=[
        'biopython == 1.73',
        'pandas == 0.23',
        'networkx == 2.2',
	]
)
