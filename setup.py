#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 6, 0):  # 3.7.1 used
    IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
    sys.exit('Upgrade your Python to >= 2.7.7 or 3.6+; version=(%s.%s.%s)' % (IMAJOR, MINOR1, MINOR2))


packages = find_packages() # exclude=['ez_setup', 'examples', 'tests'] + exclude_words

IS_TRAVIS = 'TRAVIS' in os.environ
#is_rtd = 'READTHEDOCS' in os.environ

install_requires = ['numpy']
IS_WINDOWS = 'nt' in os.name
if IS_TRAVIS and not IS_WINDOWS:
    #install_requires.append('python-coveralls')
    install_requires.append('codecov')
    #install_requires.append('coverage')


# get package metadata
import pyatmos
setup(
    name='pyatmos',
    version=pyatmos.__version__,
    description=pyatmos.__desc__,
    long_description=pyatmos.__long__,
    classifiers=[
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords='',
    #'>2.7.6,!=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*',
    python_requires='>3.5',
    author=pyatmos.__author__,
    author_email=pyatmos.__email__,
    url=pyatmos.__website__,
    license=pyatmos.__license__,
    packages=packages,
    include_package_data=True,
    zip_safe=True,
    install_requires=install_requires,
    #{'': ['license.txt']}
    #package_data={'': ['*.png']},
    #data_files=[(icon_path, icon_files2)],
    package_data={
        # https://pythonhosted.org/setuptools/setuptools.html#including-data-files
        # If any package contains *.png files, include them:
        '': ['*.png'],
        #'mypkg': ['data/*.dat'],
    },
    entry_points={
        'console_scripts': [
        ]
    },
    test_suite='pyatmos.utils.test_atmosphere',
)
