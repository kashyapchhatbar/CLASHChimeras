import logging
import os
import sys

from setuptools import setup

logging.basicConfig(level=logging.DEBUG,
                    format='[%(asctime)s] - %(name)s - %(levelname)s - '
                    '%(message)s')
logger = logging.getLogger('CLASHChimeras')

major, minor = sys.version_info.major, sys.version_info.minor
if major == 3 and minor < 4:
    logger.warn('This package is tested with Python 3.4')
    logger.warn('However it should work with Python version 3.x {}'.format(
        ' but should you face any problems, please install Python >= 3.4'
    ))
elif major < 3:
    logger.error('Please install Python version >= 3 to use this package')
    sys.exit(2)

requirements = open(os.path.join(os.path.dirname(
    __file__), 'requirements.txt')).readlines()
version = open(os.path.join(os.path.dirname(
    __file__), 'VERSION')).readline().rstrip()

setup(
    name="CLASHChimeras",
    version=version,
    author="Kashyap Chhatbar",
    author_email="kashyap.c@ed.ac.uk",
    description="Python package to find chimeras in CRAC/CLASH and"
    " HITS-CLIP datasets",
    install_requires=requirements,
    url='https://github.com/kashyapchhatbar/CLASHChimeras',

    classifiers=['Development Status :: 4 - Beta',

                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Medical Science Apps.',

                 'License :: OSI Approved :: MIT License',

                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',

                 'Operating System :: POSIX :: Linux',
                 'Operating System :: MacOS',

                 'Programming Language :: Python :: 3 :: Only',
                 ],

    keywords='clash chimeras hybrids hits-clip bioinformatics',

    license='MIT',

    py_modules=['clashchimeras.align',
                'clashchimeras.download',
                'clashchimeras.find',
                'clashchimeras.initialize',
                'clashchimeras.log',
                'clashchimeras.methods',
                'clashchimeras.parsers',
                'clashchimeras.runners'],

    scripts=['scripts/align-for-chimeras',
             'scripts/download-for-chimeras',
             'scripts/find-chimeras']

)
