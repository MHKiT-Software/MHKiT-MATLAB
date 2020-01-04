from setuptools import setup, find_packages
from distutils.core import Extension

DISTNAME = 'mhkit_python_utils'
VERSION = '0.1.0'
PACKAGES = ['mhkit_python_utils']
EXTENSIONS = []
DESCRIPTION = ''
LONG_DESCRIPTION = open('README.md').read()
AUTHOR = 'mhkit developers'
MAINTAINER_EMAIL = ''
LICENSE = 'Revised BSD'
URL = ''
CLASSIFIERS=['Development Status :: 0 - Alpha',
             'Programming Language :: Python :: 3',
             'Topic :: Scientific/Engineering',
             'Intended Audience :: Science/Research',
             'Operating System :: OS Independent',
            ],
DEPENDENCIES = ['pandas']


setup(name=DISTNAME,
      version=VERSION,
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      maintainer_email=MAINTAINER_EMAIL,
      license=LICENSE,
      url=URL,
      classifiers=CLASSIFIERS,
      zip_safe=False,
      install_requires=DEPENDENCIES,
      scripts=[],
      include_package_data=True
  )
