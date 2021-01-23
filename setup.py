#!/usr/bin/env python
# Always prefer setuptools over distutils

import sys
from setuptools import setup, find_packages

DS = 'PyAttributeProcessor'
author = 'Sergi Rubio, srubio@cells.es'
url = "https://github.com/alba-synchrotron/pyattributeprocessor"
long_description = url

# If the version is updated automatically with bumpversion
# do not update manually
__version = open(DS+'/VERSION').read().strip()
#if defined raw, it should match with .bumpversion file!
#__version = '4.5.3'

scripts = ['./bin/'+DS,]


__doc__ = """
Generic Device Server setup.py file copied from fandango/scripts/setup.ds.py

To install as system package:

  python setup.py install
  
To build src package:

  python setup.py sdist
  
To install as local package, just run:

  mkdir /tmp/builds/
  python setup.py install --root=/tmp/builds
  /tmp/builds/usr/bin/$DS -? -v4

To tune some options:

  RU=/opt/control
  python setup.py egg_info --egg-base=tmp install --root=$RU/files --no-compile \
    --install-lib=lib/python/site-packages --install-scripts=ds

# windows installer:
# python setup.py bdist_wininst

-------------------------------------------------------------------------------
"""

## All the following defines are OPTIONAL
setup_requirements = []
package_dir = {
    #DS: '.', ## For setup.py located in DS folder 
    ##'DS/tools': './tools', #(with submodules)
}
packages = package_dir.keys() or find_packages()

## Additional files, remember to edit MANIFEST.in to include them in sdist
#package_data = {'': [
  ##'VERSION',
  ##'./tools/icon/*',
  ##'./tools/*ui',
#]}

setup(
    name='tangods-'+DS.lower(),
    description='%s Tango Device Server'%DS,
    version=__version,
    author=author,
    author_email='controls-software@cells.es',
    url = 'https://git.cells.es/controls/'+DS,
    #maintainer=author,
    #maintainer_email=author,    
    #download_url=url,    
    packages=packages,
    #package_data = package_data,      
    include_package_data=True,
    license="GPLv3",
    long_description=long_description,    
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python', # :: 3.5',
        'Topic :: Communications',
        'Topic :: Software Development :: Libraries',
    ],    
    #entry_points={'console_scripts': ['%s=%s.%s:main' % (DS,DS,DS),]},
    install_requires=['pytango','fandango'],
    #python_requires='>=3.5',
    #package_dir= package_dir,          
    scripts = scripts,
)
