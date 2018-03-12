===========
FemtoFitter
===========

Tools for fitting femtoscopic correlation functions


About
-----

The art and science of femtoscopy relies on the accurate fitting of correlation functions to data
gathered by various collider experiments. 
This package hopes to provide both a standard and easy-to-use base on which fits may be performed.

There are two components of this library:

  libFemtoFitter.so
    A C++ library which interfaces with ROOT for reading and writing .root files and performing
    fits with TMinuit objects, the standard for fitting femtoscopic analysis data.

  femtofitter
    A python package which wraps libFemtoFitter with an easy to use interface


Installation
------------

This library can be built from source or installed via the python package manager ``pip`` which bundles
both the libFemtoFitter.so library and femtofitter python package.

Requirements
~~~~~~~~~~~~

Building requires ROOT6 with the following features: python cpp11/cpp14 threads minuit2.

You can check your root installation by running the following command; if it returns all 'yes'
it is correctly configured.

.. code:: bash

    $ root-config --has-python --has-threads --has-minuit2 --has-cpp14

Other Requirements:

  cmake
    Like ROOT itself, femtofitter is configured and built with the program cmake_.
    Minimum version is 3.6

  python
    This should be installed by your operating system; python3.5+ is recommended

  pip
    apt install `python3-pip` 
    There are python <-> ROOT support libraries that will be installed automatically.

From Source
~~~~~~~~~~~

.. code:: bash

    # standard build
    $ make

    # explicit in-source cmake build
    $ cmake . && make


Using pip
~~~~~~~~~

Rather than installing from source, pre-built binaries are provided using the
python-package-index (https://pypi.python.org/pypi/femtofitter) and easily downloadable using the
pip installer (installed via your operating system's installer such as ``sudo apt install python3-pip``

Note: If you are using alternative python installers, such as conda_, please use the standard 
instructions rather than installing potentially conflicting files with ``pip``.

.. code:: bash

    # install from stored binaries (may not work on incompatible machines)
    $ pip install --user femtofitter

    # automatically build from source
    $ pip install --user --no-binary=all femtofitter


Usage
-----

While libFemtoFitter.so *can* be used directly by users' macros and other libraries/executables,
the expected usage is via the ``femtofitter`` executable provided with this package.
This program is designed to simply be handed the names of input files, the classes which shall be used to fit,
and output destinations and it performs all indicated operations.

While flexibility is built into the design, there are limits to the structure of datafiles;
this will be elaborated on more in the documentation. 

The standard usage is to create a configuration file which has three responsibilities:

  * lists the input files and how to group histogram objects in 'analysis-units'
  * indicates which fitters to be used and set all options for this fitter
  * specify the output root-file and how to organize 'output-units'

  
.. code:: yaml

   ---
   - input:
       file: data.root
       matching:
         num: foo/bar
     output:
       file: fit-results-%(timestamp).root
         
     fits:
     # Standard gaussian-fit, no coulomb
     - name: Gauss3D_NoCoulomb
       class: Gauss3D
       coulomb: False

     # Standard gaussian-fit including CoulombFactor provided by
     # the CoulombInterp3D class created with file CoulombInterpData.root
     - name: Gauss3D
       class: Gauss3D
       coulomb:
         class: CoulombInterp3D
         file: CoulombInterpData.root


Development
-----------

Issues, suggestions, and patches are welcome.
Please make sure that submitted code follows the style guidelines enforced by
clang-format and editorconfig files in this directory, and that all unittests
pass when merging into the master branch.


License
-------

This code is released under the conditions of the LGPL-2.1_ free-software license, the
contents of which are provided in the LICENSE file of this repository.
Copyright is held by Andrew Kubera (mailto:andrew.michael.kubera@cern.ch).


.. _conda: https://conda.io/docs/
.. _LGPL-2.1: https://opensource.org/licenses/LGPL-2.1

