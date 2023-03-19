============
Installation
============


pytransport is developed for Python 3. The developers use 3.8 to 3.11. It can
be install through pip (with internet access): ::

  pip install pytransport

There is also a feature of root-numpy: ::

  pip install pytransport[root-numpy]

which depends on root-numpy for some extra functionality.

Requirements
------------

pytransport depends on the following Python packages available through pip:

 * matplotlib
 * numpy
 * scipy
 * pymadx
 * pybdsim


Local Installation
------------------

Although on pip, for development purposes you may wish to use pytransport from a
copy of the source code. It is possible to clone the git repository and use
pip to `point` at the local set of files, or generally install that set of
files as a once off.

We have provided a simple Makefile in the main pybdsim directory that has
a small set of 'rules' that should help with the relevant pip commands. pip
is used even though pybdsim is found from the local set of files.

To install pybdsim, simply run ``make install`` from the root pybdsim
directory.::

  cd /my/path/to/repositories/
  git clone http://bitbucket.org/jairhul/pytransport
  cd pytransport
  make install

Alternatively, run ``make develop`` from the same directory to ensure
that any local changes are picked up.
