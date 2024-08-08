Installation
============

As of version ``4.15.0.0``, ``PACK`` uses autotools_.  This means that the
typical ``./configure && make && make install`` routine should work on just
about any modern \*nix distribution, assuming autotools is installed.

.. _autotools: http://www.gnu.org/software/hello/manual/automake/Autotools-Introduction.html

MPI
***

``PACK`` supports the use of MPI_ for improved performance on multi-core or
multi-cpu based setups.  By default when ``./configure`` is run, an attempt
will be made to determine if a compatible MPI library is installed.  If you do
not want ``PACK`` to be compiled with MPI support, you can force a serial
version to be compiled by typing ``./configure --with-mpi=no``.

.. _MPI: http://en.wikipedia.org/wiki/Message_Passing_Interface