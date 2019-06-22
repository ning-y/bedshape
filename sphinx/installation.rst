Installation
============

First, clone the git repository.

.. code-block:: bash

    $ git clone https://github.com/ningyuansg/bedshape.git
    $ cd bedshape

Then, run Make. This requires `pipenv`_, and whatever shapemapper2 requires to
compile its C++ binaries.

.. code-block:: bash

    $ make all

Finally, run bedshape using the executable *bedshape* script.

.. code-block:: bash

    $ bedshape

.. _pipenv: https://github.com/pypa/pipenv
