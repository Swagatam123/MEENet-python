.. _install:

Installation Guide
==================

Instructions for installing the python library ``MEENet`` and the command line tool of the same name.

.. _require:

Requirements
------------

To run MEENet you will need python >=3.6.0 (we have no plans to support python<=3.5).
We recommend to use `anaconda <https://www.continuum.io/downloads>`_ and the ``conda`` command to install dependencies limiting the use of ``pip`` to libraries not available on conda's main channel. (``pip`` can be used too but its dependency-managing might be less robust).

.. _dependencies:

Dependencies
-----------

Following are the python library dependenccies that needs to be installed as a prerequisite to run our python package.

.. note::
   Node2Vec
   Networkx
   Matplotlib
   pandas
   sklearn
   numpy
   scanpy
   annoy
   seaborn
   umap 

Run the below conda command to install all the dependeices mentioned above. You can also use pip install the python libraries.

::

	conda install Node2vec Networkx Matplotlib pandas sklearn scanpy annoy seaborn umap


Environment Setup
----------------

Once all the dependencies are installed and `git` is installed and configured on your system, run the below command to get the 	python files.

::
	
	git clone https://github.com/Swagatam123/MEENet-python

Download the OpenOrd jar file from this location and make sure it is in the same location where your cloned project is present. Your environment 
setup is ready to run the python file. Please follow the `Turoial<https://swagatam123.github.io/MEENet-python/Tutorial.html>` page to run the code
