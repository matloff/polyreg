cfanalytics: Downloading, analyzing and visualizing CrossFit data
=================================================================

.. image:: https://travis-ci.org/raybellwaves/cfanalytics.svg?branch=master
   :target: https://travis-ci.org/raybellwaves/cfanalytics
.. .. image:: https://ci.appveyor.com/api/projects/status/github/raybellwaves/cfanalytics?svg=true&passingText=passing&failingText=failing&pendingText=pending
.. ..  :target: https://ci.appveyor.com/project/raybellwaves/cfanalytics
.. .. image:: https://coveralls.io/repos/github/raybellwaves/cfanalytics/badge.svg?branch=master
.. ..  :target: https://coveralls.io/github/raybellwaves/cfanalytics?branch=master
.. image:: https://img.shields.io/pypi/v/cfanalytics.svg
   :target: https://pypi.python.org/pypi/cfanalytics
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1208106.svg
   :target: https://doi.org/10.5281/zenodo.1208106
   
**cfanalytics** is an open source project and Python package that aims to provide analyzes to 
CrossFit® workouts. The goal is to enhance data-driven performance of athletes.

.. image:: https://github.com/raybellwaves/cfanalytics/blob/master/Plots/Men_Rx_2018_Overall_rank_P0.1.png
.. image:: https://github.com/raybellwaves/cfanalytics/blob/master/Plots/Women_Rx_2018_Miami_Overall_rankP0.png

Installing
----------

The version numbers 0.1.X are all development versions as I chip away at the project on here in my spare time. I occasionally build the package if i've finished something large. But you can grab the least-buggy version of this package by typing:

``pip install git+https://github.com/raybellwaves/cfanalytics``

Make sure you have the optional dependencies installed below for this to work.

The psudo-stable packages are here.

``conda install -c conda-forge cfanalytics``

``pip install cfanalytics``

As a precautionary note, it has been developed entirely on Mac OSX using Python installed with `anaconda <https://anaconda.org/anaconda/python>`__. Therefore, use in Windows and Linux at your 
own peril.

As good practice, I recommend installing with anaconda/miniconda and in a new enviroment:

.. parsed-literal:: 
 
    $ conda create -n cfa python=3.6
    $ source activate cfa
    $ conda install -c matplotlib cartopy joblib netcdf4
    $ pip install motionless
    $ pip install git+https://github.com/fmaussion/salem.git
    $ pip install git+https://github.com/raybellwaves/cfanalytics
    $ # eventually: conda install -c conda-forge cfanalytics matplotlib cartopy

You can type ``source deactivate`` when finished. You can also check which environments you have created by typing ``conda info --envs``. 
To remove an environment type ``conda remove --name cfa --all``.

Optional dependencies (fot plotting)
------------------------------------

- `matplotlib <https://github.com/matplotlib/matplotlib>`__
- `cartopy <https://github.com/SciTools/cartopy>`__
- `salem <https://github.com/fmaussion/salem>`__ and its dependencies

Examples
--------

See examples `here <https://github.com/raybellwaves/cfanalytics/tree/master/Examples>`__.

Documentation
-------------

The documentation is hosted at http://cfanalytics.readthedocs.io/

Projects using this data
------------------------

- `Unofficial 2018 Open Results <http://www.rpresidente.com.br/Open2018/Index>`__
- `View a map of the results <https://app.powerbi.com/view?r=eyJrIjoiNmJmODk0MGUtNjVmNi00ZWYxLTg3NjgtOTQ5ZWFlYzFmYjJiIiwidCI6IjQ2YzUxNzhlLWEwZjQtNGY0ZC04YzQwLTk1OThlM2QxMTg2MCIsImMiOjN9>`__

Acknowledgements
----------------

- Thanks to posts on `r/crossfit <https://www.reddit.com/r/crossfit/>`__. e.g. `here <https://www.reddit.com/r/crossfit/comments/5uikq8/2017_open_data_analysis/>`__, I worked out how to download data from the `CrossFit® open <https://games.crossfit.com/leaderboard/open/2017?division=1&region=0&scaled=0&sort=0&occupation=0&page=1>`__. 
- ``Cfopendata`` is a very minor adaptation from `captamericadevs/CFOpenData <https://github.com/captamericadevs/CFOpenData>`__. Who smartly developed code to download CrossFit® open data using `aiohttp <https://github.com/aio-libs/aiohttp>`__. 
- ``Affliatelist`` used this `answer <https://stackoverflow.com/questions/33618324/web-scraping-google-map-website-is-it-possible-to-scrape>`__ as well as this `answer <https://stackoverflow.com/questions/49211863/scrape-latitude-and-longitude-of-address-obtained-from-mapbox>`__ on SO. 
- ``Cfplot().regionplot()`` was made possible by the work of the `cartopy <https://github.com/SciTools/cartopy>`__ developers and developers of `Natural Earth <http://www.naturalearthdata.com/>`__.
- ``Cfplot().cityplot()`` was made possible by the work by `fmaussion <https://github.com/fmaussion>`__ in `salem <https://github.com/fmaussion/salem>`__
- `xarray <https://github.com/pydata/xarray>`__ developers. Whose package template I used for this package as well as the package itself.

Disclaimer
----------

This project is not affiliated with CrossFit, Inc.
