# Dev of python tool for RMEDtoolbox

* only modules in environment.yml can be used
* follow MetOffice lint stands using run_pylint.py
* expected to be used via rose suites and occuationally as command line tool

## Coming soon Example RMED scripts

<hr>

# Code dev

* Data is: PP format

<hr>

# RMED environment

Normally ran with in RMED suite:

A worked example [here](https://code.metoffice.gov.uk/trac/rmed/wiki/suites/eval/vn0.0_worked_eg).

Specifically: Example 4: Produce diagnostic plots (area-averages, histograms and cell-statistics)

`rosie copy u-at789@81389` shows how it interacts with suite.
Currently this only works on monsoon.

## Requirements:

* MONSOON ACCESS (for interaction with rose suite...)
* JASMIN ACCESS for worked examples
