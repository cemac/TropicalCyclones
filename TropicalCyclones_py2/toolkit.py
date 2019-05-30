# -*- coding: utf-8 -*-
"""Toolkit

.. module:: Toolkit
    :platform: Unix
    :synopis: Standalone tools

.. moduleauthors: John Ashcroft, Helen Burns CEMAC (UoL)

.. date: May 2019.

.. description: This module was developed by CEMAC as part of the WCSSP
                Project. Shared funtions accross all plotting tools

   :copyright: © 2019 University of Leeds.
   :license: MIT.

Example:
    To use::

Memebers:

.. CEMAC_TropicalCyclones:
   https://github.com/cemac/TropicalCyclones
"""
import iris


def annotate(axs, str_format, xy):
    """annotate
    Description:
        Create cube of WindSpeeds for all ensemble members
    Args:
        axs (fig axes): figure axes
        str_format (str): Regex string
        xy:
    Return:
        Adds annoation to axs
    """
    # Add initial time, valid time etc.
    bbox_args = dict(boxstyle="round", fc="0.8")
    axs.annotate(str_format, xy=xy,
                 xycoords='figure fraction', xytext=(40, 20),
                 textcoords='offset points', ha="right", va="top",
                 bbox=bbox_args, fontsize=16)


def box_constraint(minlat, maxlat, minlon, maxlon):
    """box_constraint
    Description:
        Contrain cube using lat lons
    Args:
        minlat (int): minium latitue
        minlon (int): minimum longitude
        maxlat (int): maximum latitue
        maxlon (int): maximum longitude
    Return:
        tc_box_constraint (iris Constraints)
    """
    # Create constraint to extract data from cube over a certain region
    longitude_constraint1 = iris.Constraint(
        longitude=lambda cell: cell > minlon)
    longitude_constraint2 = iris.Constraint(
        longitude=lambda cell: cell < maxlon)
    latitude_constraint1 = iris.Constraint(latitude=lambda cell: cell > minlat)
    latitude_constraint2 = iris.Constraint(latitude=lambda cell: cell < maxlat)
    tc_box_constraint = (longitude_constraint1 & longitude_constraint2 &
                         latitude_constraint1 & latitude_constraint2)
    return tc_box_constraint


def extracter(fload, minlon, maxlon, minlat, maxlat):
    """extracter
    Description:
    Args:
        fload (iris cube): loaded file
        minlon (int): minimum longitude
        maxlon (int): maximum longitude
        minlat (int): minimum latitude
        maxlat (int): maximum latitude
    Return:
        ir (iris cube): Contrained iris cube
    """
    cube = fload.extract(iris.Constraint(longitude=lambda cell: minlon < cell <
                                         maxlon, latitude=lambda cell: minlat
                                         < cell < maxlat))
    return cube
