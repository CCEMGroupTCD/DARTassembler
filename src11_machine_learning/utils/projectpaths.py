#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:18:34 2021

@author: Timo Sommer

Contains the path to the source directory of this paper.
"""
import os 
import inspect
import src11_machine_learning
from pathlib import Path

source_dir = '/Users/timosommer/PhD/projects/RCA/projects/DART/'

def projectpath(*relpaths):
    return Path(source_dir, *relpaths)
