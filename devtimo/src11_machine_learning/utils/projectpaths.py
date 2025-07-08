#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:18:34 2021

@author: Timo Sommer

Contains the path to the source directory of this paper.
"""
import os 
import inspect
from pathlib import Path

source_dir = '//'

def projectpath(*relpaths):
    return Path(source_dir, *relpaths)
