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

source_dir = '../..'

def projectpath(*relpaths):
    return os.path.join(source_dir, *relpaths)
