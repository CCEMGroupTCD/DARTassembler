# Overview of DART Assembler Rework (Jan 2025)

## Introduction
Having recieved the feedback from the reviewers we are going to to attempt to rework the assembler to address the issues raised. 
The main issues raised were:
- lack of novel geometries
- lack of obvious expandability 

## Changes
- we will need to implement a number of new geometries and change the rotation algrorithm to allow for more complex shapes

## Short term Goals
- change the rotation algorithm to use the metal centres 
- incorporate more diverse ligands
- implement a bimettalic geometry

## Shorter term Goals
- implement bindetate, monodentate and haptic pentadentate ligand rotations
- assemble a molecule
- assemble a bimettalic molecule 

## Notes
– We will need to add a filter to target haptic ligands specifically

– Get ligands ---> generate complex ---> generate isomers ---> generate conformers for each isomer ---> choose lowest energy

