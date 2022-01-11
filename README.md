# PanfluteTippingPoints

## Introduction
This repository contains data for the manuscript: Yuya Karita, David T. Limmer, Oskar Hallatschek, "Scale-dependent tipping points of bacterial colonization resistance", bioRxiv, 2021.
Contact yuya.karita[at]berkeley.edu if you have any questions.

## Data descriptions
Most of the data were taken by Olympus IX81 microscope with a 10x objective every 20 minutes, and 1 um = 1.55 pixel, unless otherwise specified.
### Magnification
Data whose file name contains "4x" or "20x" were taken by a 4x or 20x objective respectively.
### Microscope
For confocal images, Zeiss LSM 700 was used. For "DiffusivityMeasurement/Self_jammed" and "TempChange/30to22", Nikon Eclipse Ti was used.
### Timelapse
"JammingDynamics" was taken every 10 minutes, "PIV" was taken every 3 minutes, and "DiffusivityMeasurement/Self_gaseous" was taken every 30 seconds.

## Codes descriptions
"growth_p.c" is a C code for doing the agent based Brownian dynamics simulation with input file "input.txt" 
Compiled with gcc -o growth_p growth_p.c -O3 and ran with ./growth_p
