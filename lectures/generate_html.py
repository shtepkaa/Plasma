#!/usr/bin/python
import os

print "Generating the html files"

os.system("ipython nbconvert index.ipynb")
os.system("ipython nbconvert Cheat_sheat.ipynb")
os.system("ipython nbconvert Charged_particles.ipynb")