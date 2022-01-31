# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 14:59:50 2022

@author: masenpitts
"""

from urllib.request import urlopen

link = "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=2022&MONTH=01&FROM=2412&TO=2412&STNM=72572"

datafile = urlopen(link)
myfile = datafile.read()
header = 3
counter = 0

data = myfile.splitlines()
for eachline in data:
    datalines = eachline.split()
    if len(datalines) == 11:
        if counter >= header:
            print(datalines)
        counter = counter + 1