# -*- coding: utf-8 -*-
"""
bcelldb_plotting.py

Module that contains standard funtions for plotting data from the igdb.
"""

import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random
import sys
import colorsys
from datetime import datetime as dt

conf = bcelldb.get_config()

# Get log information for a plot
# date, user, database, software version
# config file variables

def plot_log (object_name, args, db, short='T'):
    
    """
    Arguments:
        - args: sys.args 
        - verbose
    Return string for logging how plot was generated.
    Short: For titles in PDF.
    Long: also logging script that was executed to generate the plot
    """
    
    date = dt.now().strftime('%Y-%m-%d %H:%M:%S')

    log_string = "%s generated from database %s on %s.\nCommand: \"%s\"\n" % (object_name, db, date, " ".join(args))
    
    if (short == 'F'):
        
        file_name = args[0]
        file_content = open(filename, 'r').read()
        
        log_string = log_string + file_content
    
    return log_string



# Define colors

colortxt = 'seg_families_colors.txt'
color_file = open(colortxt, 'r').readlines()

color_dict = {}

for line in color_file:
    if (line != "\n" and line[0] != '#'):
        entries = line[:-1].split(',')
        try:
            color_dict[entries[0]] = colorsys.hls_to_rgb(float(entries[1]), float(entries[2]), float(entries[3])) 
        except ValueError:
            next

def get_color (name):
    try:
        return color_dict[name]
    except KeyError:
        try:
            print "WARNING COLOR CODE"
            print name + " does not have a color definition in " + colortxt
            print "Assigning random color to " + name
            color_dict[name] = random_colors(1)[0]
            return color_dict[name]
        except TypeError:
            print 'None Type for one of the coloridentifiers. Choosing black color.'
            return 'black'

 
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16)/255 for i in range(0, lv, lv // 3))

# N random colors
def random_colors (N):
    colors = []
    for n in range(N):
        rgb = tuple(random.random_sample(3))
        colors.append(rgb)
    return colors
