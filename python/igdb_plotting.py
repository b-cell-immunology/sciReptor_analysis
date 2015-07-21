# -*- coding: utf-8 -*-
"""
bcelldb_plotting.py

Module that contains standard funtions for plotting data from the igdb.
"""

import bcelldb_init as bcelldb
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

#for line in color_file:
#    if (line != "\n" and line[0] != '#'):
#        try :
#            entries = line[:-1].split('\t')
#            color_dict[entries[0]] = (int(entries[1]), int(entries[2]), int(entries[3]))
#        except IndexError:
#            pass
 
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16)/255 for i in range(0, lv, lv // 3))





#hls_list = []  
#plt.figure()
#for i,color in zip(range(len(color_list)), color_list):
#    rgb = hex_to_rgb(color)
#    hls = colorsys.rgb_to_hls(rgb[0], rgb[1], rgb[2])
#    #print rgb
#    #print hls
#    hls_list.append(hls)
#    rgb2 = colorsys.hls_to_rgb(hls[0],hls[1],hls[2])
#    #plt.scatter(i,i,color=rgb, s= 300)    
#    #plt.scatter(i,i/2,color=rgb2, s= 300)  
#    
#hls_list.sort()
#
#ext_hls_list = []
#for hls in hls_list:
#    ext_hls_list.append(hls)
#    hue = hls[0]
#    li = hls[1]
#    sat = hls[2]
#    for s in np.arange(0.3,1,0.2):
#        ext_hls_list.append((hue,li,s))
#    for l in np.arange(0.2,0.8,0.15):
#        ext_hls_list.append((hue,l,sat))
#        
#
#for i,hls in zip(range(len(ext_hls_list)), ext_hls_list):
#    rgb2 = colorsys.hls_to_rgb(hls[0],hls[1],hls[2])  
#    #plt.scatter(i*2,i/2,color=rgb2, s= 300, label = hls)  
#
#out = open('/home/katharina/Desktop/hls_list.txt','w')
#for line in ext_hls_list:
#    out.write(str(line)+ '\n')
#out.close()


# N random colors
def random_colors (N):
    colors = []
    for n in range(N):
        rgb = tuple(random.random_sample(3))
        colors.append(rgb)
    return colors
