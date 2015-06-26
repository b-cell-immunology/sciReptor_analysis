# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/katharina/.spyder2/.temp.py
"""

import numpy as np
#import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import MySQLdb as mysql
import bcelldb_init as bcelldb
from datetime import datetime as dt
import itertools as itt
import igdb_queries as igdbq
import igdb_plotting as igdbplt
from mpl_toolkits.axes_grid1 import ImageGrid

# get configuration using bcelldb_init
conf = bcelldb.get_config()


event_infile = 'healthy_donors.events'

# database
db = 'healthy'
experiment_id = 'D01'
locus = 'H'
plate_barcode = 'D01X0002C0'


# connect to database via ~/.my.conf settings
dab = mysql.connect(db=db,read_default_file="~/.my.cnf", read_default_group='mysql_igdb')
cursor = dab.cursor()

# read in event file
event_names, event_statements = igdbq.read_eventfile(event_infile, db)

def get_null_events (event_statement):
    # get all the events from one experiment
    null_events_statement = event_statement
    cursor.execute(null_events_statement)
    null_events = cursor.fetchall()
    return null_events

def get_positive_events (event_statement):
    # get the event_ids where sequences where amplified
    events_statement = "SELECT heavy_light.event_id FROM heavy_light \
        WHERE event_id in (%s) " % (event_statement)   
    cursor.execute(events_statement)
    events = cursor.fetchall()
    return events

def get_channels (plate_barcode):
    channel_statement = "SELECT marker_name FROM flow \
        JOIN flow_meta \
        ON flow.channel_id = flow_meta.channel_id \
        JOIN event ON event.event_id = flow.event_id \
        WHERE plate_barcode = '%s' and marker_name != 'None' GROUP BY marker_name;" % (plate_barcode)
    cursor.execute(channel_statement)
    channels = cursor.fetchall()
    return channels

channels = get_channels(plate_barcode)
# channels koennen auch ueber kommandozeile uebergeben werden

for combi in itt.combinations(channels, 2):
    # INITIATE PLOTTING INSTANCE
    F = plt.figure(1,(9.5, 5.5))
    
    grid = ImageGrid(F, 111,
              nrows_ncols = (1, len(event_names)),
              direction="row",
              axes_pad = 0.05,
              add_all=True,
              share_all = True,
              )
    channel1, channel2 = combi
    grid[0].set_ylabel(channel2[0])
    
    for event_name, event_statement, ax in zip(event_names, event_statements, grid):
        
        null_events = get_null_events(event_statement)
        
        value_channel_dict = {}
        null_value_channel_dict = {}
        
        # plot all events (grey)    
        for event in null_events:
            for channel in [channel1[0],channel2[0]]:
                try:
                    query_statement = "SELECT value FROM flow \
                    JOIN flow_meta \
                    ON flow.channel_id = flow_meta.channel_id \
                    WHERE event_id = %d and marker_name='%s';" % (int(event[0]),channel)
                    cursor.execute(query_statement)
                    value = cursor.fetchall()
                    value = value[0][0]
                    print value
                    
                    try:
                        null_value_channel_dict[channel].append(value)
                    except KeyError:
                        null_value_channel_dict[channel] = []
                        null_value_channel_dict[channel].append(value)
        
                except IndexError:
                    next
        
        events = get_positive_events(event_statement)        
        
        for event in events:
            for channel in [channel1[0],channel2[0]]:
                try:
                    query_statement = "SELECT value FROM flow \
                    JOIN flow_meta \
                    ON flow.channel_id = flow_meta.channel_id \
                    WHERE event_id = %d and marker_name='%s';" % (int(event[0]),channel)
                    cursor.execute(query_statement)
                    value = cursor.fetchall()
                    value = value[0][0]
                    print value
                    
                    try:
                        value_channel_dict[channel].append(value)
                    except KeyError:
                        value_channel_dict[channel] = []
                        value_channel_dict[channel].append(value)
                        
                    # determine corresponding isotype
                    if channel == channel1[0]:
                        isotype = igdbq.get_H_isotype(int(event[0]), cursor)
                        if isotype:
                            isotype = isotype[0][0]
                            try:
                                color = igdbplt.get_color(isotype)
                                # some are IGKC???
                            except KeyError:
                                    color = 'black'
                                    print event
                                    print isotype
                        else: color = 'white'
                            
                        try:
                            value_channel_dict['color'].append(color)
                            print 'appended color'
                        except KeyError:
                            value_channel_dict['color'] = []
                            value_channel_dict['color'].append(color)
                            print 'key error'
                    
                
                except IndexError:
                    next
        
        for x,y in zip(null_value_channel_dict[channel1[0]], null_value_channel_dict[channel2[0]]):
            if (x>=0 and y>=0):
    		ax.scatter(np.arcsinh(x), np.arcsinh(y), color = 'lightgrey', s = 80)
            #plt.scatter(x, y, color = 'lightgrey', s = 80)
        
        
        
        #plt.figure()
        for x,y,c in zip(value_channel_dict[channel1[0]], value_channel_dict[channel2[0]], value_channel_dict['color']):
            if (x>=0 and y>=0):	
    		ax.scatter(np.arcsinh(x), np.arcsinh(y), color = c, s = 20)
            #plt.scatter(x, y, color = c, s = 20)
            #plt.scatter(x, y, color = c)
        ax.set_xlabel(channel1[0] + "\n" + event_name)
        ticks = [0,10, 100, 10**3,10**4,10**5]
        tick_labels = ["0","1E+01","1E+02","1E+03", "1E+04", "1E+05"]
        ax.set_xticks(np.arcsinh(ticks))
        ax.set_xticklabels(tick_labels, size = 7, rotation = 90)
        ax.set_yticks(np.arcsinh(ticks))
        ax.set_yticklabels(tick_labels, size = 7)
        ax.set_xlim(0,np.arcsinh(10**5))
        ax.set_ylim(0,np.arcsinh(10**5))
        plt.show()
    
    plt.savefig(experiment_id+'cluster_'+channel1[0]+'_'+channel2[0]+'.pdf')
        
    plt.close()
    
#        for key,color in sorted(isotype_color_dict.items()):
#            plt.figtext(0.1,position_count,key, color = color, size=14, weight='bold', style = 'italic')
#            position_count = position_count + 0.04
#        plt.xlabel(channel1, size = 14)
#        plt.ylabel(channel2, size = 14)    
    #ighm = [color == 'darkorchid' for color in value_channel_dict['color']]
    #for seqevent, igmbool in zip(events,ighm):
     #   if igmbool is True:
     #       isotype_statement = "SELECT name FROM " + database + ".constant_segments \
     #                   where seq_id = %d;" % (int(seqevent[0]))
     #       cursor.execute(isotype_statement)
     #       print cursor.fetchall()
     #       print seqevent
