# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import MySQLdb as mysql
from datetime import datetime as dt
import itertools as itt
import igdb_queries as igdbq
import igdb_plotting as igdbplt
from mpl_toolkits.axes_grid1 import ImageGrid
import argparse

# parse arguments

parser = argparse.ArgumentParser()

parser.add_argument("event_infile", 
                    type = str, 
                    help="File containing different SQL queries yielding a list of event_ids.")
parser.add_argument("-d", "--database", 
                    type = str, 
                    help="manual input of database scheme")
parser.add_argument("-l", "--locus", type=str, 
                    help="locus H, K, L", 
                    choices=['H','K','L'])
# need to find out which channels to select. Take them from exemplary plate with barcode
parser.add_argument("-pb", "--platebarcode", type=str, 
                    help="plate barcode to select channels that will be displayed")
parser.add_argument("-o", "--outputdir", type=str, 
                    help="directory for pdf output") 
parser.add_argument("-m", "--mutation", 
                   help="Size according to mutation count",
                   action="store_true")

args = parser.parse_args()

db = args.database
locus = args.locus
plate_barcode = args.platebarcode
event_infile = args.event_infile

# connect to database via ~/.my.conf settings
dab = mysql.connect(db=db,read_default_file="~/.my.cnf", read_default_group='mysql_igdb')
cursor = dab.cursor()

# generate heavy_light table
igdbq.create_temp_heavy_light(cursor)

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
    
def arcsinh_fct (x):
    x = np.array(x)
    y = np.arcsinh(x/10.)
    return y

channels = get_channels(plate_barcode)

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
                        else: color = 'white'
                            
                        try:
                            value_channel_dict['color'].append(color)
                        except KeyError:
                            value_channel_dict['color'] = []
                            value_channel_dict['color'].append(color)
                        
                        factor = 1
                        if args.mutation == True:
                            factor = igdbq.get_mutation_count(int(event[0]), cursor)
                            try:
                                value_channel_dict['size'].append(factor)
                            except KeyError:
                                value_channel_dict['size'] = []
                                value_channel_dict['size'].append(factor)
                    
                
                except IndexError:
                    next
        
        for x,y in zip(null_value_channel_dict[channel1[0]], null_value_channel_dict[channel2[0]]):
            #if (x>=0 and y>=0):
    		ax.scatter(arcsinh_fct(x), arcsinh_fct(y), color = 'lightgrey', s = 70, alpha=0.3)
        
        
        
        for x,y,c,s in zip(value_channel_dict[channel1[0]], value_channel_dict[channel2[0]], value_channel_dict['color'], value_channel_dict['size']):
            #if (x>=0 and y>=0):	
            ax.scatter(arcsinh_fct(x), arcsinh_fct(y), color = c, s = s/len(event_names))
        ax.set_xlabel(channel1[0] + "\n" + event_name)
        ticks = [-100, 0,10, 100, 10**3,10**4,10**5]
        tick_labels = ["-1E+02","0","1E+01","1E+02","1E+03", "1E+04", "1E+05"]
        ax.set_xticks(arcsinh_fct(ticks))
        ax.set_xticklabels(tick_labels, size = 7, rotation = 90)
        ax.set_yticks(arcsinh_fct(ticks))
        ax.set_yticklabels(tick_labels, size = 7)
        ax.set_xlim(arcsinh_fct(-10**2),arcsinh_fct(10**5))
        ax.set_ylim(arcsinh_fct(-10**2),arcsinh_fct(10**5))
        plt.show()

    plt.tight_layout()    
    plt.savefig(args.outputdir + '/flow_'+args.event_infile[:-7] + '_' +channel1[0]+'_'+channel2[0]+'.pdf')
        
    plt.close()

# drop temporary heavy_light table
igdbq.drop_temp_heavy_light(cursor)

