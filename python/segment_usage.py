# -*- coding: utf-8 -*-
"""

=======

@author: katharina

Plot V,J or constant segment family or gene frequencies for H, K or L. Plot as histogram or stacked barchart.
Restrict to different populations/patients by using different sets of event_ids.

usage: segment_usage.py [-h] [-d DATABASE] [-n] [-p {stacked,hist}] [-l {H,K,L}]
                   [-c | -f {V,J} | -g {V,J}]
                   event_infile

positional arguments:
  event_infile      File containing different SQL queries yielding a list of
                    event_ids. Currently the function only supports event_statements which
                    require the pattern match (db, db, db, db).
                    The format resebled fasta format: 
                    
                    +>event_name
                    event_statement

optional arguments:
  -h, --help            show this help message and exit
  
  -d DATABASE, --database DATABASE
                        Optional indication, which database the event statements refer to.
                        If not indicated, database from config file is taken.
                        
  -n, --normalize       Option to activate normalization of data (relative frequencies).
  
  -p {stacked,hist}, --plotstyle {stacked,hist}
                        Choose stacked bargraph or histogram.
                        
  -l {H,K,L}, --locus {H,K,L}
                        Immunoglogulin locus.
                        
  -c, --constant        Show data for constant segments.
  -f {V,J}, --families {V,J}
                        Group data by families of a given segment type V or J.
  -g {V,J}, --genes {V,J}
                        Group data by families of a given segment type V or J.
                        
"""

import bcelldb_init as bcelldb
import igdb_plotting as igplt
import numpy as np
import numpy.random as random
import MySQLdb as mysql
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import colorsys
import argparse
import igdb_queries as igdbq
import seaborn

parser = argparse.ArgumentParser()
parser.add_argument("event_infile", 
                    type = str, 
                    help="File containing different SQL queries yielding a list of event_ids.")
parser.add_argument("-d", "--database", 
                    type = str, 
                    help="optional manual input of database, otherwise taken from config")
parser.add_argument("-n", "--normalize", 
                    help="normalization to total count", 
                    action="store_true")
parser.add_argument("-p", "--plotstyle", type=str, 
                    help="plot stacked or hist", 
                    choices=['stacked','hist'])
parser.add_argument("-l", "--locus", type=str, 
                    help="locus H, K, L", 
                    choices=['H','K','L'])
group = parser.add_mutually_exclusive_group()
group.add_argument("-c", "--constant", 
                   help="resolve constant segments",
                   action="store_true")
group.add_argument("-f", "--families", type=str, 
                   help="resolve families V or J", 
                   choices=['V','J'])
group.add_argument("-g", "--genes", type=str, 
                   help="resolve genes V or J", 
                   choices=['V','J'])

args = parser.parse_args()

# get config variables
conf = bcelldb.get_config()
if args.database:
    db = args.database
else:
    db = conf['dabase']
lib = 'library_scireptor'

if args.constant:
    resolve = 'constant'
    segment = ""
elif args.families:
    resolve = 'families'
    segment = args.families
elif args.genes:
    resolve = 'genes'
    segment = args.genes


# connect to database via ~/.my.conf settings
connection = mysql.connect(db=db,read_default_file="~/.my.cnf", read_default_group='mysql_igdb')
cursor = connection.cursor()

# generate event list. Will later on be generated by another program and taken up by pickle (or called as module).

event_names, event_statements = igdbq.read_eventfile(args.event_infile, db)

def get_gene_list (event_statement):
    if resolve == 'genes':
        
        if args.plotstyle == 'hist':
            order = "ORDER BY cnt DESC;"
        elif args.plotstyle == 'stacked':
            order = "ORDER BY seg_family, seg_gene ASC;"
        
        gene_statement = "SELECT COUNT(seg_family) as cnt, seg_family, seg_gene \
            FROM %s.VDJ_segments \
            JOIN %s.sequences ON sequences.seq_id = VDJ_segments.seq_id \
            AND sequences.consensus_rank = 1 \
            JOIN %s.VDJ_library on VDJ_library.VDJ_id = VDJ_segments.VDJ_id \
            JOIN %s.event ON event.event_id = sequences.event_id \
            where igblast_rank=1 and VDJ_segments.type = '%s' AND VDJ_segments.locus = '%s' \
            AND event.event_id IN (%s) \
            GROUP BY concat(seg_family, seg_gene) " % (db, db, lib, db, segment, args.locus, event_statement)
        gene_statement = gene_statement + order            
        
            
    elif resolve == 'families':

        if args.plotstyle == 'hist':
            order = "ORDER BY cnt DESC;"
        elif args.plotstyle == 'stacked':
            order = "ORDER BY seg_family ASC;"        
        
        gene_statement = "SELECT COUNT(seg_family) as cnt, seg_family \
            FROM %s.VDJ_segments \
            JOIN %s.sequences ON sequences.seq_id = VDJ_segments.seq_id \
            AND sequences.consensus_rank = 1 \
            JOIN %s.VDJ_library on VDJ_library.VDJ_id = VDJ_segments.VDJ_id \
            JOIN %s.event ON event.event_id = sequences.event_id \
            where igblast_rank=1 and VDJ_segments.type = '%s' AND VDJ_segments.locus = '%s' \
            AND event.event_id IN (%s) \
            GROUP BY seg_family " % (db, db, lib, db,segment, args.locus, event_statement)
        gene_statement = gene_statement + order  
            
    elif resolve == 'constant':

        if args.plotstyle == 'hist':
            order = "ORDER BY cnt DESC;"
        elif args.plotstyle == 'stacked':
            order = "ORDER BY constant_segments.name ASC;"         
        
        gene_statement = "SELECT COUNT(constant_segments.name) as cnt, constant_segments.name \
            FROM %s.constant_segments \
            JOIN %s.sequences ON sequences.seq_id = constant_segments.seq_id \
            AND sequences.consensus_rank = 1 \
            JOIN %s.event ON event.event_id = sequences.event_id \
            where sequences.locus = '%s' \
            AND event.event_id IN (%s) \
            GROUP BY constant_segments.name " % (db, db, db, args.locus, event_statement)
        gene_statement = gene_statement + order
            
    else: 
        print "Parameter 'resolve' needs to be either 'genes' or 'families' or 'constant'. Exiting....\n"
        exit

    cursor.execute(gene_statement)
    gene_rows = cursor.fetchall()

    gene_heights = []
    gene_labels = []
    
    for gene in gene_rows:
        gene_heights.append(gene[0])
        if resolve == 'genes':
            gene_labels.append(gene[1] + '-' + gene[2])
        if resolve == 'families':
            gene_labels.append(gene[1])
        if resolve == 'constant':
            gene_labels.append(gene[1])
    
    return gene_heights, gene_labels 
    
if args.plotstyle == 'hist':
    for event_statement, event_name in zip(event_statements, event_names):
        gene_heights, gene_labels = get_gene_list(event_statement)
        plt.figure(figsize=(0.3*len(gene_labels), 7))
        
        if args.normalize == True:
            norm_fact = 1./sum(gene_heights)
            gene_heights = [height*norm_fact for height in gene_heights]
            norm = "norm"
            plt.ylabel("Relative Frequencies")
        else:
            norm_fact = 1
            norm = "abs"
            plt.ylabel("Absolute frequencies")
    
        positions = np.arange(0,len(gene_heights),1)    
        plt.bar(positions, gene_heights, color = 'grey')
        ticks = plt.xticks(positions + 0.4, gene_labels, rotation = 90, fontsize = 12)
        ttl = plt.title(event_name + "\n" + igplt.plot_log('Segment usage', sys.argv))
        plt.savefig("%s_%s_%s_%s_%s_%s_%s" % (args.event_infile, event_name, resolve, segment, args.locus, args.plotstyle, norm) 
                    + '.pdf', bbox_extra_artists=(ttl,), bbox_inches='tight')
        
elif args.plotstyle == 'stacked':
    plt.figure(figsize=(7, 0.7*len(event_statements))) 
    label_list = []
    for event_statement, i in zip(event_statements, range(len(event_statements))):
        #print event_statement
        gene_heights, gene_labels = get_gene_list(event_statement)
        if args.normalize == True:
            norm_fact = 1./sum(gene_heights)
            plt.xlim(0,1)
            norm = "norm"
        else:
            norm_fact = 1
            norm = "abs"
        palette = igplt.random_colors(200)
        bottom = 0
        for height, label in zip(gene_heights, gene_labels):
            # make redundant labels dissapear from legend
            if label not in label_list:
                label_list.append(label)
                leg_label = label
            else: leg_label = ""
            plt.bar(left = bottom, height = 0.8, width = height*norm_fact, orientation='horizontal', label = leg_label, bottom = i, color = igplt.get_color(label) + (0.8,))
            
            bottom = height*norm_fact + bottom
           
        
        # format legend
        if len(label_list) > 15:
            ncol = 2
        else: ncol = 1
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = ncol)
    plt.yticks(np.arange(len(event_names))+0.4, event_names)
    lbl = plt.xlabel('Counts')
    plt.ylim(-0.1, len(event_names)-0.1)
    ttl = plt.title(igplt.plot_log('Segment usage', sys.argv, db))
    plt.savefig("%s_%s_%s_%s_%s_%s_%s" % (db, args.event_infile, resolve, segment, args.locus, args.plotstyle, norm) + '.pdf', bbox_extra_artists=(lgd,ttl,lbl,), bbox_inches='tight')
else:
    print "Plot option must be 'hist' or 'stacked'\n"
    exit