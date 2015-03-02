# -*- coding: utf-8 -*-
"""
Module for common igdb database queries, functions for standardized reading of db info etc....

@author: katharina
"""

def read_eventfile (eventfile, db):
    """
    Args:
    eventfile   File containing event names and event queries in the following format:
                +>event_name1
                event_query1
                +>event_name1
                event_query1
                Event query is completed by %(db,db,db,db)                
    db          Database.
    
    Returns:
    event_names         List of names given to each set of events.
    event_statements    List of database query command for each set of events.
    
    """
    infile = open(eventfile, 'r').readlines()
    event_names = []
    event_statements = []
    event_statement = ''

    for line in infile:
        if line[:2] == '+>':
            event_names.append(line[2:-1])
            # update the event statement list, when new event comes in, but only from second identifier on
            if not len(event_names) == 1:
                event_statement = event_statement % (db,db,db,db)
                event_statements.append(event_statement)
                event_statement = ''
        elif line != '\n':
            event_statement = event_statement + line
    # add last event statement
    event_statement = event_statement % (db,db,db,db)
    event_statements.append(event_statement)
    
    return event_names, event_statements