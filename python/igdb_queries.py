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
    


def create_temp_heavy_light (db_cursor):
    create_statement = "CREATE TABLE IF NOT EXISTS heavy_light AS ( \
        SELECT sequences.event_id, sequences.seq_id as H_seq_id, sequences.locus as H_locus, \
        max_light.seq_id as KL_seq_id, max_light.locus as KL_locus \
        FROM \
        (SELECT s1.event_id, s1.locus, s1.seq_id FROM \
        (SELECT sequences.event_id, sequences.locus, sequences.seq_id, n_seq FROM sequences \
        JOIN consensus_stats \
        ON consensus_stats.sequences_seq_id = sequences.seq_id \
        AND (sequences.locus = 'K' OR sequences.locus = 'L') \
        AND sequences.consensus_rank = 1) as s1 \
        LEFT JOIN \
        (SELECT sequences.event_id, sequences.locus, sequences.seq_id, n_seq FROM sequences \
        JOIN consensus_stats \
        ON consensus_stats.sequences_seq_id = sequences.seq_id \
        AND (sequences.locus = 'K' OR sequences.locus = 'L') \
        AND sequences.consensus_rank = 1) as s2 \
        ON s1.event_id = s2.event_id \
        AND s1.n_seq < s2.n_seq \
        WHERE s2.event_id is NULL \
        ) as max_light \
        JOIN sequences ON sequences.event_id = max_light.event_id \
        AND sequences.locus = 'H' AND sequences.consensus_rank = 1 \
        );"
    db_cursor.execute(create_statement)
    
def drop_temp_heavy_light (db_cursor):
    drop_statement = "DROP TABLE heavy_light;"
    db_cursor.execute(drop_statement)

def get_H_isotype (event_id, db_cursor):
    isotype_statement = "SELECT constant_segments.name FROM constant_segments \
        JOIN sequences ON sequences.seq_id = constant_segments.seq_id \
        WHERE locus = 'H' AND consensus_rank = 1 AND event_id = %d;" % (event_id)
    db_cursor.execute(isotype_statement)
    isotype = db_cursor.fetchall()
    return isotype

def get_mutation_count (event_id, cursor):
    mutation_statement = "SELECT sum(replacement) + sum(silent) \
        FROM mutations\
        JOIN sequences ON sequences.seq_id = mutations.seq_id AND sequences.consensus_rank = 1 \
        WHERE event_id = %d;" % (event_id)
    cursor.execute(mutation_statement)
    try:
        mutation = cursor.fetchall()[0][0]
    except IndexError:
        mutation = 0
    return int(mutation)

def get_mutation_count_seqid (seq_id, cursor):
    mutation_statement = "SELECT sum(replacement) + sum(silent) \
        FROM mutations \
        JOIN sequences ON sequences.seq_id = mutations.seq_id \
        WHERE sequences.seq_id = %d;" % (seq_id)
    cursor.execute(mutation_statement)
    try:
        mutation = cursor.fetchall()[0][0]
        return int(mutation)
    except TypeError:
        return 0
    