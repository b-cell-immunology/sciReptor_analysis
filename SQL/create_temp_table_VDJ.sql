# Name:			create_temp_table_VDJ.sql
# Verson:		0.1.1  (2015-03-18)
# Author(s):	Christian Busse
# Maintainer:	Christian Busse (christian.busse@dkfz-heidelberg.de)
# Licence:		AGPL3
# Provides:		Temporary tables containing the connected V(D)J joints for all three loci.
#


# Create table for heavy chain
#
DROP TABLE IF EXISTS temp_table_H_VDJ, temp_table_H_VJ, temp_table_H_D, temp_table_H_C;

CREATE TEMPORARY TABLE temp_table_H_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='H' AND igblast_rank=1) AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='H' AND igblast_rank=1) AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
);

ALTER TABLE temp_table_H_VJ ADD COLUMN d_segment_name VARCHAR(20) AFTER v_segment_name;
ALTER TABLE temp_table_H_VJ ADD COLUMN c_segment_name VARCHAR(20) AFTER j_segment_name;

CREATE TEMPORARY TABLE temp_table_H_D AS (
	SELECT 
		seq_id,
		`name` AS d_segment_name,
		locus
	FROM VDJ_segments
    WHERE type='D'
		AND locus='H'
		AND igblast_rank=1
);

CREATE TEMPORARY TABLE temp_table_H_C AS (
	SELECT 
		seq_id,
		`name` AS c_segment_name
	FROM constant_segments
);

SET SQL_SAFE_UPDATES=0;
UPDATE temp_table_H_VJ, temp_table_H_D SET temp_table_H_VJ.d_segment_name = temp_table_H_D.d_segment_name
WHERE temp_table_H_VJ.seq_id = temp_table_H_D.seq_id;

UPDATE temp_table_H_VJ, temp_table_H_C SET temp_table_H_VJ.c_segment_name = temp_table_H_C.c_segment_name
WHERE temp_table_H_VJ.seq_id = temp_table_H_C.seq_id;
SET SQL_SAFE_UPDATES=1;

DROP TABLE IF EXISTS temp_table_H_D, temp_table_H_C;

ALTER TABLE temp_table_H_VJ RENAME temp_table_H_VDJ;


# Create tables for kappa chain
#
DROP TABLE IF EXISTS temp_table_K_VJ, temp_table_K_C;

CREATE TEMPORARY TABLE temp_table_K_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='K' AND igblast_rank=1) AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='K' AND igblast_rank=1) AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
);

ALTER TABLE temp_table_K_VJ ADD COLUMN c_segment_name VARCHAR(20) AFTER j_segment_name;

CREATE TEMPORARY TABLE temp_table_K_C AS (
	SELECT 
		seq_id,
		`name` AS c_segment_name
	FROM constant_segments
);

SET SQL_SAFE_UPDATES=0;

UPDATE temp_table_K_VJ, temp_table_K_C SET temp_table_K_VJ.c_segment_name = temp_table_K_C.c_segment_name
WHERE temp_table_K_VJ.seq_id = temp_table_K_C.seq_id;

SET SQL_SAFE_UPDATES=1;

DROP TABLE IF EXISTS temp_table_K_C;


# Create tables for lambda chain
#
DROP TABLE IF EXISTS temp_table_L_VJ, temp_table_L_C;

CREATE TEMPORARY TABLE temp_table_L_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='L' AND igblast_rank=1) AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='L' AND igblast_rank=1) AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
);

ALTER TABLE temp_table_L_VJ ADD COLUMN c_segment_name VARCHAR(20) AFTER j_segment_name;

CREATE TEMPORARY TABLE temp_table_L_C AS (
	SELECT 
		seq_id,
		`name` AS c_segment_name
	FROM constant_segments
);

SET SQL_SAFE_UPDATES=0;

UPDATE temp_table_L_VJ, temp_table_L_C SET temp_table_L_VJ.c_segment_name = temp_table_L_C.c_segment_name
WHERE temp_table_L_VJ.seq_id = temp_table_L_C.seq_id;

SET SQL_SAFE_UPDATES=1;

DROP TABLE IF EXISTS temp_table_L_C;
