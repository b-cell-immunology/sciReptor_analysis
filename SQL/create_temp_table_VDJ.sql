# Create temporary tables containing the connected V(D)J joints.
#

# === Heavy ===

DROP TABLE IF EXISTS temp_table_H_VDJ, temp_table_H_VJ, temp_table_H_D, temp_table_H_C;

CREATE TEMPORARY TABLE temp_table_H_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.igblast_rank AS igblast_rank,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='H') AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='H') AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
	AND v_segment.igblast_rank = j_segment.igblast_rank
	AND v_segment.locus = j_segment.locus
);

ALTER TABLE temp_table_H_VJ ADD COLUMN d_segment_name VARCHAR(20) AFTER v_segment_name;
ALTER TABLE temp_table_H_VJ ADD COLUMN c_segment_name VARCHAR(20) AFTER j_segment_name;

CREATE TEMPORARY TABLE temp_table_H_D AS (
	SELECT 
		seq_id,
		`name` AS d_segment_name,
		igblast_rank,
		locus
	FROM VDJ_segments
    WHERE type='D'
	AND locus='H'
);

CREATE TEMPORARY TABLE temp_table_H_C AS (
	SELECT 
		seq_id,
		`name` AS c_segment_name
	FROM constant_segments
);

SET SQL_SAFE_UPDATES=0;
UPDATE temp_table_H_VJ, temp_table_H_D SET temp_table_H_VJ.d_segment_name = temp_table_H_D.d_segment_name
WHERE temp_table_H_VJ.seq_id = temp_table_H_D.seq_id
AND temp_table_H_VJ.igblast_rank = temp_table_H_D.igblast_rank
AND temp_table_H_VJ.locus = temp_table_H_D.locus;

UPDATE temp_table_H_VJ, temp_table_H_C SET temp_table_H_VJ.c_segment_name = temp_table_H_C.c_segment_name
WHERE temp_table_H_VJ.seq_id = temp_table_H_C.seq_id;
SET SQL_SAFE_UPDATES=1;

DROP TABLE IF EXISTS temp_table_H_D, temp_table_H_C;

ALTER TABLE temp_table_H_VJ RENAME temp_table_H_VDJ;

# === Kappa ===

DROP TABLE IF EXISTS temp_table_K_VJ, temp_table_K_C;

CREATE TEMPORARY TABLE temp_table_K_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.igblast_rank AS igblast_rank,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='K') AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='K') AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
	AND v_segment.igblast_rank = j_segment.igblast_rank
	AND v_segment.locus = j_segment.locus
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

# === Lambda ===

DROP TABLE IF EXISTS temp_table_L_VJ, temp_table_L_C;

CREATE TEMPORARY TABLE temp_table_L_VJ AS (
	SELECT 
		v_segment.seq_id AS seq_id,
		v_segment.name AS v_segment_name,
		j_segment.name AS j_segment_name,
		v_segment.igblast_rank AS igblast_rank,
		v_segment.locus AS locus
	FROM 
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='V' AND locus='L') AS v_segment
	INNER JOIN
	(SELECT seq_id, `name`, igblast_rank, locus FROM VDJ_segments WHERE type='J' AND locus='L') AS j_segment
	ON v_segment.seq_id = j_segment.seq_id
	AND v_segment.igblast_rank = j_segment.igblast_rank
	AND v_segment.locus = j_segment.locus
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
