# Name			:	Stored Procedures for segment analysis
# Author		:	Christian Busse
# Maintainer	:	Christian Busse (christian.busse@dkfz-heidelberg.de)
# Version:		:	0.1.1
# Date			:	2015-05-02
# License		:	AGPLv3
# Description	:	This script generates several stored procedures to assist in the analysis of segment associations.
#					Currently these are 'sub_temp_table_segment_association', which combines V, D and segments of all
#					given consensus sequences and 'create_segment_view' which further combines this data to create a
#					table containing associated/paired sequences based on identical event_id. Importantly,
#					'create_segments_view' does consider the functionality status of the sequences and will fall back
#					to the secondary sequence of a given locus if the primary is non-functional.
# Notes			:	ATTENTION: 1) 'create_segment_view' does its functionality assessment *ONLY* based on the presence
#					of stop codons in the CDR3. 2) Runtime of the procedure seems to scale exponentially with the amount
#					of data. While small data sets complete within minutes, large data sets can take several hours.
#					3) Be aware that 'sequences.consensus_rank' is NULL for Sanger sequences (this is handled correctly
#					in the current implementation).
#

DELIMITER $$
CREATE PROCEDURE `sub_temp_table_segment_association`()
BEGIN
# Create temporary tables containing the connected V(D)J joints.
#

# === Heavy ===

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

# === Kappa ===

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

# === Lambda ===

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
END$$


CREATE PROCEDURE `create_segment_association`()
BEGIN
# Generate a plain table showing the associated segments
# Requires: sub_temp_table_association function
#
CALL sub_temp_table_segment_association;
CREATE TABLE derived_segment_association ENGINE=MYISAM AS (
	SELECT
		`event`.event_id,
		`event`.plate,
		`event`.well,
		donor.donor_identifier,
		sample.tissue,
		sort.population,
		sort.antigen,
		temp_associated.igh_segment_v,
		temp_associated.igh_segment_d,
		temp_associated.igh_segment_j,
		temp_associated.igh_cdr3,
		temp_associated.igh_segment_c,
		temp_associated.igk_segment_v,
		temp_associated.igk_segment_j,
		temp_associated.igk_cdr3,
		temp_associated.igk_segment_c,
		temp_associated.igl_segment_v,
		temp_associated.igl_segment_j,
		temp_associated.igl_cdr3,
		temp_associated.igl_segment_c
	FROM (
		SELECT 
			temp_heavy.event_id			AS igh_segment_event_id,
			temp_heavy.v_segment_name	AS igh_segment_v,
			temp_heavy.d_segment_name	AS igh_segment_d,
			temp_heavy.j_segment_name	AS igh_segment_j,
			temp_heavy.CDR3				AS igh_cdr3,
			temp_heavy.c_segment_name	AS igh_segment_c,
			temp_kappa.v_segment_name	AS igk_segment_v,
			temp_kappa.j_segment_name	AS igk_segment_j,
			temp_kappa.CDR3				AS igk_cdr3,
			temp_kappa.c_segment_name	AS igk_segment_c,
			temp_lambda.v_segment_name	AS igl_segment_v,
			temp_lambda.j_segment_name	AS igl_segment_j,
			temp_lambda.CDR3			AS igl_cdr3,
			temp_lambda.c_segment_name	AS igl_segment_c
		FROM (
			SELECT 
				sequences.event_id,
				temp_table_H_VDJ.v_segment_name,
				temp_table_H_VDJ.d_segment_name,
				temp_table_H_VDJ.j_segment_name,
				CDR.prot_seq AS CDR3,
				temp_table_H_VDJ.c_segment_name
			FROM temp_table_H_VDJ
			INNER JOIN (
				# The following query basically returns a seq_id and a CDR3 protein sequence, however it also contains
				# an additional logic that the sequence with the lower consensus_rank will be given preferences unless it
				# does contain a stop codon in CDR3. This is implemented via a LEFT JOIN which is one way to select for
				# group-wise extrema in MySQL (comparable to the OVER ... PARTION BY statement in PostgreSQL).
				#
				SELECT
					temp_table_functionality_1.seq_id,
					temp_table_functionality_1.prot_seq
				FROM (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='H'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_1
				LEFT JOIN (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='H'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_2
				ON temp_table_functionality_1.event_id = temp_table_functionality_2.event_id 
					AND temp_table_functionality_1.consensus_rank > temp_table_functionality_2.consensus_rank
				WHERE temp_table_functionality_2.event_id IS NULL
			) AS CDR
			INNER JOIN sequences
			ON temp_table_H_VDJ.seq_id = CDR.seq_id
				AND temp_table_H_VDJ.seq_id = sequences.seq_id
		) AS temp_heavy
		LEFT OUTER JOIN (
			SELECT 
				sequences.event_id,
				temp_table_K_VJ.v_segment_name,
				temp_table_K_VJ.j_segment_name,
				CDR.prot_seq AS CDR3,
				temp_table_K_VJ.c_segment_name
			FROM temp_table_K_VJ
			INNER JOIN (
				SELECT
					temp_table_functionality_1.seq_id,
					temp_table_functionality_1.prot_seq
				FROM (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='K'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_1
				LEFT JOIN (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='K'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_2
				ON temp_table_functionality_1.event_id = temp_table_functionality_2.event_id 
					AND temp_table_functionality_1.consensus_rank > temp_table_functionality_2.consensus_rank
				WHERE temp_table_functionality_2.event_id IS NULL
			) AS CDR
			INNER JOIN sequences
			ON temp_table_K_VJ.seq_id = CDR.seq_id
				AND temp_table_K_VJ.seq_id = sequences.seq_id
		) AS temp_kappa
		ON temp_heavy.event_id = temp_kappa.event_id
		LEFT OUTER JOIN (
			SELECT 
				sequences.event_id,
				temp_table_L_VJ.v_segment_name,
				temp_table_L_VJ.j_segment_name,
				CDR.prot_seq AS CDR3,
				temp_table_L_VJ.c_segment_name
			FROM temp_table_L_VJ
			INNER JOIN (
				SELECT
					temp_table_functionality_1.seq_id,
					temp_table_functionality_1.prot_seq
				FROM (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='L'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_1
				LEFT JOIN (
					SELECT 
						sequences.event_id,
						sequences.seq_id,
						sequences.consensus_rank,
						CDR_FWR.prot_seq,
						CDR_FWR.stop_codon 
					FROM sequences
					INNER JOIN CDR_FWR
					ON sequences.seq_id=CDR_FWR.seq_id
					WHERE locus='L'
						AND region='CDR3'
						AND stop_codon=0
				) AS temp_table_functionality_2
				ON temp_table_functionality_1.event_id = temp_table_functionality_2.event_id 
					AND temp_table_functionality_1.consensus_rank > temp_table_functionality_2.consensus_rank
				WHERE temp_table_functionality_2.event_id IS NULL
			) AS CDR
			INNER JOIN sequences
			ON temp_table_L_VJ.seq_id = CDR.seq_id
				AND temp_table_L_VJ.seq_id = sequences.seq_id
		) AS temp_lambda
		ON temp_heavy.event_id = temp_lambda.event_id
		WHERE temp_kappa.v_segment_name IS NOT NULL 
			OR temp_lambda.v_segment_name IS NOT NULL
	) AS temp_associated 
	INNER JOIN `event`
	INNER JOIN sort
	INNER JOIN sample
	INNER JOIN donor
	ON `event`.event_id = temp_associated.igh_segment_event_id
		AND sample.donor_id = donor.donor_id
		AND `event`.sort_id = sort.sort_id
		AND sort.sample_id = sample.sample_id
);

END$$

DELIMITER ;
