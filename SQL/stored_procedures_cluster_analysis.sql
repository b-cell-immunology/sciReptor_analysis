DELIMITER $$
CREATE PROCEDURE `create_cluster`()
BEGIN

DROP TABLE IF EXISTS derived_cluster_meta, derived_cluster_comp, derived_cluster;
DROP TEMPORARY TABLE IF EXISTS temp_table_cluster_comp;

# DROP FUNCTION IF EXISTS LEN_CDR3;
# DELIMITER //
# CREATE FUNCTION LEN_CDR3(cdr3 varchar(100))
# RETURNS int(10) unsigned
# BEGIN
#   DECLARE cdr3_len int(10) unsigned;
#   IF cdr3 IS NULL THEN SET cdr3_len = NULL;
#   ELSE SET cdr3_len = LENGTH(cdr3);
#   END IF;
#   RETURN cdr3_len;
# END; //
# DELIMITER ;

CREATE TABLE derived_cluster_meta (
	`cluster_comp_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
	`nodes` int(10) unsigned NOT NULL,
	`donor_identifier` varchar(45) NOT NULL,
	`locus` char(1) NOT NULL,
	`segment_v` varchar(20),
	`segment_j` varchar(20),
	`cdr3_length` int(10) unsigned,
	UNIQUE KEY `cluster_comp_id_UNIQUE` (`cluster_comp_id`)
) AS
SELECT
	COUNT(DISTINCT event_id) AS nodes,
	donor_identifier,
	CAST('H' AS char(1)) AS locus,
	igh_segment_v AS segment_v,
	igh_segment_j AS segment_j,
	LENGTH(igh_cdr3) AS cdr3_length
FROM derived_segment_association
GROUP BY donor_identifier, segment_v, segment_j, cdr3_length
UNION
SELECT
	COUNT(DISTINCT event_id) AS nodes,
	donor_identifier,
	CAST('K' AS char(1)) AS locus,
	igk_segment_v AS segment_v,
	igk_segment_j AS segment_j,
	LENGTH(igk_cdr3) AS cdr3_length
FROM derived_segment_association
GROUP BY donor_identifier, segment_v, segment_j, cdr3_length
UNION
SELECT
	COUNT(DISTINCT event_id) AS nodes,
	donor_identifier,
	CAST('L' AS char(1)) AS locus,
	igl_segment_v AS segment_v,
	igl_segment_j AS segment_j,
	LENGTH(igl_cdr3) AS cdr3_length
FROM derived_segment_association
GROUP BY donor_identifier, segment_v, segment_j, cdr3_length;


CREATE TEMPORARY TABLE temp_table_cluster_comp (
	`event_id` int(10) unsigned NOT NULL,
	`cluster_comp_heavy_id` int(10) unsigned NOT NULL,
	`cluster_comp_kappa_id` int(10) unsigned,
	`cluster_comp_lambda_id` int(10) unsigned
) AS
SELECT
	segview.event_id,
	cluster_comp_heavy.cluster_comp_id AS cluster_comp_heavy_id,
	cluster_comp_kappa.cluster_comp_id AS cluster_comp_kappa_id,
	cluster_comp_lambda.cluster_comp_id AS cluster_comp_lambda_id
FROM (
	SELECT
		*,
		LENGTH(igh_cdr3) AS igh_cdr3_length,
		LENGTH(igk_cdr3) AS igk_cdr3_length,
		LENGTH(igl_cdr3) AS igl_cdr3_length
	FROM derived_segment_association
) AS segview
INNER JOIN (
	SELECT * FROM derived_cluster_meta WHERE locus='H'
) AS cluster_comp_heavy
ON segview.donor_identifier = cluster_comp_heavy.donor_identifier
	AND segview.igh_segment_v = cluster_comp_heavy.segment_v
	AND segview.igh_segment_j = cluster_comp_heavy.segment_j
	AND segview.igh_cdr3_length = cluster_comp_heavy.cdr3_length
LEFT OUTER JOIN (
	SELECT * FROM derived_cluster_meta WHERE locus='K'
) AS cluster_comp_kappa
ON segview.donor_identifier = cluster_comp_kappa.donor_identifier
	AND segview.igk_segment_v = cluster_comp_kappa.segment_v
	AND segview.igk_segment_j = cluster_comp_kappa.segment_j
	AND segview.igk_cdr3_length = cluster_comp_kappa.cdr3_length
LEFT OUTER JOIN (
	SELECT * FROM derived_cluster_meta WHERE locus='L'
) AS cluster_comp_lambda
ON segview.donor_identifier = cluster_comp_lambda.donor_identifier
	AND segview.igl_segment_v = cluster_comp_lambda.segment_v
	AND segview.igl_segment_j = cluster_comp_lambda.segment_j
	AND segview.igl_cdr3_length = cluster_comp_lambda.cdr3_length
;

CREATE TABLE derived_cluster_comp (
	`cluster_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
	`nodes` int(10) unsigned NOT NULL,
	`cluster_comp_heavy_id` int(10) unsigned NOT NULL,
	`cluster_comp_kappa_id` int(10) unsigned,
	`cluster_comp_lambda_id` int(10) unsigned,
	UNIQUE KEY `cluster_id_UNIQUE` (`cluster_id`)
) AS
SELECT
	COUNT(DISTINCT event_id) AS nodes,
	cluster_comp_heavy_id,
	cluster_comp_kappa_id,
	cluster_comp_lambda_id
FROM temp_table_cluster_comp
GROUP BY cluster_comp_heavy_id, cluster_comp_kappa_id, cluster_comp_lambda_id;

CREATE TABLE derived_cluster (
	`event_id` int(10) unsigned NOT NULL,
	`cluster_id` int(10) unsigned NOT NULL
) AS
SELECT
	temp_table_cluster_comp.event_id AS event_id,
	derived_cluster_comp.cluster_id AS cluster_id
FROM temp_table_cluster_comp
INNER JOIN derived_cluster_comp
ON temp_table_cluster_comp.cluster_comp_heavy_id = derived_cluster_comp.cluster_comp_heavy_id
AND (
	temp_table_cluster_comp.cluster_comp_kappa_id = derived_cluster_comp.cluster_comp_kappa_id
	OR (
		temp_table_cluster_comp.cluster_comp_kappa_id IS NULL
		AND derived_cluster_comp.cluster_comp_kappa_id IS NULL
	)
) AND (
	temp_table_cluster_comp.cluster_comp_lambda_id = derived_cluster_comp.cluster_comp_lambda_id
	OR (
		temp_table_cluster_comp.cluster_comp_lambda_id IS NULL
		AND derived_cluster_comp.cluster_comp_lambda_id IS NULL
	)
);


END$$

DELIMITER ;
