SELECT DISTINCT t1.source_organism_org_id, COUNT(t1.source_organism_org_id) AS count FROM
	(SELECT DISTINCT structure_id, source_organism_org_id
	 FROM structure_object WHERE source_organism_org_id IS NOT NULL) t1
GROUP BY source_organism_org_id;