select distinct t1.source_organism_org_id, count(t1.source_organism_org_id)
as count from
(
select distinct structure_id, source_organism_org_id
from structure_object where source_organism_org_id is not null
) t1
group by source_organism_org_id;