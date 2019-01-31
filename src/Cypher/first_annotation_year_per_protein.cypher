// Get first annotation year of each human proteoform
MATCH (re:ReferenceEntity{databaseName:"UniProt"})<-[:referenceEntity]-(pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:literatureReference]->(lr)
WITH DISTINCT re.identifier as id, lr.year as year ORDER BY id, year
WITH id, collect(year) as years
RETURN id, head(years);