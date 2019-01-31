MATCH (pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:hasModifiedResidue]->(tm:TranslationalModification)-[:created]-(tmie:InstanceEdit), (psimod:PsiMod)<-[:psiMod]-(tm) 
WHERE psimod.displayName =~ ".*phospho.*" 
RETURN pe.displayName as physicalEntity, psimod.identifier as psimod, tm.displayName as modification, tmie.dateTime as date 
ORDER BY date, psimod
