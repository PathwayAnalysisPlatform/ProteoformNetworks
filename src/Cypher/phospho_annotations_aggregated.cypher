MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:hasModifiedResidue]->(tm:TranslationalModification)-[:created]-(tmie:InstanceEdit), (psimod:PsiMod)<-[:psiMod]-(tm) 
WHERE psimod.displayName =~ '.*phospho.*' 
RETURN psimod.identifier as psimod, tm.displayName as modification, tmie.dateTime as date, count(pe) as times
ORDER BY date, psimod
