MATCH (pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:hasModifiedResidue]->(tm:TranslationalModification)-[:created]-(tmie:InstanceEdit), (psimod:PsiMod)<-[:psiMod]-(tm) 
WHERE psimod.displayName =~ ".*gluc.*"
OR psimod.displayName =~ ".*gly.*"
OR psimod.displayName =~ ".*fucosyl.*"
OR psimod.displayName =~ ".*xylosyl.*"
OR psimod.displayName =~ ".*galactos.*"
RETURN psimod.identifier as psimod, tm.displayName as modification, tmie.dateTime as date, count(pe) as times
ORDER BY date, psimod
