// Get creation date of each modification that have a PSIMOD ontology id
MATCH (psimod:PsiMod)<-[:psiMod]-(tm:TranslationalModification)-[:created]-(tmie:InstanceEdit)
RETURN DISTINCT psimod.displayName as psimod, tm.displayName as modification, tmie.dateTime as concreteModDate 
ORDER BY concreteModDate, psimod