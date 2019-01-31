MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(psimod:PsiMod) 
RETURN psimod.identifier as psimod, psimod.displayName as name, count({pe: pe.displayName, modification: tm.displayName}) as times
ORDER BY times DESC
