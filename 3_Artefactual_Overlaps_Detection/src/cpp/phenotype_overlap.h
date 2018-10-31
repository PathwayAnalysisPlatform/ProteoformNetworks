#ifndef PHENOTYPE_OVERLAP_H_
#define PHENOTYPE_OVERLAP_H_

#include "entity.h"

#include <map>
#include <set>
#include <string>

bool IsModified(const std::string& proteoform);

std::set<std::string> GetTraits();

void CreatePhenotypeRatiosFile();

std::set<std::string> GetPhenotypeVertices(Entity entity, const std::string& trait);

std::set<std::string> GetPhenotypeVerticesFromEdgesFile(Entity entity, const std::string& trait);

std::set<std::string> GetVertexOverlap(Entity entity, const std::string& one_trait, const std::string& other_trait);

std::set<std::string> GetVertexOverlap(Entity entity, const std::string& one_trait, const std::string& other_trait, bool verbose);

void CreatePhenotypePairsFile(double min_modified_percentage, Entity entity_type);

void CheckNumVerticesForTraits();

std::multimap<std::string, std::string> GetPairsOverlappingAtGenes(std::set<std::string> genes, std::multimap<std::string, std::string> candidate_pairs);

std::set<std::string> GetProteoformsWithAccession(const std::string& protein_accession, std::set<std::string> candidate_proteoforms);

void ShowOverlapsAndSplittingProteoforms(std::multimap<std::string, std::string> pairs);

#endif /* PHENOTYPE_OVERLAP_H_ */
