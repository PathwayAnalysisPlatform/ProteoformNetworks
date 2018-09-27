#ifndef PHENOTYPE_OVERLAP_H_
#define PHENOTYPE_OVERLAP_H_

#include "entity.h"

#include <set>
#include <string>

void CreatePhenotypesRatiosFile();

std::set<std::string> GetPhenotypeVertices(Entity entity, std::string& trait);

std::set<std::string> GetPhenotypeOverlap(Entity entity, const std::string& one_trait, const std::string& other_trait);

void CreatePhenotypePairsFile(float min_modified_percentage);

#endif /* PHENOTYPE_OVERLAP_H_ */
