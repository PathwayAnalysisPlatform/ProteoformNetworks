#ifndef PATHWAY_OVERLAP_H_
#define PATHWAY_OVERLAP_H_

#include "entity.h"

#include <set>
#include <string>

void CreatePathwayRatiosFile();

std::set<std::string> GetPathwayVertices(Entity entity, std::string& trait);

std::set<std::string> GetPathwayOverlap(Entity entity, const std::string& one_trait, const std::string& other_trait);

void CreatePathwayPairsFile(float min_modified_percentage);

#endif /* PATHWAY_OVERLAP_H_ */
