#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <stdlib.h>

#include <set>
#include <vector>
#include <string>
#include <sstream>

namespace GENOPS
{

typedef int gene_t;
typedef std::vector<gene_t> CGenome;
typedef std::set<gene_t> CGenePool;

};

#endif

