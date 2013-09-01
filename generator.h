#ifndef GENERATOR_H
#define GENERATOR_H

#include <stdlib.h>

#include <vector>
#include <set>
#include <algorithm>

#include <iostream>

#include "basic_types.h"
#include "differ.h"

namespace GENOPS
{

class CGenerator
{
  private:
    std::vector<gene_t> m_genes;

  public:
    CGenerator(const CGenePool &gene_pool)
    : m_genes(gene_pool.begin(), gene_pool.end())
    {
    }

    CGenome Create(size_t size) const;

    CGenome AddMutation(const CGenome &input) const;
    CGenome DeleteMutation(const CGenome &input) const;
    CGenome ReplaceMutation(const CGenome &input) const;

    CGenome Crossover(const CGenome &father, 
                      const CGenome &mother) const;
};

};

#endif

