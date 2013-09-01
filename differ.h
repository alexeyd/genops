#ifndef DIFFER_H
#define DIFFER_H

#include <map>
#include <list>
#include <stdexcept>

#include "dtl/dtl.hpp"
#include "basic_types.h"

namespace GENOPS
{


struct CGenomeDiff
{
  size_t m_pos;
  size_t m_remove_size;
  CGenome m_add_genes;

  CGenomeDiff()
  : m_pos(0), m_remove_size(0)
  {
  }
};


class CDiffer
{
  private:
    typedef std::list<gene_t> CGenes;
    typedef std::map<size_t, CGenes> CMappedGenome;
    CMappedGenome m_mapped_genome;

  public:
    static double CalcSimilarity(const CGenome &a, const CGenome &b);
    static std::vector<CGenomeDiff> GenerateDiffs(const CGenome &a,
                                                  const CGenome &b);


    explicit CDiffer(const CGenome &original);

    void Apply(const CGenomeDiff &diff);
    CGenome Assemble() const;
};

};

#endif

