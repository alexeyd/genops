#include "generator.h"

namespace GENOPS
{

CGenome CGenerator::Create(size_t size) const
{
  CGenome output;

  for(size_t i = 0; i < size; ++i)
  {
    output.push_back( m_genes[rand() % m_genes.size()] );
  }

  return output;
}


CGenome CGenerator::AddMutation(const CGenome &input) const
{
  CGenome output = input;
  gene_t new_gene = m_genes[rand() % m_genes.size()];
  size_t pos = rand() % (output.size() + 1);

  CGenome::iterator iter = output.begin();

  while(pos)
  {
    ++iter;
    --pos;
  }

  output.insert(iter, new_gene);
  return output;
}


CGenome CGenerator::DeleteMutation(const CGenome &input) const
{
  CGenome output = input;
  if(output.empty())
  {
    return output;
  }

  size_t pos = rand() % output.size();

  CGenome::iterator iter = output.begin();

  while(pos)
  {
    ++iter;
    --pos;
  }

  output.erase(iter);
  return output;
}


CGenome CGenerator::ReplaceMutation(const CGenome &input) const
{
  CGenome output = input;
  if(output.empty())
  {
    return output;
  }

  gene_t new_gene = m_genes[rand() % m_genes.size()];
  size_t pos = rand() % output.size();

  CGenome::iterator iter = output.begin();

  while(pos)
  {
    ++iter;
    --pos;
  }

  (*iter) = new_gene;
  return output;
}


CGenome CGenerator::Crossover(const CGenome &father, 
                              const CGenome &mother) const
{
  std::vector<CGenomeDiff> diffs = CDiffer::GenerateDiffs(father, mother);

  if(diffs.empty())
  {
    return father;
  }

  CDiffer differ(father);

  for(size_t i = 0; i < diffs.size(); ++i)
  {
    if(rand() % 2 == 1)
    {
      differ.Apply(diffs[i]);
    }
  }

  return differ.Assemble();
}

};

