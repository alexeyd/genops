#include "differ.h"

namespace GENOPS
{

double CDiffer::CalcSimilarity(const CGenome &a, const CGenome &b)
{
  dtl::Diff<gene_t> genome_diff(a, b);
  genome_diff.compose();

  CGenome lcs_genome = genome_diff.getLcsVec();

  size_t total_size = a.size();
  if(b.size() > total_size)
  {
    total_size = b.size();
  }

  return static_cast<double>(lcs_genome.size()) /
         static_cast<double>(total_size);
}


std::vector<CGenomeDiff> CDiffer::GenerateDiffs(const CGenome &a,
                                                const CGenome &b)
{
  std::vector<CGenomeDiff> diffs;

  dtl::Diff<gene_t> genome_diff(a, b);
  genome_diff.compose();

  const dtl::Ses<gene_t> &ses = genome_diff.getSes();
  const dtl::Ses<gene_t>::sesElemVec &ses_elems = ses.getSequence();

  size_t pos = 0;
  bool diff_active = false;
  CGenomeDiff diff;

  for(size_t i = 0; i < ses_elems.size(); ++i)
  {
    switch(ses_elems[i].second.type)
    {
      case dtl::SES_DELETE:
        if(!diff_active)
        {
          diff_active = true;
          diff.m_pos = pos;
        }

        ++(diff.m_remove_size);
        ++pos;
        break;

      case dtl::SES_ADD:
        if(!diff_active)
        {
          diff_active = true;
          diff.m_pos = pos;
        }

        diff.m_add_genes.push_back(ses_elems[i].first);
        break;

      case dtl::SES_COMMON:
        if(diff_active)
        {
          diff_active = false;
          diffs.push_back(diff);

          diff.m_pos = 0;
          diff.m_remove_size = 0;
          diff.m_add_genes.clear();
        }
        ++pos;
        break;

      default:
      {
        std::ostringstream error_stream;
        error_stream<<__FILE__<<":"<<__LINE__<<" "
                    <<"Unknown SES type: "
                    <<ses_elems[i].second.type;

        throw std::runtime_error(error_stream.str());
      }
    }
  }

  if(diff_active)
  {
    diffs.push_back(diff);
  }

  return diffs;
}


CDiffer::CDiffer(const CGenome &original)
{
  for(size_t i = 0; i < original.size(); ++i)
  {
    m_mapped_genome[i].push_back(original[i]);
  }
}


void CDiffer::Apply(const CGenomeDiff &diff)
{
  if(diff.m_remove_size)
  {
    for(size_t i = diff.m_pos; i < diff.m_pos + diff.m_remove_size; ++i)
    {
      m_mapped_genome[i].pop_front();
    }
  }

  if(!diff.m_add_genes.empty())
  {
    CGenes &genes = m_mapped_genome[diff.m_pos];
    genes.insert(genes.begin(), diff.m_add_genes.begin(),
                 diff.m_add_genes.end());
  }
}


CGenome CDiffer::Assemble() const
{
  CGenome output;

  for(CMappedGenome::const_iterator genes = m_mapped_genome.begin();
      genes != m_mapped_genome.end(); ++genes)
  {
    for(CGenes::const_iterator gene = genes->second.begin();
        gene != genes->second.end(); ++gene)
    {
      output.push_back(*gene);
    }
  }

  return output;
}

};

