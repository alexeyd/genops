#include "population.h"

namespace GENOPS
{

void CPopulation::DeleteOldGenomes()
{
  CGenomeMap::iterator genome_iter = m_genomes.begin();
  while(genome_iter != m_genomes.end())
  {
    if(genome_iter->second.m_ttl == 0)
    {
      genome_iter = m_genomes.erase(genome_iter);
      ++m_deleted;
    }
    else
    {
      --(genome_iter->second.m_ttl);
      ++genome_iter;
    }
  }
}


void CPopulation::CreateNewGenomes()
{
  size_t max_rand_size = GetMaxRandSize();

  while(m_genomes.size() < m_config.m_min_population_size)
  {
    size_t genome_size = (rand() % max_rand_size) + 1;
    CGenome new_genome = m_generator.Create(genome_size);
 
    if(m_genomes.find(new_genome) == m_genomes.end())
    {
      ++m_created;
      m_genomes[new_genome].m_ttl = m_config.m_genome_ttl;
    }
  }
}


std::vector<size_t> CPopulation::SubmitTournaments()
{
  std::vector<size_t> tournament_ids;
  for(size_t i = 0; i < m_config.m_operations_per_tick; ++i)
  {
    size_t collected = 0;
    size_t cursor = 0;
    size_t desired_cursor = 0;

    size_t tournament_id = m_tournament_funcs.m_create_tournament();
    std::set<size_t> select_set;

    CGenomeMap::const_iterator iter = m_genomes.begin();

    while(collected < m_config.m_tournament_size)
    {
      do
      {
        desired_cursor = rand() % m_genomes.size();
      }
      while(select_set.insert(desired_cursor).second == false);

      while(desired_cursor < cursor)
      {
        --iter;
        --cursor;
      }

      while(desired_cursor > cursor)
      {
        ++iter;
        ++cursor;
      }

      m_tournament_funcs.m_add_genome(tournament_id, iter->first);
      ++collected;
    }

    tournament_ids.push_back(tournament_id);
  }

  return tournament_ids;
}


void CPopulation::Mutation(const CGenome &genome)
{
  size_t mutation_count = 1;

  if(!genome.empty())
  {
    mutation_count = (rand() % genome.size()) + 1;
  }

  CGenome mutant = genome;

  for(size_t i = 0; i < mutation_count; ++i)
  {
    if(CDiffer::CalcSimilarity(genome, mutant) <= 
         (1.0 - m_config.m_max_mutation_percent))
    {
      break;
    }

    int mutation_id = rand() % 3;
    if(mutation_id == 0)
    {
      mutant = m_generator.AddMutation(mutant);
    }
    else if(mutation_id == 1)
    {
      mutant = m_generator.DeleteMutation(mutant);
    }
    else
    {
      mutant = m_generator.ReplaceMutation(mutant);
    }
  }


  if(m_genomes.find(mutant) == m_genomes.end())
  {
    m_genomes[mutant].m_ttl = m_config.m_genome_ttl;
    ++m_mutated;
  }
}


void CPopulation::UpdateAverageRand(const CGenome &genome)
{
  m_average_size_array[m_average_size_cursor] = genome.size() * 2;
  ++m_average_size_cursor;
  if(m_average_size_cursor >= m_average_size_array.size())
  {
    m_average_size_cursor = 0;
  }
}

struct CCrossCand
{
  const CGenome *m_father;
  const CGenome *m_mother;
  double m_similarity;

  bool operator < (const CCrossCand &other) const
  {
    return m_similarity > other.m_similarity;
  }
};


void CPopulation::Tick()
{
  m_deleted = 0;
  m_created = 0;
  m_mutated = 0;
  m_crossed = 0;

  DeleteOldGenomes();
  CreateNewGenomes();

  std::vector<size_t> tournament_ids = SubmitTournaments();
  std::vector<CGenome> leaders;

  for(size_t i = 0; i < m_config.m_operations_per_tick; ++i)
  {
    leaders.push_back(m_tournament_funcs.m_get_result(tournament_ids[i]));
    m_genomes[leaders.back()].m_ttl = m_config.m_genome_ttl;
    UpdateAverageRand(leaders.back());
    Mutation(leaders.back());
  }

  m_tick_similarity = 0.0;
  size_t denom = 0;

  std::vector<CCrossCand> cross_cands;

  for(size_t i = 0; i < leaders.size(); ++i)
  {
    const CGenome &father = leaders[i];
    for(size_t j = i + 1; j < leaders.size(); ++j)
    {
      const CGenome &mother = leaders[j];

      double similarity = CDiffer::CalcSimilarity(father, mother);
      m_tick_similarity += similarity;
      ++denom;

      if(similarity >= m_config.m_crossover_similarity)
      {
        CCrossCand cross_cand;
        cross_cand.m_father = &father;
        cross_cand.m_mother = &mother;
        cross_cand.m_similarity = similarity;
        cross_cands.push_back(cross_cand);
      }
    }
  }

  std::sort(cross_cands.begin(), cross_cands.end());
  
  for(size_t i = 0; i < cross_cands.size(); ++i)
  {
    CGenome child = m_generator.Crossover(*(cross_cands[i].m_father),
                                          *(cross_cands[i].m_mother));

    if(m_genomes.find(child) == m_genomes.end())
    {
      m_genomes[child].m_ttl = m_config.m_genome_ttl;
      ++m_crossed;

      if(m_crossed >= m_config.m_operations_per_tick)
      {
        break;
      }
    }
  }

  m_tick_similarity /= static_cast<double>(denom);
}

std::vector<CGenome> CPopulation::GetGenomes() const
{
  std::vector<CGenome> result;
  for(CGenomeMap::const_iterator genome_iter = m_genomes.begin();
      genome_iter != m_genomes.end(); ++genome_iter)
  {
    result.push_back(genome_iter->first);
  }

  return result;
}

};

