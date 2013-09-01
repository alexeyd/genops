#ifndef POPULATION_H
#define POPULATION_H

#include <stdlib.h>
#include <math.h>

#include <functional>
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <memory>


#include "basic_types.h"
#include "differ.h"
#include "generator.h"

namespace GENOPS
{

struct CGenomeProps
{
  size_t m_ttl;

  CGenomeProps()
  : m_ttl(0)
  {
  }
};

typedef std::map<CGenome,CGenomeProps> CGenomeMap;

struct CTournamentFuncs
{
  std::function<size_t()> m_create_tournament;
  std::function<void(size_t,const CGenome&)> m_add_genome;
  std::function<CGenome(size_t)> m_get_result;
};


class CPopulation
{
  public:
    struct CConfig
    {
      size_t m_tournament_size;
      size_t m_operations_per_tick;
      size_t m_genome_ttl;

      size_t m_random_genome_size_ma_period;
      size_t m_random_genome_size_ma_initial_value;
      size_t m_min_population_size;

      double m_max_mutation_percent;
      double m_crossover_similarity;

      CConfig()
      : m_tournament_size(20),
        m_operations_per_tick(100),
        m_genome_ttl(5),
        m_random_genome_size_ma_period(1000),
        m_random_genome_size_ma_initial_value(200),
        m_min_population_size(2000),
        m_max_mutation_percent(0.4),
        m_crossover_similarity(0.7)
      {
      }
    };

  private:
    CConfig m_config;
    CGenerator m_generator;
    CTournamentFuncs m_tournament_funcs;
    CGenomeMap m_genomes;

    size_t m_deleted;
    size_t m_created;
    size_t m_mutated;
    size_t m_crossed;
    double m_tick_similarity;

    size_t m_average_size_cursor;
    std::vector<size_t> m_average_size_array;

    void DeleteOldGenomes();
    void CreateNewGenomes();
    std::vector<size_t> SubmitTournaments();
    void Mutation(const CGenome &genome);
    void UpdateAverageRand(const CGenome &genome);

  public:
    CPopulation(const CConfig &config,
                const CGenePool &gene_pool,
                const CTournamentFuncs &tournament_funcs)
    : m_config(config),
      m_generator(gene_pool),
      m_tournament_funcs(tournament_funcs),
      m_deleted(0),
      m_created(0),
      m_mutated(0),
      m_crossed(0),
      m_tick_similarity(0.0),
      m_average_size_cursor(0),
      m_average_size_array(config.m_random_genome_size_ma_period,
                           config.m_random_genome_size_ma_initial_value)
    {
    }

    void Tick();
    std::vector<CGenome> GetGenomes() const;

    size_t GetDeleted() const { return m_deleted; }
    size_t GetCreated() const { return m_created; }
    size_t GetMutated() const { return m_mutated; }
    size_t GetCrossed() const { return m_crossed; }
    double GetTickSimilarity() const { return m_tick_similarity; }

    size_t GetSize() const { return m_genomes.size(); }

    size_t GetMaxRandSize() const
    {
      size_t max_rand_size = 0;
      for(size_t i = 0; i < m_average_size_array.size(); ++i)
      {
        max_rand_size += m_average_size_array[i];
      }
      max_rand_size /= m_average_size_array.size();
      return max_rand_size;
    }
};

};

#endif

