#ifndef _SYNTHETIC_LETHALITY_
#define _SYNTHETIC_LETHALITY_
#include <vector>
#include <memory>
#include "simulator.hpp"

using namespace std;

class SyntheticLethalTaxon: public Taxon{
	struct Gene{
		Gene(const string);
		Gene(const vector<bool>);
		double similarity(const Gene, const Gene, const vector<double> condonweight) const; // aa similarity
		/*
			= SUM BLOSUM62((3kmer in Gene1 -> aa), (3kmer in Gene2 -> aa)) * condonweight
			special case 1st aa != M -> similarity = -OO
						last aa != Stop -> Similarity = -OO
						other than last aa == Stop -> Similarity = -OO
		*/
		Gene mutate() const; // mutate one base
		Gene mutate(const double mutationRate) const;
		Gene mutate(bool& mutated, const double mutationRate, const double pressure) const;
		
		const vector<bool> seqbinary;
		/*
		00 -> A; 01 -> C; 10 -> G; 11 -> T
		00011011 -> ACGT
		*/
	};
	
	struct Model{
		struct Parameters;
		struct Impl;
		
		Model(const Parameters);
		const shared_ptr<Gene> geneTarget(int condition, int geneId) const;
		/*
		Random: requirments:
			NOW ASSUME condition = 0
			same condition && same geneid -> same target
			diff condition && same geneid -> similar target (TBD)
			1st aa = M
			last aa = stop
			other than last != stop
		*/
		const vector<double> condonweight(int condition, int geneId) const;
		/*
		
		*/
		const vector<double> mutationRates(int condition, const vector<double> similarities) const;
		const vector<double> mutationProbabilities(int condition, const vector<double> similarities) const;
		double occupany(int condition, const vector<double> similarities) const;
		double birthRate(int condition, const vector<double> similarities) const;
		double naturalDeathRate(int condition, const vector<double> similarities) const;
		double occupancyDeathRateFactor(int condition, const vector<double> similarities) const;
		
		unique_ptr<Impl> pImpl;
	};
	
	
public:
	SyntheticLethalTaxon(const vector<shared_ptr<Gene> >, Model&, int condition);
	~SyntheticLethalTaxon();
	Taxon* mutateLifeTime(int conditon) const;
	pair<Taxon*, Taxon*> mutateBirth(int conditon) const;
	
	bool operator == (const Taxon &other) const;
	double getOccupancy(int condition) const;
	double getBirthRate(int condition) const;
	double getNaturalDeathRate(int condition) const;
	double getOccupancyDeathRateFactor(int condition) const;
	double getLifeTimeMutationRate(int condition) const;
	double getBirthMutationProbability(int condition) const;
	size_t hash() const;
	void switchCondition(int condition);
	
private:
	const vector<shared_ptr<Gene> > pGenes;
	const Model& model;
	int previousCondition;
	vector<double> mutationRates, mutationProbabilities;
	double occupancy, birthRate, naturalDeathRate, occupancyDeathRateFactor, lifeTimeMutationRate, birthMutationProbability;
};

#endif
