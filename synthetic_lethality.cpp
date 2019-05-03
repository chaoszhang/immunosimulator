#include <vector>
#include <memory>
#include "synthetic_lethality.hpp"

using namespace std;

uniform_int_distribution<int> uniform_binary(0, 1);
uniform_real_distribution<double> uniform_real(0, 1);

void mutateBase(vector<bool> &seqb, int index){
	bool b1 = seqb[index * 2], b2 = seqb[index * 2 + 1];
	bool nb1 = b1, nb2 = b2;
	while (nb1 == b1 && nb2 == b2){
		nb1 = uniform_binary(generator);
		nb2 = uniform_binary(generator);
	}
	seqb[index * 2] = nb1, seqb[index * 2 + 1] = nb2;
}

struct SyntheticLethalTaxon::Model::Parameters{
	// TODO
};

struct SyntheticLethalTaxon::Model::Impl{
	Impl(const Parameters p, const default_random_engine& g){
		//TODO
	}
	//TODO
};

SyntheticLethalTaxon::Model::Model(const Parameters p): pImpl(new Impl(p, generator)){}

SyntheticLethalTaxon::Gene::Gene(const vector<bool> b): seqbinary(b){}

SyntheticLethalTaxon::Gene SyntheticLethalTaxon::Gene::mutate() const{
	vector<bool> seqb = seqbinary;
	mutateBase(seqb, floor(uniform_real(generator) * seqbinary.size() / 2));
	return SyntheticLethalTaxon::Gene(seqb);
}

SyntheticLethalTaxon::Gene SyntheticLethalTaxon::Gene::mutate(const double mutationRate) const{
	vector<bool> seqb = seqbinary;
	int L = seqbinary.size() / 2;
	for (int i = 0; i < L; i++){
		if (uniform_real(generator) < mutationRate) mutateBase(seqb, i);
	}
	return SyntheticLethalTaxon::Gene(seqb);
}

SyntheticLethalTaxon::Gene SyntheticLethalTaxon::Gene::mutate(bool& mutated, const double mutationRate, const double pressure) const{
	// r' = r / (1 - p)
	vector<bool> seqb = seqbinary;
	int L = seqbinary.size() / 2;
	double q = 1 - mutationRate;
	vector<double> p(L);
	p[L - 1] = q * pressure;
	for (int i = L - 2; i >= 0; i--){
		p[i] = q * p[i + 1];
	}
	for (int i = 0; i < L; i++){
		if (mutated){
			if (uniform_real(generator) < mutationRate) mutateBase(seqb, i);
		}
		else {
			if (uniform_real(generator) * (1 - p[i]) < mutationRate){
				mutateBase(seqb, i);
				mutated = true;
			}
		}
	}
	return SyntheticLethalTaxon::Gene(seqb);
}

void SyntheticLethalTaxon::switchCondition(int condition){
	if (condition == previousCondition) return;
	// TODO
	
	lifeTimeMutationRate = 0;
	for (double r: mutationRates) lifeTimeMutationRate += r;
	double q = 1;
	for (double p: mutationProbabilities) q *= 1 - p;
	birthMutationProbability = 1 - q * q;
}

int main(int argc, char** argv){
	return 0;
}
