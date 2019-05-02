#ifndef _SIMULATOR_
#define _SIMULATOR_
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <random>

#define NO_EVENT 0
#define REGULAR_BIRTH_EVENT 1
#define BIRTH_MUTATION_EVENT 2
#define LIFE_TIME_MUTATION_EVENT 3
#define DEATH_EVENT 4

using namespace std;

default_random_engine generator;
uniform_real_distribution<double> uniform(0.0, 1.0);
exponential_distribution<double> exponential(1.0);

class Taxon{
public:
	virtual bool operator == (const Taxon &) const {return false;}
	virtual double getOccupancy(int conditon) const {return 0;}
	virtual double getBirthRate(int conditon) const {return 0;}
	virtual double getNaturalDeathRate(int conditon) const {return 0;}
	virtual double getOccupancyDeathRateFactor(int conditon) const {return 0;}
	virtual double getLifeTimeMutationRate(int conditon) const {return 0;}
	virtual double getBirthMutationProbability(int conditon) const {return 0;}
	virtual Taxon* mutateLifeTime(int conditon) const {return new Taxon(*this); }
	virtual pair<Taxon*, Taxon*> mutateBirth(int conditon) const {return pair<Taxon*, Taxon*>(new Taxon(*this), new Taxon(*this)); }
	virtual size_t hash() const {return (size_t) this;}
	virtual ~Taxon(){}
};

class TaxonPointerEq{
public:
	bool operator()(const Taxon* a, const Taxon* b) const { return (*a) == (*b); }
};

class TaxonPointerHash{
public:
	size_t operator()(const Taxon* a) const { return a->hash(); }
};

class Population{
	double occupancy, birthRate, naturalDeathRate, occupancyDeathRateFactor, lifeTimeMutationRate, birthMutationProbability; 
	long long size;
	double timeOfEmergence;
	Population *ancestor = nullptr;
	int nDescendants = 0, condition;
	Taxon* t;
	
public:
	
	Population(Taxon* t, Population *p = nullptr, long long size = 1, int condition = 0, double timeOfEmergence = 0.0): condition(condition), t(t), size(size), timeOfEmergence(timeOfEmergence){
		assignAncestor(p);
		forceUpdateCondition(condition);
	}
	
	void assignAncestor(Population *p){
		if (ancestor != nullptr){
			ancestor->nDescendants--;
		}
		if (p != nullptr){
			p->nDescendants++;
		}
		ancestor = p;
	}
	
	void forceUpdateCondition(int c){
		condition = c;
		occupancy = t->getOccupancy(condition);
		birthRate = t->getBirthRate(condition);
		naturalDeathRate = t->getNaturalDeathRate(condition);
		occupancyDeathRateFactor = t->getOccupancyDeathRateFactor(condition);
		lifeTimeMutationRate = t->getLifeTimeMutationRate(condition);
		birthMutationProbability = t->getBirthMutationProbability(condition);
	}
	
	void updateCondition(int c){
		if (condition == c) return;
		forceUpdateCondition(c);
	}
	
	~Population(){
		assignAncestor(nullptr);
		delete t;
	}
	
	Taxon* getTaxon() const{
		return t;
	}
	
	Population* getAncestor() const{
		return ancestor;
	}
	
	bool hasDescendant() const{
		return nDescendants;
	}
	
	int getBatchSize(double precision) const{
		return floor(precision * size + 1);
	}
	
	long long getSize() const{
		return size;
	}
	
	void adjustSize(long long n){
		size += n;
	}
	
	double getTimeOfEmergence() const{
		return timeOfEmergence;
	}
	
	double getOccupancy() const{
		return occupancy * size;
	}
	
	double getOccupancyIndependentEventRate(double precision) const{
		return ((birthRate * birthMutationProbability + lifeTimeMutationRate)
			+ (birthRate * (1 - birthMutationProbability) + naturalDeathRate) / getBatchSize(precision)) * size;
	}
	
	double getOccupancyEventRateFactor(double precision) const{
		return occupancyDeathRateFactor * size / getBatchSize(precision);
	}
	
	int randomEvent(double totalOccupancy, double precision){
		if (precision == 0){
			double _r = uniform(generator) * (birthRate + lifeTimeMutationRate + naturalDeathRate + occupancyDeathRateFactor * totalOccupancy);
			if (_r < birthRate + lifeTimeMutationRate){
				if (_r < birthRate){
					if (_r < birthRate * birthMutationProbability) return BIRTH_MUTATION_EVENT;
					else return REGULAR_BIRTH_EVENT;
				}
				else return LIFE_TIME_MUTATION_EVENT;
			}
			else return DEATH_EVENT;
		}
		else {
			int n = getBatchSize(precision);
			double _r = uniform(generator) * ((birthRate * birthMutationProbability + lifeTimeMutationRate) * n
				+ (birthRate * (1 - birthMutationProbability) + naturalDeathRate + occupancyDeathRateFactor * totalOccupancy));
			if (_r < (birthRate * birthMutationProbability + lifeTimeMutationRate) * n){
				if (_r < lifeTimeMutationRate * n) return LIFE_TIME_MUTATION_EVENT;
				else return BIRTH_MUTATION_EVENT;
			}
			else{
				double b = birthRate * (1 - birthMutationProbability), d = naturalDeathRate + occupancyDeathRateFactor * totalOccupancy;
				double maxbd = max(b, d);
				b /= maxbd;
				d /= maxbd;
				if (b == d){
					if (uniform(generator) * n < 1){
						if (uniform(generator) < 0.5) return REGULAR_BIRTH_EVENT;
						else return DEATH_EVENT;
					}
					else return NO_EVENT;
				}
				else {
					if (uniform(generator) * abs(pow(b, n) - pow(d, n)) * (b + d) < (pow(b, n) + pow(d, n)) * abs(b - d)){
						if (uniform(generator) * (pow(b, n) + pow(d, n)) < pow(b, n)) return REGULAR_BIRTH_EVENT;
						else return DEATH_EVENT;
					}
					else return NO_EVENT;
				}
			}
		}
	}
};

class PopulationPointerEq{
public:
	bool operator()(const Population* a, const Population* b) const { return *(a->getTaxon()) == *(b->getTaxon()); }
};

class PopulationPointerHash{
public:
	size_t operator()(const Population* a) const { return a->getTaxon()->hash(); }
};

class IntervalTree{
	vector<vector<double> > arr;
	int n = 0;
public:
	IntervalTree(){
		arr.emplace_back(1);
	}
	
	double getSum() const{
		return arr.back()[0];
	}
	
	void clear(){
		for (vector<double> &a: arr){
			for (double &v: a){
				v = 0;
			}
		}
	}
	
	void set(int i, double v){
		arr[0][i] = v;
		for (int j = 1; j < arr.size(); j++){
			i /= 2;
			arr[j][i] = arr[j - 1][2 * i] + arr[j - 1][2 * i + 1];
		}
	}
	
	void append(double v){
		if (n == arr[0].size()){
			for (vector<double> &a: arr){
				a.resize(a.size() * 2);
			}
			arr.emplace_back(1, arr.back()[0]);
		}
		set(n++, v);
	}
	
	int sample() const{
		int i = 0;
		for (int j = arr.size() - 1; j >= 1; j--){
			if (arr[j][i] * uniform(generator) < arr[j - 1][2 * i]) i *= 2;
			else i = 2 * i + 1;
		}
		return i;
	}
};

class TimeCompare{
public:
	bool operator()(const tuple<Population*, Population*, long long> a, const tuple<Population*, Population*, long long> b) const{
		if (get<1>(a)->getTimeOfEmergence() == get<1>(b)->getTimeOfEmergence()) return PopulationPointerEq()(get<0>(b), get<1>(b));
		else return (get<1>(a)->getTimeOfEmergence() < get<1>(b)->getTimeOfEmergence());
	}
};

class Simulator{
	const double precision;
	int nPop = 0, condition = 0;
	double time = 0.0;
	double totalOccupancyIndependentEventRate = 0.0, totalOccupancyEventRateFactor = 0.0, totalOccupancy = 0.0;
	vector<double> occupancy;
	vector<Population*> popList;
	vector<int> emptySpot;
	IntervalTree occupancyIndependentIntervalTree, occupancyFactorIntervalTree;
	unordered_map<Population*, int, PopulationPointerHash, PopulationPointerEq> popPos;
	unordered_map<Taxon*, Population*, TaxonPointerHash, TaxonPointerEq> popSet;
	
	void clear(){
		totalOccupancyIndependentEventRate = totalOccupancyEventRateFactor = totalOccupancy = 0.0;
		for (double &v: occupancy) v = 0;
		occupancyIndependentIntervalTree.clear();
		occupancyFactorIntervalTree.clear();
	}
	
	int allocatePos(Population* pop){
		nPop++;
		if (emptySpot.size()){
			int pos = emptySpot.back();
			emptySpot.pop_back();
			popList[pos] = pop;
			popPos[pop] = pos;
			return pos;
		}
		else {
			int pos = popList.size();
			popList.push_back(pop);
			occupancy.push_back(0);
			occupancyIndependentIntervalTree.append(0);
			occupancyFactorIntervalTree.append(0);
			popPos[pop] = pos;
			return pos;
		}
	}
	
	void freePos(int pos){
		nPop--;
		Population* pop = popList[pos];
		popList[pos] = nullptr;
		popPos.erase(pop);
		emptySpot.push_back(pos);
		garbageCollect(pop);
	}
	
	void garbageCollect(Population* pop){
		if (pop->hasDescendant() == false && pop->getSize() == 0){
			Population* p = pop->getAncestor();
			popSet.erase(pop->getTaxon());
			delete pop;
			if (p != nullptr) garbageCollect(p);
		}
	}
	
	void updateInfo(int pos){
		Population* pop = popList[pos];
		pop->updateCondition(condition);
		setOccupancy(pos, pop->getOccupancy());
		setOccupancyIndependentEventRate(pos, pop->getOccupancyIndependentEventRate(precision));
		setOccupancyEventRateFactor(pos, pop->getOccupancyEventRateFactor(precision));
	}
	
	bool proceedToNextEvent(){
		if (nPop == 0) return false;
		double r = uniform(generator) * (totalOccupancyIndependentEventRate + totalOccupancyEventRateFactor * totalOccupancy);
		int i = (r < totalOccupancyIndependentEventRate) ? occupancyIndependentIntervalTree.sample() : occupancyFactorIntervalTree.sample();
		Population* pop = popList[i];
		while (pop == nullptr){
			cerr << "Simulator::proceedToNextEvent(): pop null pointer exception; retry.\n";
			cerr << "@i = " << i << endl;
			r = uniform(generator) * (totalOccupancyIndependentEventRate + totalOccupancyEventRateFactor * totalOccupancy);
			i = (r < totalOccupancyIndependentEventRate) ? occupancyIndependentIntervalTree.sample() : occupancyFactorIntervalTree.sample();
			pop = popList[i];
		}
		int event = pop->randomEvent(totalOccupancy, precision);
		if (event == REGULAR_BIRTH_EVENT){
			adjustPopSize(pop, pop->getBatchSize(precision));
		}
		if (event == BIRTH_MUTATION_EVENT){
			pair<Taxon*, Taxon*> newT = pop->getTaxon()->mutateBirth(condition);
			createPop(get<0>(newT), pop->getTaxon(), 1, time);
			createPop(get<1>(newT), pop->getTaxon(), 1, time);
			adjustPopSize(pop, -1);
		}
		if (event == LIFE_TIME_MUTATION_EVENT){
			Taxon* newT = pop->getTaxon()->mutateLifeTime(condition);
			createPop(newT, pop->getTaxon(), 1, time);
			adjustPopSize(pop, -1);
		}
		if (event == DEATH_EVENT){
			adjustPopSize(pop, -pop->getBatchSize(precision));
		}
		return true;
	}
	
public:
	
	int getNumPop() const{
		return nPop;
	}
	
	double getTime() const{
		return time;
	}
	
	double getTotalOccupancy() const{
		return totalOccupancy;
	}
	
	void setOccupancy(int i, double v){
		totalOccupancy += v - occupancy[i];
		occupancy[i] = v;
	}
	
	double getTotalOccupancyIndependentEventRate() const{
		return totalOccupancyIndependentEventRate;
	}
	
	void setOccupancyIndependentEventRate(int i, double v){
		occupancyIndependentIntervalTree.set(i, v);
		totalOccupancyIndependentEventRate = occupancyIndependentIntervalTree.getSum();
	}
	
	double getTotalOccupancyEventRateFactor() const{
		return totalOccupancyEventRateFactor;
	}
	
	void setOccupancyEventRateFactor(int i, double v){
		occupancyFactorIntervalTree.set(i, v);
		totalOccupancyEventRateFactor = occupancyFactorIntervalTree.getSum();
	} 
	
	Simulator(vector<tuple<Taxon*, Taxon*, int, double> > startPop, double precision = 0, double time = 0, int condition = 0): condition(condition), precision(precision), time(time){
		for (tuple<Taxon*, Taxon*, int, double> &e: startPop){
			createPop(get<0>(e), nullptr, get<2>(e), get<3>(e));
		}
		for (tuple<Taxon*, Taxon*, int, double> &e: startPop){
			if (get<1>(e) != nullptr) popSet[get<0>(e)]->assignAncestor(popSet[get<1>(e)]);
		}
	}
	
	void createPop(Taxon* t, Taxon* p, long long size, double timeOfEmergence){
		Population* parent = (p != nullptr) ? popSet[p] : nullptr;
		Population* pop;
		if (popSet.count(t)){
			pop = popSet[t];
			delete t;
		} 
		else {
			pop = new Population(t, parent, 0, condition, time);
			popSet[t] = pop;
		}
		adjustPopSize(pop, size);
	}
	
	void adjustPopSize(Population* pop, long long size){
		int pos = (popPos.count(pop)) ? popPos[pop] : allocatePos(pop);
		pop->adjustSize(size);
		updateInfo(pos);
		if (pop->getSize() == 0) freePos(pos);
	}
	
	void simulate(double duration){
		double endTime = time + duration;
		while (nPop != 0){
			double step = exponential(generator) / (totalOccupancyIndependentEventRate + totalOccupancyEventRateFactor * totalOccupancy);
			if (time + step > endTime) break;
			time += step;
			proceedToNextEvent();
		}
		time = endTime;
	}
	
	void switchCondition(int c){
		clear();
		condition = c;
		for (int i = 0; i < popList.size(); i++){
			if (popList[i] == nullptr) continue;
			updateInfo(i);
		}
	}
	
	vector<Population*> getPops() const{
		vector<Population*> result;
		for (Population* p: popList){
			if (p != nullptr) result.push_back(p);
		}
		return result;
	}
	
	vector<tuple<Population*, Population*, long long> > sample(int sampleSize){
		IntervalTree sampleIntervalTree;
		for (Population* p: popList){
			if (p == nullptr) sampleIntervalTree.append(0);
			else sampleIntervalTree.append(p->getSize());
		}
		unordered_map<int, long long> sampleCnt;
		for (int i = 0; i < sampleSize; i++){
			sampleCnt[sampleIntervalTree.sample()]++;
		}
		priority_queue<tuple<Population*, Population*, long long>, vector<tuple<Population*, Population*, long long> >, TimeCompare> q;
		vector<tuple<Population*, Population*, long long> > result;
		for (const pair<int, long long> &e: sampleCnt){
			q.push(make_tuple(popList[e.first], popList[e.first], e.second));
		}
		while (!q.empty()){
			tuple<Population*, Population*, long long> a = q.top();
			q.pop();
			if (!q.empty() && PopulationPointerEq()(get<0>(a), get<0>(q.top()))){
				long long cnt = get<2>(q.top());
				q.pop();
				q.push(make_tuple(get<0>(a), get<1>(a), get<2>(a) + cnt));
				continue;
			}
			else {
				Population* p = get<1>(a)->getAncestor();
				if (p != nullptr) {
					q.push(make_tuple(p, p, 0));
					result.push_back(make_tuple(get<0>(a), p, get<2>(a)));
				}
				else if (get<2>(a) != 0) result.push_back(make_tuple(get<0>(a), get<0>(a), get<2>(a)));
			}
		}
		return result;
	}
	
	vector<tuple<Population*, Population*, long long> > sampleCondensed(int sampleSize){
		long long popCnt = 0;
		IntervalTree sampleIntervalTree;
		for (Population* p: popList){
			if (p == nullptr) sampleIntervalTree.append(0);
			else {
				sampleIntervalTree.append(p->getSize());
				popCnt += p->getSize();
			}
		}
		unordered_map<int, long long> sampleCnt;
		for (int i = 0; i < sampleSize; i++){
			sampleCnt[sampleIntervalTree.sample()]++;
		}
		priority_queue<tuple<Population*, Population*, long long>, vector<tuple<Population*, Population*, long long> >, TimeCompare> q;
		vector<tuple<Population*, Population*, long long> > result;
		for (const pair<int, long long> &e: sampleCnt){
			q.push(make_tuple(popList[e.first], popList[e.first], e.second));
		}
		while (!q.empty()){
			tuple<Population*, Population*, long long> a = q.top();
			q.pop();
			if (!q.empty() && PopulationPointerEq()(get<1>(a), get<1>(q.top()))){
				tuple<Population*, Population*, long long> b = q.top();
				q.pop();
				if (get<0>(a) == get<1>(a)){
					result.push_back(make_tuple(get<0>(b), get<1>(b), get<2>(b)));
					q.push(a);
				}
				else if (get<0>(b) == get<1>(b)){
					result.push_back(make_tuple(get<0>(a), get<1>(a), get<2>(a)));
					q.push(b);
				}
				else {
					result.push_back(make_tuple(get<0>(a), get<1>(a), get<2>(a)));
					result.push_back(make_tuple(get<0>(b), get<1>(b), get<2>(b)));
					q.push(make_tuple(get<1>(a), get<1>(a), 0));
				}
			}
			else {
				Population* p = get<1>(a)->getAncestor();
				if (p != nullptr) q.push(make_tuple(get<0>(a), p, get<2>(a)));
				else {
					if (!PopulationPointerEq()(get<0>(a), get<1>(a))) result.push_back(make_tuple(get<0>(a), get<1>(a), get<2>(a)));
				}
			}
		}
		return result;
	}
	
	void debugInfo() const{
		int popCnt = 0;
		for (Population* p: popList){
			if (p == nullptr) continue;
			popCnt += p->getSize();
		}
		cerr << "Taxon_count: " << nPop << "\tCurrent_time: " << time << "\tTotal_occupancy: " << totalOccupancy << "\tCondition: " << condition << "\tTotal_population: " << popCnt << endl;
	}
	
	void debugInfo(string getInfo(Taxon&)) const{
		debugInfo();
		for (Population* p: popList) {
			if (p == nullptr) continue;
			cerr << "Taxon_info: " << getInfo(*(p->getTaxon())) << endl;
			cerr << "Population_size: " << p->getSize() << "\tTime_Of_Emergence: " << p->getTimeOfEmergence() << endl;
		}
		cerr << endl;
	}
	
	double debugSum(double getValue(Taxon&)) const{
		double result = 0;
		for (Population* p: popList) {
			if (p == nullptr) continue;
			result += getValue(*(p->getTaxon())) * p->getSize();
		}
		return result;
	}
};

class ExampleTaxon: public Taxon{
	double occupancy = 0.0002, birthRate = 1, naturalDeathRate = 0.5, occupancyDeathRateFactor = 1, lifeTimeMutationRate = 0.0, birthMutationProbability = 0.0;
public:
	int type = 0;
	ExampleTaxon(int t){
		if (t == 0){
			type = 0;
			occupancy = 0.00001;
			birthRate = 1;
			naturalDeathRate = 0;
			occupancyDeathRateFactor = 1;
			lifeTimeMutationRate = 0.2;
			birthMutationProbability = 0;
		}
		else {
			type = t;
			occupancy = 0.00001;
			birthRate = 1;
			naturalDeathRate = 0;
			occupancyDeathRateFactor = 1;
			lifeTimeMutationRate = 0;
			birthMutationProbability = 0.1;
		}
	}
	
	Taxon* mutateLifeTime(int conditon) const {return new ExampleTaxon(1 + type); }
	pair<Taxon*, Taxon*> mutateBirth(int conditon) const { return {new ExampleTaxon(1 + type), new ExampleTaxon(1 + type)};}
	
	bool operator == (const Taxon &other) const {return type == ((const ExampleTaxon*) &other)->type;}
	double getOccupancy(int condition) const {return occupancy * (1 + 99 * condition);}
	double getBirthRate(int condition) const {return birthRate;}
	double getNaturalDeathRate(int condition) const {return naturalDeathRate;}
	double getOccupancyDeathRateFactor(int condition) const {return occupancyDeathRateFactor;}
	double getLifeTimeMutationRate(int condition) const {return lifeTimeMutationRate;}
	double getBirthMutationProbability(int condition) const {return birthMutationProbability;}
	size_t hash() const {return type;}
};
#endif
