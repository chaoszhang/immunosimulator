#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include "simulator.hpp"

using namespace std;

const int BLOSUM[23][23] = {
{8, -3, -4, -5, -2, -2, -3, -1, -4, -4, -4, -2, -3, -5, -2, 1, -1, -6, -5, -2, -4, -2, -2}, 
{-3, 10, -2, -5, -8, 0, -2, -6, -1, -7, -6, 3, -4, -6, -5, -3, -3, -7, -5, -6, -4, -1, -3}, 
{-4, -2, 11, 1, -5, -1, -2, -2, 0, -7, -7, -1, -5, -7, -5, 0, -1, -8, -5, -7, 5, -2, -3}, 
{-5, -5, 1, 10, -8, -2, 2, -4, -3, -8, -8, -3, -8, -8, -5, -2, -4, -10, -7, -8, 6, 0, -4}, 
{-2, -8, -5, -8, 14, -7, -9, -7, -8, -3, -5, -8, -4, -4, -8, -3, -3, -7, -6, -3, -7, -8, -5}, 
{-2, 0, -1, -2, -7, 11, 2, -5, 1, -6, -5, 2, -2, -6, -4, -2, -3, -5, -4, -5, -2, 5, -2}, 
{-3, -2, -2, 2, -9, 2, 10, -6, -2, -7, -7, 0, -5, -8, -4, -2, -3, -8, -7, -5, 0, 7, -3}, 
{-1, -6, -2, -4, -7, -5, -6, 9, -6, -9, -8, -5, -7, -8, -6, -2, -5, -7, -8, -8, -3, -5, -4}, 
{-4, -1, 0, -3, -8, 1, -2, -6, 13, -7, -6, -3, -5, -4, -5, -3, -4, -5, 1, -7, -2, -1, -4}, 
{-4, -7, -7, -8, -3, -6, -7, -9, -7, 8, 2, -6, 1, -2, -7, -5, -3, -6, -4, 4, -8, -7, -3}, 
{-4, -6, -7, -8, -5, -5, -7, -8, -6, 2, 8, -6, 3, 0, -7, -6, -4, -5, -4, 0, -8, -6, -3}, 
{-2, 3, -1, -3, -8, 2, 0, -5, -3, -6, -6, 10, -4, -6, -3, -2, -3, -8, -5, -5, -2, 0, -3}, 
{-3, -4, -5, -8, -4, -2, -5, -7, -5, 1, 3, -4, 12, -1, -5, -4, -2, -4, -5, 0, -7, -4, -3}, 
{-5, -6, -7, -8, -4, -6, -8, -8, -4, -2, 0, -6, -1, 11, -7, -5, -5, 0, 4, -3, -7, -7, -4}, 
{-2, -5, -5, -5, -8, -4, -4, -6, -5, -7, -7, -3, -5, -7, 12, -3, -4, -8, -7, -6, -5, -4, -4}, 
{1, -3, 0, -2, -3, -2, -2, -2, -3, -5, -6, -2, -4, -5, -3, 9, 2, -7, -5, -4, -1, -2, -2}, 
{-1, -3, -1, -4, -3, -3, -3, -5, -4, -3, -4, -3, -2, -5, -4, 2, 9, -7, -5, -1, -2, -3, -2}, 
{-6, -7, -8, -10, -7, -5, -8, -7, -5, -6, -5, -8, -4, 0, -8, -7, -7, 17, 2, -5, -9, -7, -6}, 
{-5, -5, -5, -7, -6, -4, -7, -8, 1, -4, -4, -5, -5, 4, -7, -5, -5, 2, 12, -5, -6, -6, -4}, 
{-2, -6, -7, -8, -3, -5, -5, -8, -7, 4, 0, -5, 0, -3, -6, -4, -1, -5, -5, 8, -7, -5, -3}, 
{-4, -4, 5, 6, -7, -2, 0, -3, -2, -8, -8, -2, -7, -7, -5, -1, -2, -9, -6, -7, 6, 0, -4}, 
{-2, -1, -2, 0, -8, 5, 7, -5, -1, -7, -6, 0, -4, -7, -4, -2, -3, -7, -6, -5, 0, 6, -2}, 
{-2, -3, -3, -4, -5, -2, -3, -4, -4, -3, -3, -3, -3, -4, -4, -2, -2, -6, -4, -3, -4, -2, -3}
};

const unordered_map<char, int> AA2INT({{'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4}, {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13}, {'P', 14}, {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19}, {'B', 20}, {'Z', 21}, {'X', 22}});
const unordered_map<string, char> CODON({
{"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'}, 
{"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
{"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'}, 
{"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
{"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'}, 
{"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
{"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'}, 
{"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
{"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'}, 
{"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
{"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'}, 
{"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
{"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'}, 
{"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
{"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'}, 
{"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
});

struct Arguments{
	double NORMAL_BIRTH_RATE = 0;
	double NONSENSE_DEATH_RATE = 1000;
	double NORMAL_NONMEMORY_DEATH_RATE = 1000;
	double NORMAL_DEATH_RATE = 1.0 / 402;
	double NORMAL_DEATH_FACTOR = 0.04 / 100000;
	double INFACTED_BIRTH_RATE = 4;
	double INFACTED_DEATH_RATE = 1;
	double INFACTED_BASE_DEATH_FACTOR = 0.000003;
	double AFFINITY_FACTOR = 0.1;
	double AFFINITY_ACTIVATION_FACTOR = 0.05;
	double NORMAL_MUTATION_RATE = 2.7e-9;
	double HYPERMUTATION_RATE = 3.2e-3;
	double MUTATION_TO_MEMORY_CELL_RATE = 0.01;
	double MUTATION_FROM_MEMORY_CELL_RATE_FACTOR = 0.4;
	double AFFINITY_THRESHOLD = 10000;
	string STARTING_DNA = "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCAGCTTTCGTAGTTATTGGATGACCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGGCCAACATATACCAAGATGGAAGTGAGCGACTTTATGGGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTTTCTCCAAATGAACAGCCTGAGAGTCGAGGACACGGCTGTGTACTACTGTGCGAGGCGAGCGAGCTACGGTGATTACGCGGTCCAAGTTAACCCCTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA";
	double TIME_INTERVAL = 0.1;
	double TIME_INTERVAL0 = 10;
	int CDR1_START = 30;
	int CDR1_END = 35;
	int CDR2_START = 49;
	int CDR2_END = 65;
	int CDR3_START = 97;
	int CDR3_END = 114;
	double NONE_CDR_MODIFIER = 0.5;
	int SAMPLE_SIZE = 1000;
	int NUM_TARGETS = 10;
	
	vector<string> HTARGET;
	vector<double> SWITCHING_TIME;
		
	Arguments(string filename){
		ifstream fin(filename);
		fin >> NORMAL_BIRTH_RATE;
		fin >> NONSENSE_DEATH_RATE;
		fin >> NORMAL_NONMEMORY_DEATH_RATE;
		fin >> NORMAL_DEATH_RATE;
		fin >> NORMAL_DEATH_FACTOR;
		fin >> INFACTED_BIRTH_RATE;
		fin >> INFACTED_DEATH_RATE;
		fin >> INFACTED_BASE_DEATH_FACTOR;
		fin >> AFFINITY_FACTOR;
		fin >> AFFINITY_ACTIVATION_FACTOR;
		fin >> NORMAL_MUTATION_RATE;
		fin >> HYPERMUTATION_RATE;
		fin >> MUTATION_TO_MEMORY_CELL_RATE;
		fin >> MUTATION_FROM_MEMORY_CELL_RATE_FACTOR;
		fin >> AFFINITY_THRESHOLD;
		fin >> STARTING_DNA;
		fin >> TIME_INTERVAL;
		fin >> TIME_INTERVAL0;
		fin >> CDR1_START;
		fin >> CDR1_END;
		fin >> CDR2_START;
		fin >> CDR2_END;
		fin >> CDR3_START;
		fin >> CDR3_END;
		fin >> NONE_CDR_MODIFIER;
		fin >> SAMPLE_SIZE;
		fin >> NUM_TARGETS;
		
		for (int i = 0; i < NUM_TARGETS; i++){
			string s;
			fin >> s;
			HTARGET.push_back(s);
			double f;
			fin >> f;
			SWITCHING_TIME.push_back(f);
		}
		
		cerr << NORMAL_BIRTH_RATE << endl;
		cerr << NONSENSE_DEATH_RATE << endl;
		cerr << NORMAL_NONMEMORY_DEATH_RATE << endl;
		cerr << NORMAL_DEATH_RATE << endl;
		cerr << NORMAL_DEATH_FACTOR << endl;
		cerr << INFACTED_BIRTH_RATE << endl;
		cerr << INFACTED_DEATH_RATE << endl;
		cerr << INFACTED_BASE_DEATH_FACTOR << endl;
		cerr << AFFINITY_FACTOR << endl;
		cerr << AFFINITY_ACTIVATION_FACTOR << endl;
		cerr << NORMAL_MUTATION_RATE << endl;
		cerr << HYPERMUTATION_RATE << endl;
		cerr << MUTATION_TO_MEMORY_CELL_RATE << endl;
		cerr << MUTATION_FROM_MEMORY_CELL_RATE_FACTOR << endl;
		cerr << AFFINITY_THRESHOLD << endl;
		cerr << STARTING_DNA << endl;
		cerr << TIME_INTERVAL << endl;
		cerr << TIME_INTERVAL0 << endl;
		cerr << CDR1_START << endl;
		cerr << CDR1_END << endl;
		cerr << CDR2_START << endl;
		cerr << CDR2_END << endl;
		cerr << CDR3_START << endl;
		cerr << CDR3_END << endl;
		cerr << NONE_CDR_MODIFIER << endl;
		cerr << SAMPLE_SIZE << endl;
		cerr << NUM_TARGETS << endl;
		for (int i = 0; i < NUM_TARGETS; i++){
			cerr << HTARGET[i] << endl;
			cerr << SWITCHING_TIME[i] << endl;
		}
	}
} const ARG("arguments.txt");

class FiveMerModel {
	unordered_map<string, array<double, 4> > p5, p4l, p4r, p3l, p3r;

public:	
	FiveMerModel(string fileName, double mutationRate){
		unordered_map<string, array<int, 4> > c5, c4l, c4r, c3l, c3r;
		int mutationCnt = 0, totalCnt = 0;
		ifstream fin(fileName);
		string seq, subseq;
		while (fin >> seq){
			int a, c, g, t;
			fin >> a >> c >> g >> t;
			c5[seq][0] += a; c5[seq][1] += c; c5[seq][2] += g; c5[seq][3] += t;
			subseq = seq.substr(0, 4); c4l[subseq][0] += a; c4l[subseq][1] += c; c4l[subseq][2] += g; c4l[subseq][3] += t;
			subseq = seq.substr(1, 4); c4r[subseq][0] += a; c4r[subseq][1] += c; c4r[subseq][2] += g; c4r[subseq][3] += t;
			subseq = seq.substr(0, 3); c3l[subseq][0] += a; c3l[subseq][1] += c; c3l[subseq][2] += g; c3l[subseq][3] += t;
			subseq = seq.substr(2, 3); c3r[subseq][0] += a; c3r[subseq][1] += c; c3r[subseq][2] += g; c3r[subseq][3] += t;
			totalCnt += a + c + g + t;
			char b = seq[2];
			if (b != 'A') mutationCnt += a;
			if (b != 'C') mutationCnt += c;
			if (b != 'G') mutationCnt += g;
			if (b != 'T') mutationCnt += t;
		}
		double r = mutationRate * totalCnt / mutationCnt;
		cnt2p(p5, c5, 2, r, mutationRate);
		cnt2p(p4l, c4l, 1, r, mutationRate);
		cnt2p(p4r, c4r, 2, r, mutationRate);
		cnt2p(p3l, c3l, 0, r, mutationRate);
		cnt2p(p3r, c3r, 2, r, mutationRate);
	}
	
	static void cnt2p(unordered_map<string, array<double, 4> > &pm, unordered_map<string, array<int, 4> > &cm, int pos, double r, double mutationRate){
		for (pair<const string, array<int, 4> > &p: cm){
			string seq = p.first;
			int a = p.second[0], c = p.second[1], g = p.second[2], t = p.second[3];
			double s = a + c + g + t;
			char b = seq[pos];
			if (a + c + g + t == 0) {
				if (b == 'A') pm[seq][0] = 1 - mutationRate; else pm[seq][0] = mutationRate / 3;
				if (b == 'C') pm[seq][1] = 1 - mutationRate; else pm[seq][1] = mutationRate / 3;
				if (b == 'G') pm[seq][2] = 1 - mutationRate; else pm[seq][2] = mutationRate / 3;
				if (b == 'T') pm[seq][3] = 1 - mutationRate; else pm[seq][3] = mutationRate / 3;
			}
			else {
				if (b == 'A') pm[seq][0] = 1 - (1 - a / s) * r; else pm[seq][0] = a / s * r;
				if (b == 'C') pm[seq][1] = 1 - (1 - c / s) * r; else pm[seq][1] = c / s * r;
				if (b == 'G') pm[seq][2] = 1 - (1 - g / s) * r; else pm[seq][2] = g / s * r;
				if (b == 'T') pm[seq][3] = 1 - (1 - t / s) * r; else pm[seq][3] = t / s * r;
			}
		}
	}
	
	vector<array<double, 4> > mutationProbabilities(string seq){
		int n = seq.length();
		vector<array<double, 4> > res;
		res.push_back(p3l[seq.substr(0, 3)]);
		res.push_back(p4l[seq.substr(0, 4)]);
		for (int i = 2; i < n - 2; i++){
			res.push_back(p5[seq.substr(i - 2, 5)]);
		}
		res.push_back(p4r[seq.substr(n - 4, 4)]);
		res.push_back(p3r[seq.substr(n - 3, 3)]);
		return res;
	}
	
} model("kmerFreq.txt", ARG.HYPERMUTATION_RATE);

class BCell: public Taxon{
public:
	static std::hash<string> str_hash;
	
	static int char2int(char c){
		if (c == 'A') return 0;
		if (c == 'C') return 1;
		if (c == 'G') return 2;
		return 3;
	}
	
	static char complement(char c){
		if (c == 'A') return 'T';
		if (c == 'C') return 'G';
		if (c == 'G') return 'C';
		return 'A';
	}
	
	static string complement(const string s){
		string res;
		for (char c: s){
			res += complement(c);
		}
		return res;
	}
	
	static int countC(const string s){
		int cnt = 0;
		for (char c: s){
			if (c == 'C') cnt++;
		}
		return cnt;
	}
	
	static int findC(const string s, double cnt){
		for (int i = 0; i < s.length(); i++){
			if (s[i] == 'C'){
				cnt -= 1;
				if (cnt <= 0) return i;
			}
		}
		return -1;
	}
	
	static char mutate(char c){
		char mc = c;
		while (mc == c){
			double r = uniform(generator);
			if (r < 0.25) mc = 'A';
			else if (r < 0.5) mc = 'C';
			else if (r < 0.75) mc = 'G';
			else mc = 'T';
		}
		return mc;
	}
	
	static char mutate(char c, array<double, 4> p){
		double s = 0;
		if (c != 'A') s += p[0];
		if (c != 'C') s += p[1];
		if (c != 'G') s += p[2];
		if (c != 'T') s += p[3];
		double r = uniform(generator) * s;
		if (c != 'A'){
			if (r < p[0]) return 'A';
			else r -= p[0];
		}
		if (c != 'C'){
			if (r < p[1]) return 'C';
			else r -= p[1];
		}
		if (c != 'G'){
			if (r < p[2]) return 'G';
			else r -= p[2];
		}
		return 'T';
	}
	
	static string translate(string dna){
		string aa;
		for (int i = 0; i < dna.length(); i += 3){
			aa += CODON.at(dna.substr(i, 3));
		}
		return aa;
	}
	
	static double score(const string a, const string b){
		int cdr = 0, fwr = 0;
		for (int i = 0; i < ARG.CDR1_START; i++){
			fwr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR1_START; i < ARG.CDR1_END; i++){
			cdr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR1_END; i < ARG.CDR2_START; i++){
			fwr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR2_START; i < ARG.CDR2_END; i++){
			cdr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR2_END; i < ARG.CDR3_START; i++){
			fwr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR3_START; i < ARG.CDR3_END; i++){
			cdr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		for (int i = ARG.CDR3_END; i < a.length(); i++){
			fwr += BLOSUM[AA2INT.at(a[i])][AA2INT.at(b[i])];
		}
		return cdr + fwr * ARG.NONE_CDR_MODIFIER;
	}
	
	static vector<bool> encodeDNA(const string dna){
		vector<bool> res;
		for (char c: dna){
			if (c == 'A') {
				res.push_back(false);
				res.push_back(false);
			}
			else if (c == 'C') {
				res.push_back(false);
				res.push_back(true);
			}
			else if (c == 'G') {
				res.push_back(true);
				res.push_back(false);
			}
			else {
				res.push_back(true);
				res.push_back(true);
			}
		}
		return res;
	}
	
	static string decodeDNA(const vector<bool> code){
		string res;
		for (int i = 0; i < code.size(); i += 2){
			if (code[i]){
				if (code[i + 1]) res += 'T';
				else res += 'G';
			}
			else {
				if (code[i + 1]) res += 'C';
				else res += 'A';
			}
		}
		return res;
	}
	
	const vector<bool> HF;
	const bool isMemoryCell;
	const bool nonSense;
	
	BCell(const string hf, bool memory = false): HF(encodeDNA(hf)), isMemoryCell(memory),
		nonSense(translate(hf).find('*') != string::npos){}
	BCell(const vector<bool> hf, bool memory = false): HF(hf), isMemoryCell(memory),
		nonSense(translate(decodeDNA(hf)).find('*') != string::npos){}
	
	virtual BCell* newBCell(const string hf, bool memory = false, bool copy = false) const{
		return new BCell(hf, memory);
	}
	
	virtual BCell* newBCell(const vector<bool> hf, bool memory = false, bool copy = false) const{
		return new BCell(hf, memory);
	}
	
	virtual bool operator == (const Taxon &other) const override{
		BCell &o = (BCell &) other;
		return HF == o.HF && isMemoryCell == o.isMemoryCell;
	}
	
	virtual size_t hash() const override{
		return str_hash(decodeDNA(HF)) * 33 + isMemoryCell;
	}
	
	double getBirthRate(int condition) const override{
		if (condition == 0) return ARG.NORMAL_BIRTH_RATE;
		else if (nonSense) return 0;
		else if (isMemoryCell) {
			const string HTARGET = ARG.HTARGET[condition - 1];
			return ARG.MUTATION_FROM_MEMORY_CELL_RATE_FACTOR * exp(-(score(HTARGET, HTARGET) - score(translate(decodeDNA(HF)), HTARGET)) * ARG.AFFINITY_ACTIVATION_FACTOR);
		}
		else return ARG.INFACTED_BIRTH_RATE;
	}
	
	double getNaturalDeathRate(int condition) const override{
		if (condition == 0) {
			if (isMemoryCell) return ARG.NORMAL_DEATH_RATE;
			else return ARG.NORMAL_NONMEMORY_DEATH_RATE;
		}
		else if (nonSense) return ARG.NONSENSE_DEATH_RATE;
		else if (isMemoryCell) return ARG.NORMAL_DEATH_RATE;
		else return ARG.INFACTED_DEATH_RATE;
	}
	
	double getLifeTimeMutationRate(int condition) const override{
		if (!isMemoryCell) return ARG.MUTATION_TO_MEMORY_CELL_RATE;
		return 0;
	}
	
	Taxon* mutateLifeTime(int conditon) const override{
		return newBCell(HF, !isMemoryCell);
	}
	
	double getBirthMutationProbability(int condition) const override{
		string h = decodeDNA(HF);
		if (condition == 0){
			return 1 - pow(1 - ARG.NORMAL_MUTATION_RATE, 2 * h.length());
		}
		else if (isMemoryCell) {
			return 1;
		}
		else {
			vector<array<double, 4> > p = model.mutationProbabilities(h);
			double t = 1;
			for (int i = 0; i < h.length(); i++){
				t *= p[i][char2int(h[i])];
			}
			return 1 - t * t;
		}
	}

	pair<Taxon*, Taxon*> mutateBirth(int condition) const override{
		string h = decodeDNA(HF), h1 = h, h2 = h;
		bool mutated = false;
		vector<double> p1(h1.length()), p2(h2.length());
		
		if (condition == 0){
			p2[h2.length() - 1] = 1 - ARG.NORMAL_MUTATION_RATE;
			for (int i = h2.length() - 2; i >= 0; i--){
				p2[i] = p2[i + 1] * (1 - ARG.NORMAL_MUTATION_RATE);
			}
			p1[h1.length() - 1] = p2[0] * (1 - ARG.NORMAL_MUTATION_RATE);
			for (int i = h1.length() - 2; i >= 0; i--){
				p1[i] = p1[i + 1] * (1 - ARG.NORMAL_MUTATION_RATE);
			}
			for (int i = 0; i < h.length(); i++){
				if (mutated){
					if (uniform(generator) < ARG.NORMAL_MUTATION_RATE){
						h1[i] = mutate(h[i]);
						mutated = true;
					}
				}
				else {
					if (uniform(generator) * (1 - p1[i]) < ARG.NORMAL_MUTATION_RATE){
						h1[i] = mutate(h[i]);
						mutated = true;
					}
				}
			}
			for (int i = 0; i < h.length(); i++){
				if (mutated){
					if (uniform(generator) < ARG.NORMAL_MUTATION_RATE){
						h2[i] = mutate(h[i]);
						mutated = true;
					}
				}
				else {
					if (uniform(generator) * (1 - p2[i]) < ARG.NORMAL_MUTATION_RATE){
						h2[i] = mutate(h[i]);
						mutated = true;
					}
				}
			}
		}
		else if (isMemoryCell){
			return {newBCell(h, true, true), newBCell(h, false, false)};
		}
		else {
			vector<array<double, 4> > p = model.mutationProbabilities(h);
			p2[h2.length() - 1] = p[h2.length() - 1][char2int(h[h2.length() - 1])];
			for (int i = h2.length() - 2; i >= 0; i--){
				p2[i] = p2[i + 1] * p[i][char2int(h[i])];
			}
			p1[h1.length() - 1] = p2[0] * p[h1.length() - 1][char2int(h[h1.length() - 1])];
			for (int i = h1.length() - 2; i >= 0; i--){
				p1[i] = p1[i + 1] * p[i][char2int(h[i])];
			}
			for (int i = 0; i < h.length(); i++){
				if (mutated){
					if (uniform(generator) < 1 - p[i][char2int(h[i])]){
						h1[i] = mutate(h[i], p[i]);
						mutated = true;
					}
				}
				else {
					if (uniform(generator) * (1 - p1[i]) < 1 - p[i][char2int(h[i])]){
						h1[i] = mutate(h[i], p[i]);
						mutated = true;
					}
				}
			}
			for (int i = 0; i < h.length(); i++){
				if (mutated){
					if (uniform(generator) < 1 - p[i][char2int(h[i])]){
						h2[i] = mutate(h[i], p[i]);
						mutated = true;
					}
				}
				else {
					if (uniform(generator) * (1 - p2[i]) < 1 - p[i][char2int(h[i])]){
						h2[i] = mutate(h[i], p[i]);
						mutated = true;
					}
				}
			}
		}
		return {newBCell(h1, isMemoryCell, !mutated), newBCell(h2, isMemoryCell, !mutated)};
	}
	
	double getOccupancyDeathRateFactor(int condition) const override{
		if (condition == 0) {
			if (isMemoryCell) return ARG.NORMAL_DEATH_FACTOR;
			else return 0;
		}
		if (nonSense) return ARG.NONSENSE_DEATH_RATE;
		if (isMemoryCell) return 0;
		else {
			const string HTARGET = ARG.HTARGET[condition - 1];
			return ARG.INFACTED_BASE_DEATH_FACTOR * exp((score(HTARGET, HTARGET) - score(translate(decodeDNA(HF)), HTARGET)) * ARG.AFFINITY_FACTOR);
		}
	}
	
	double getOccupancy(int condition) const override{
		if (condition == 0 && isMemoryCell) return 1;
		if (condition == 0 && !isMemoryCell) return 0;
		if (nonSense) return 0;
		if (isMemoryCell) return 0;
		else {
			const string HTARGET = ARG.HTARGET[condition - 1];
			return exp(-(score(HTARGET, HTARGET) - score(translate(decodeDNA(HF)), HTARGET)) * ARG.AFFINITY_FACTOR);	
		}
	}
};
std::hash<string> BCell::str_hash = std::hash<string>();

class BCellHomoplasy: public BCell{
	int id;
	
public:
	static int cnt;
	
	BCellHomoplasy(const string hf, bool memory = false, bool copy = false, int id = 0): BCell(hf, memory), id(copy ? id : cnt++){
		//cerr << getOccupancy(1) << " " << getOccupancyDeathRateFactor(1) << " " << getBirthMutationProbability(1) << endl;
		//system("pause");
	}
	BCellHomoplasy(const vector<bool> hf, bool memory = false, bool copy = false, int id = 0): BCell(hf, memory), id(copy ? id : cnt++){
		//cerr << getOccupancy(1) << " " << getOccupancyDeathRateFactor(1) << " " << getBirthMutationProbability(1) << endl;
		//system("pause");
	}
	
	BCell* newBCell(const string hf, bool memory = false, bool copy = false) const override{
		return new BCellHomoplasy(hf, memory, copy, id);
	}
	
	BCell* newBCell(const vector<bool> hf, bool memory = false, bool copy = false) const override{
		return new BCellHomoplasy(hf, memory, copy, id);
	}
	
	bool operator == (const Taxon &other) const override{
		BCellHomoplasy &o = (BCellHomoplasy &) other;
		return HF == o.HF && isMemoryCell == o.isMemoryCell && id == o.id;
	}
	
	size_t hash() const override{
		return str_hash(decodeDNA(HF)) * 33 + id;
	}
};
int BCellHomoplasy::cnt = 0;

int distance(const string a, const string b){
	int res = 0;
	for (int i = 0; i < a.length(); i++){
		if (a[i] != b[i]) res++;
	}
	return res;
}

string getInfo(Taxon &t){
	BCell &b = (BCell &) t;
	stringstream ss;
	ss << "DNA_distance_to_origin: " << distance(BCell::decodeDNA(b.HF), ARG.STARTING_DNA);
	string s;
	getline(ss, s);
	return s;
}

double cntMemory(Taxon &t){
	BCell &b = (BCell &) t;
	return (b.isMemoryCell) ? 1 : 0;
}

int main(int argc, char** argv){
	generator.seed(7);
	vector<double> conditionTime;
	cerr << BCell::translate(ARG.STARTING_DNA) << endl;
	Simulator s(vector<tuple<Taxon*, Taxon*, int, double> >({
		make_tuple<Taxon*, Taxon*, int, double>(new BCell/*Homoplasy*/(ARG.STARTING_DNA, true), nullptr, 1, 0.0)
	}), 0.0, 0, 1);
	s.debugInfo(getInfo);
	{
	ofstream fout("year_by_year.sif");
	double endofyear = 364;
	int year = 1999;
	for (int i = 1; i <= ARG.NUM_TARGETS; i++){
		s.switchCondition(i);
		while (s.getTotalOccupancy() < ARG.AFFINITY_THRESHOLD && s.getTime() < ARG.SWITCHING_TIME[i - 1]){
			s.simulate(min(ARG.TIME_INTERVAL, ARG.SWITCHING_TIME[i - 1] - s.getTime()));
			s.debugInfo();
			cerr << "Memory Cell Count: " << s.debugSum(cntMemory) << endl;
			if(s.getNumPop() == 0){
				cerr << "All dead.";
				return 0;
			}
		}
		s.switchCondition(0);
		while (s.getTime() < ARG.SWITCHING_TIME[i - 1]){
			s.simulate(min(ARG.TIME_INTERVAL0, ARG.SWITCHING_TIME[i - 1] - s.getTime()));
			s.debugInfo();
			cerr << "Memory Cell Count: " << s.debugSum(cntMemory) << endl;
		}
		conditionTime.push_back(s.getTime());
		if (s.getTime() > endofyear){
			vector<tuple<Population*, Population*, long long> > sample = s.sampleCondensed(20);
			for (auto &e: sample){
				fout << get<1>(e)->getTaxon() << " " << year << ((get<2>(e)) ? "_p_" : "_n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
			}
			endofyear += 365;
			year++;
		}
	}
	}
	{
		ofstream fout("1.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("2.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("3.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("4.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("5.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("6.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("7.sif");
		for (auto &e: s.sampleCondensed(ARG.SAMPLE_SIZE)){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	vector<tuple<Population*, Population*, long long> > sample = s.sampleCondensed(ARG.SAMPLE_SIZE);
	{
		ofstream fout("sample.sif");
		for (auto &e: sample){
			fout << get<1>(e)->getTaxon() << " " << ((get<2>(e)) ? "p_" : "n_") << BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF) << " " << get<0>(e)->getTaxon() << endl;
		}
	}
	{
		ofstream fout("sample_mutation.sif");
		for (auto &e: sample){
			string p = BCell::decodeDNA(((BCell*) get<1>(e)->getTaxon())->HF);
			string c = BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF);
			
			fout << p << " _";
			for (int i = 0; i < p.size(); i++){
				if (p[i] != c[i]) fout << i << p[i] << c[i];
			}
			fout << "_ " << c << endl;
		}
	}
	{
		ofstream fout("sample_mutation_position.sif");
		for (auto &e: sample){
			string p = BCell::decodeDNA(((BCell*) get<1>(e)->getTaxon())->HF);
			string c = BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF);
			
			fout << p << " _";
			for (int i = 0; i < p.size(); i++){
				if (p[i] != c[i]) fout << i;
			}
			fout << "_ " << c << endl;
		}
	}
	{
		ofstream fout("sample_condition.sif");
		for (auto &e: sample){
			string p = BCell::decodeDNA(((BCell*) get<1>(e)->getTaxon())->HF);
			string c = BCell::decodeDNA(((BCell*) get<0>(e)->getTaxon())->HF);
			
			fout << p << " ";
			for (int i = 0; i < conditionTime.size(); i++){
				if (get<0>(e)->getTimeOfEmergence() < conditionTime[i]){
					fout << i + 1;
					break;
				}
			}
			fout << " " << c << endl;
		}
	}
	return 0;
}
