#pragma once

#include <vector>
#include <algorithm>
#include "util.h"
#include "TSPMap.h"
#include "MainScene.h"
//#define MutationRate 0.75

enum SelectionType { Roulette, Tournament, ALT_Tournament, SelRandom};
enum CrossoverType {PMX, OBX, PBX, CroRandom};
enum MutationType {EM, DM, SM, IM, MutRandom};
enum FitnessScaleType {Sigma, Rank, Boltzmann, FitRandom};

struct Genome{
	vector<int> mv_gene; //use vector to store the gene
	double md_fitness; //the fitness of this genome
	Genome() {}
	Genome(int n) {
		mv_gene.reserve(n);
		for (int i = 0; i < n; ++i) mv_gene.push_back(i);
		for (int i = 0; i < n; ++i) {
			int x = randInt(0, n - 1);
			int y = randInt(0, n - 1);
			if (x == y) continue;
			swap(mv_gene[x], mv_gene[y]);
		}
	}

	int at(int k) const { return mv_gene[k]; }
	int& operator[](int n) { return mv_gene[n]; }
	int size() { return mv_gene.size(); }
	bool operator<(const Genome& g) const{ return md_fitness < g.md_fitness; }
	bool operator==(const Genome& g)const { return mv_gene == g.mv_gene; }
};

class TSPGA {
public:
	const double Inf = (1 << 30);
	bool mb_running;
	bool mb_started;

	double md_sigma;
	double md_bestFitness;
	double md_worstFitness;
	double md_aveFitness;
	double md_longestRoute;
	double md_shortestRoute;
	double md_totalFitness;
	//used for fitnessScaleBoltzmann()
	double md_boltzmannTemp;
	double md_boltzmannMinTemp = 1;
	double md_boltzmannDT = 0.05;

	int mi_populationSize;
	int mi_timeOfCity;
	int mi_numOfEliteToChoose;
	int mi_bestOne; 
	int mi_generation;

	vector<Genome> mv_popu;
	Genome mg_bestGene;
	int mi_cntBestTime;
	int mi_maxBestTime;
	TSPMap *mp_map;
	MainScene* mp_scene;

	
	SelectionType ST[4] = { Roulette, Tournament, ALT_Tournament, SelRandom};
	CrossoverType CT[4] = { PMX, OBX, PBX , CroRandom};
	MutationType MT[5] = {IM, DM, SM, EM, MutRandom};
	FitnessScaleType FST[4] = { Sigma, Rank, Boltzmann ,FitRandom};

	MutationType mutationType = IM;
	SelectionType selectionType = Tournament;
	CrossoverType crossoverType = PMX;
	FitnessScaleType scaleType = Sigma;
	bool mb_mutateRandom = false;
	bool mb_selectRandom = false;
	bool mb_crossoverRandom = false;
	bool mb_fitnessRandom = false;

	TSPGA(vector<City*> &cities, MainScene* mainScene) {
		this->mp_map = new TSPMap();
		this->mp_scene = mainScene;
		mp_map->setCitys(cities);
		mi_timeOfCity = 20;
		mi_populationSize = mp_map->numberOfCities() * mi_timeOfCity;
		mi_generation = 0;
		mi_cntBestTime = 0;
		mb_running = false;
		mb_started = false;
	}

	void setCities(vector<City*> &cities) {
		mp_map->setCitys(cities);
		mi_populationSize = mp_map->numberOfCities() * mi_timeOfCity;
	}

	vector<Vec2> getBestRoute() {
		vector<Vec2> bestRoute;
		for (int i = 0, n = mg_bestGene.size(); i < n; ++i) {
			//city number start from 1 in gene
			bestRoute.push_back(mp_map->PosOfCity(mg_bestGene.at(i)));
		}
		return bestRoute;
	}


	void epoch() {
		reset();
		calculatePopulationFitness();

		fitnessScaleSwitch();

		calculateStatistics();

		//create a temporary vector for the new population
		vector<Genome> vecNewPopu=mv_popu;

		int doublePopu = mi_populationSize * 2;
		while (vecNewPopu.size() < doublePopu) {
			Genome dad = selectionSwitch(mv_popu);
			Genome mum = selectionSwitch(mv_popu);
			Genome baby1, baby2;

			crossoverSwitch(dad, mum, baby1, baby2);

			mutationSwitch(baby1, baby2);

			vecNewPopu.push_back(baby1);
			vecNewPopu.push_back(baby2);
		}

		//select from the existing population and the new generation
		mv_popu.clear();
		grabBestN(mi_numOfEliteToChoose, vecNewPopu, mv_popu);
	//	grabBestN(0, vecNewPopu, mv_popu);
		while (mv_popu.size() < mi_populationSize) {
			mv_popu.push_back(selectionSwitch(vecNewPopu));
		}
		++mi_generation;
	}

	void reset() {
		//md_sigma = 1;
		md_bestFitness = 0;
		//md_worstFitness = Inf;
		//md_longestRoute = 0;
		//md_shortestRoute = Inf;
		md_aveFitness = 0;
		md_totalFitness = 0;
		mi_numOfEliteToChoose = 0.005 * mi_populationSize;
		if (mi_numOfEliteToChoose < 1) mi_numOfEliteToChoose = 1;
		mi_bestOne = 0;

		mi_maxBestTime = mp_map->numberOfCities() * 0.5;
	}

	void calculatePopulationFitness() {
		double longest = 0;
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			mv_popu[i].md_fitness = mp_map->lengthOfRouteSquare(mv_popu[i].mv_gene);
			if (mv_popu[i].md_fitness > longest)  longest = mv_popu[i].md_fitness;
			md_totalFitness += mv_popu[i].md_fitness;
		}
		md_aveFitness = md_totalFitness / mv_popu.size();

		//calculate fitness
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			mv_popu[i].md_fitness = longest - mv_popu[i].md_fitness;
			if (mv_popu[i].md_fitness > md_bestFitness) md_bestFitness = mv_popu[i].md_fitness;
		}

		//if more then 80% of population have the same value(length of route) 
		//with the best one,then stop
		//if ((md_totalFitness / (longest*mv_popu.size())) >= 0.80) mb_running = false;

	}

	//-----------------------FitnessScale-------------------------//

	void fitnessScaleSwitch() {
		switch (scaleType) {
		case Sigma:
			fitnessScaleSigma();
			break;
		case Rank:
			fitnessScaleRank();
			break;
		case Boltzmann:
			fitnessScaleBoltzmann();
			break;
		default:
			fitnessScaleRank();
		}
		scaleType = FST[randInt(0, 2)];
	}

	//recalculate fitness according to the standard deviation
	void fitnessScaleSigma() {
		double sum = 0;
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			sum += (mv_popu[i].md_fitness - md_aveFitness) *
				   (mv_popu[i].md_fitness - md_aveFitness);
		}
		double variance = sum / mv_popu.size();
		//standard deviation
		md_sigma = sqrt(variance);

		//calculate the scaled fitness for each member
		double sigme2 = md_sigma * 2;
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			mv_popu[i].md_fitness = (mv_popu[i].md_fitness - md_aveFitness) / sigme2;
		}
	}

	
	//functor used to sort population
	struct CMPGenome {
		vector<Genome>* popu;
		CMPGenome(vector<Genome>* population):popu(population) {}
		bool operator()(const int g1, const int g2) { return popu->at(g1) < popu->at(g2);}
	};

	//recalculate fitness according the rank of the genome
	void fitnessScaleRank() {
		vector<int> v;
		v.reserve(mv_popu.size());
		for (int i = 0, n = mv_popu.size(); i < n; ++i) v.push_back(i);
		std::sort(v.begin(), v.end(), CMPGenome(&mv_popu));

		for (int i = 0, n = v.size(); i < n; ++i) mv_popu[v[i]].md_fitness = i;
	}

	// The value md_boltzmannTemp is the boltzmann temperature which is reduced
	//  each generation by a small amount. As Temp decreases the difference 
	//  spread between the high and low fitnesses increases.
	void fitnessScaleBoltzmann() {
		if (md_boltzmannTemp > md_boltzmannMinTemp) {
			md_boltzmannTemp -= md_boltzmannDT;
		}
		//first calculate the average fitness/Temp
		double divider = md_aveFitness / md_boltzmannTemp;
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			mv_popu[i].md_fitness = mv_popu[i].md_fitness / md_boltzmannTemp / divider;
		}
	}

	//calculate the statistics information after the fitness has been scaled
	void calculateStatistics() {
		md_totalFitness = 0;
		for (int i = 0, n = mv_popu.size(); i < n; ++i) {
			md_totalFitness += mv_popu[i].md_fitness;
		}
		md_aveFitness = md_totalFitness / mv_popu.size();
	//	if (fabs(md_bestFitness - md_aveFitness) < md_totalFitness * 0.2) mb_running = false;
	//	if (md_sigma / mv_popu.size() < 3000.0) mb_running = false;
	}


	//create the initial population
	void createStartingPopulation() {
		mv_popu.clear();
		int numofCities = mp_map->numberOfCities();
		for (int i = 0; i < mi_populationSize; ++i) 
			mv_popu.push_back(Genome(numofCities));
		mi_generation = 0;
		md_boltzmannTemp = mp_map->numberOfCities() * 2;
	}


	//copy the best gene n times to the new population
	void grabBestN(int n, vector<Genome> &src, vector<Genome> &dest) {
		int best = 0;
		for (int i = 1,m= src.size(); i < m; ++i) {
			if (src[best].md_fitness < src[i].md_fitness) best = i;
		}

		if (mg_bestGene == src[best]) {
			if (++mi_cntBestTime >= mi_maxBestTime) {
				mb_running = false;
				mi_cntBestTime = 0;
			}
		}
		else {
			mi_cntBestTime = 0;
			mg_bestGene = src[best];

		}
		
		for (int i = 0; i < n; ++i) {
			dest.push_back(src[best]);
		}
	}

	//---------------------selection---------------------------//
	// Roulette, Tournament, ALT_Tournament
	Genome selectionSwitch(vector<Genome>& src) {
		const int TournamentTry = 10;
		switch (selectionType) {
		case Roulette: 
			return selectRoulette(src);
		case Tournament: 
			return selectTournament(TournamentTry, src);
		case ALT_Tournament: 
			return selectALTTournament(src);
		default: 
			return selectALTTournament(src);
		}
		selectionType = ST[randInt(0,2)];
	}

	Genome selectRoulette(vector<Genome>& src) {
		double fSlice = randFloat() * md_totalFitness;
		double sum = 0;
		for (int i = 0, n = src.size(); i < n; ++i) {
			sum += src[i].md_fitness;
			if (sum >= fSlice) return src[i];
		}
	}

	//select the best one in this n trys
	Genome selectTournament(int numToTry, vector<Genome>& src) {
		int best = randInt(0, src.size()-1);
		for (int i = 0; i < numToTry; ++i) {
			int thisTry = randInt(0, src.size() - 1);
			if (src[best].md_fitness < src[thisTry].md_fitness) best = thisTry;
		}
		return src[best];
	}

	//select the better one of the randomly selected two Genome in 80% percentage,
	//and the weaker one in 20% percentage.
	Genome selectALTTournament(vector<Genome>& src) {
		int g1 = randInt(0, src.size() - 1);
		int g2 = randInt(0, src.size() - 1);
		while (g1 == g2) g2 = randInt(0, src.size() - 1);
		double r = randFloat();
		if (r < 0.80) {
			return src[g1] < src[g2] ? src[g2] : src[g1];
		}
		else {
			return src[g1] < src[g2] ? src[g1] : src[g2];
		}
	}


	//---------------------Crossover---------------------------//

	void crossoverSwitch(Genome& dad, Genome& mum, Genome &baby1, Genome &baby2) {
		switch (crossoverType) {
		case PMX:
			crossoverPMX(dad, mum, baby1, baby2);
			break;
		case PBX:
			crossoverPBX(dad, mum, baby1, baby2);
			break;
		case OBX:
			crossoverOBX(dad, mum, baby1, baby2);
			break;
		default:
			crossoverPMX(dad, mum, baby1, baby2);
			break;
		}
		if(mb_crossoverRandom) crossoverType = CT[randInt(0, 2)];
	}

	/*permutation based crossover
	*randomly choose a range, for each pair(one from mama and one from dada) of gene 
	* in the range, swap their position in genome for both baby
	*/
	void crossoverPMX(Genome& dad, Genome& mum, Genome &baby1, Genome &baby2) {
		baby1 = dad;
		baby2 = mum;
		if (dad == mum) return;
		int beg = randInt(0, dad.size()-1);
		int end = randInt(beg, dad.size() - 1);
		for (int pos = beg; pos <= end; ++pos) {
			int gene1 = dad[pos];
			int gene2 = mum[pos];
			if (gene1 == gene2) continue; 
			int pos1, pos2;
			//find and swap the position of gene1 and gene2
			for (int i = 0, n = baby1.size(); i < n; ++i) {
				if (baby1[i] == gene1) pos1 = i;
				else if (baby1[i] == gene2) pos2 = i;
			}
			swap(baby1[pos1], baby1[pos2]);
			//and for baby2
			for (int i = 0, n = baby2.size(); i < n; ++i) {
				if (baby2[i] == gene1) pos1 = i;
				else if (baby2[i] == gene2) pos2 = i;
			}
			swap(baby2[pos1], baby2[pos2]);
		}
	}

	/*order based crossover
	* randomly choose some position in the genome, apply the order of the selected
	* gene in dad to mum, and vice versus.
	*/
	void crossoverOBX(Genome& dad, Genome& mum, Genome &baby1, Genome &baby2) {
		baby1 = dad;
		baby2 = mum;
		if (mum == dad) return;
		vector<int> cities;
		vector<int> positions;
		int n = dad.size();

		int pos = randInt(0, n - 2);
		//randomly choose some cities from dad and record their positions
		while (pos < n) {
			positions.push_back(pos);
			cities.push_back(dad[pos]);
			//notice the range, the return value must >=1,otherwise the same cities 
			//maybe choosed twice, which will result in problem
			pos += randInt(1, n - pos); 
		}

		//now apply the order of the selected cities to mum to generate baby2
		pos = 0;
		for (int cit = 0; cit < n; ++cit) {
			for (int i = 0; i < cities.size(); ++i) {
				if (baby2[cit] == cities[i]) {
					baby2[cit] = cities[pos++];
					break;
				}
			}
		}

		//now grab the cooresponding cities from mum 
		cities.clear();
		for (int i = 0; i < positions.size(); ++i) {
			cities.push_back(mum[positions[i]]);
		}

		//and impose their order in dad to generate baby1
		pos = 0;
		for (int cit = 0; cit < n; ++cit) {
			for (int i = 0; i < cities.size(); ++i) {
				if (baby1[cit] == cities[i]) {
					baby1[cit] = cities[pos++];
					break;
				}
			}
		}
	}

	
	/*position based crossover
	* randomly select some positions, the cities in those positions will  be
	* copied to offsprings to the same positions£¬the rest positions will be fill
	* with cities from the other parent.
	*/
	void crossoverPBX(Genome& dad, Genome& mum, Genome &baby1, Genome &baby2) {
		if (dad == mum) {
			baby1 = dad;
			baby2 = mum;
			return;
		}
		baby1.mv_gene.reserve(mum.size());
		baby2.mv_gene.reserve(mum.size());
		baby1.mv_gene.assign(mum.size(), -1);
		baby2.mv_gene.assign(mum.size(), -1);

		//randomly select some positions, the cities in those position will  be
		//copied to offsprings to the same position
		vector<int> positions;
		for (int n = mum.size(), pos = randInt(0, n - 2); pos < n; pos += randInt(0, n - pos)) {
			positions.push_back(pos);
		}

		//copy the selected cities to offsprings to the same position
		for (int i = 0, n = positions.size(); i < n; ++i) {
			//baby1 receives selected cities from dad
			baby1[positions[i]] = dad[positions[i]];
			//baby2 receives selected cities from mum
			baby2[positions[i]] = mum[positions[i]];
		}

		//fill the rest positions with cities from the other parent.
		int c1=0, c2=0;
		for (int n = mum.size(), pos = 0; pos < n; ++pos) {
			//baby1's free positions will be filled with cities from mum
			while (c1 < n && baby1[c1] > -1) ++c1; //find a free position
			if (c1<n && find(baby1.mv_gene.begin(), baby1.mv_gene.end(), mum[pos]) 
				== baby1.mv_gene.end()) {
				baby1[c1] = mum[pos];
			}

			//baby2's free positions will be filled with cities from dad
			while (c2 < n && baby2[c2] > -1) ++c2;
			if (c2<n && find(baby2.mv_gene.begin(), baby2.mv_gene.end(), dad[pos]) 
				== baby2.mv_gene.end()) {
				baby2[c2] = dad[pos];
			}
		}
	}


	//---------------------mutation---------------------------//
	void mutationSwitch(Genome& baby1, Genome &baby2) {
		switch (mutationType) {
		case IM:
			mutateIM(baby1);
			mutateIM(baby2);
			break;
		case EM:
			mutateEM(baby1);
			mutateEM(baby2);
			break;
		case  DM:
			mutateDM(baby1);
			mutateDM(baby2);
			break;
		case SM:
			mutateSM(baby1);
			mutateSM(baby2);
			break;
		default:
			mutateIM(baby1);
			mutateIM(baby2);
		}
		//randomly change the mutation type
		mutationType = MT[randInt(0, 3)];
	}


	//randomly choose two genes and swap their position.
	void mutateEM(Genome& genome) {
		int gene1 = randInt(0, genome.size() - 1);
		int gene2 = randInt(0, genome.size() - 1);
		while (gene1 == gene2) {
			gene2 = randInt(0, genome.size() - 1);
		}
		swap(genome[gene1], genome[gene2]);
	}

	//randomly choose a gene and move it to another random place.
	void mutateIM(Genome& genome) {
		//randomly choose a gene
		vector<int>::iterator pos = genome.mv_gene.begin() + randInt(0, genome.size()-1);
		int city = *pos; //get it
		//remove it from the genome
		genome.mv_gene.erase(pos); 
		//randomly choose a place
		pos = genome.mv_gene.begin() + randInt(0, genome.size() - 1); 
		//insert it to this place
		genome.mv_gene.insert(pos, city); 
	}

	//randomly choose a gene and disturbs its order randomly
	void mutateSM(Genome& genome) {
		int beg = randInt(0, genome.size()-1);
		int end = randInt(0, genome.size()-1);
		while (beg == end) end = randInt(0, genome.size()-1);
		if (beg > end) swap(beg, end);
		//if there are noly two genes in the segment, just swap this two
		if (end-beg == 1) {
			swap(genome[beg], genome[end]);
			return;
		}
		//swap cnt times, cnt is half length of the range
		for (int i = 0, cnt = (end - beg + 1) / 2; i<cnt; ++i) {
			int x = randInt(beg, end);
			int y = randInt(beg, end);
			if (x == y) continue;
			swap(genome[x], genome[y]);
		}
	}

	//randomly choose a segment of gene and move it to another random place
	void mutateDM(Genome& genome) {
		//make sure that S < E < genome.size()
		int S = randInt(0, genome.size() - 2);
		int E = randInt(S+1, genome.size() - 1);
		vector<int>::iterator beg = genome.mv_gene.begin() + S;
		vector<int>::iterator end = genome.mv_gene.begin() + E;
		vector<int> seg(beg, end);
		genome.mv_gene.erase(beg, end);
		vector<int>::iterator insertPos = genome.mv_gene.begin() + randInt(0, genome.size()-1);
		genome.mv_gene.insert(insertPos, seg.begin(), seg.end());
	}
	

	//----------------setter getter-----------------------------//
	int getCrossoverType(){ return crossoverType; }
	void setCrossoverType(int type) { crossoverType = (CrossoverType)type; }
};