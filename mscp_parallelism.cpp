#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <math.h>
#include <random>
#include <bits/random.h>
#include "include/BasicCDS.h"
#include <chrono>
#include <unistd.h>
#include <execution>
#include <algorithm>
#include <omp.h>
#include <numeric>
#include <iterator>
using namespace std;
using namespace std::chrono;
using namespace cds;

#define PRINT 0
#define CHECK 0 
#define TEST 1
#define REP 1

//NOTE: Check should be 1 only if you want step by step logs of the procedures.

//  Structure for each element of chi(universe).
typedef struct
{
	int value;
	int repetitions = 0;
	vector<int> inSet;
} item;
// Structure used for memory
typedef struct
{
	int bits;		// total number of bits needed.
	int nW;			// number of cells/words.
	int filled = 0; // summary of elements already covered
	ulong *mem;		// bits already covered.
} memory;
// Structure with all globals parameters program.
typedef struct
{
	vector<set<int>> F_sets; // Original sets of elements.
	vector<ulong *> bF_sets; // Binary representation of F nxM
	vector<ulong *> nF_sets; // Binary representation of F mxN
	vector<int> chi;		 // Universe elements
	int chiSize;
	unordered_map<int, int> chi_map; // Map each element of chi from 0 to n. (index in P_structure to the real value)
	unordered_map<int, int> back_map;
	set<int> greedy_solution;
	set<int> k_greedy_solution;
	vector<int> vidca_solution;
	int last_visited;
	int ns;
	memory M; // bits for subsets already covered.
	memory I; // bits for elements of Chi already covered.
	item *P;  // elements of chi, their cardinality and which subset included them.
	ifstream file;
	int numThreads;
	int threshold;

} ParProg;

// Generators
void genF(ParProg *par);
void readFsets(ParProg *par, ifstream &file); // new: To read from file.
void generateMatrixNM(ParProg *par);
void generateMatrixMN(ParProg *par);
void analyzeF(ParProg *par);
void initMemory(ParProg *par);
// Checks
void checkChi(ParProg *par);
void checkP(ParProg *par);
void checkMemory(ParProg *par);
// Main methods
void preSetCover(ParProg *par);
void setCover(ParProg *par);
void greedyAlg(ParProg *par);
void succintGreedyAlg(ParProg *par);
void kGreedyAlg(ParProg *par, int k);
// Subrutines.
bool compareRep(item a, item b);
bool compareRepMost(item a, item b);
void printVec(vector<set<int>> conjuntos);
void printSuccint(vector<ulong *> bSets, int chiSize);
int intersect(set<int> conjunto, vector<int> chi);
int countBits(ulong *e, int bits);
int checkBit(ulong *e, int pos);

set<int> kCandidates(map<int, set<int>>, vector<int>, int);
pair<int, int> bestC(ParProg *par, int pos);
pair<int, int> findCandidate(ParProg *par, int index);
// Test Methods
void testMemoryFilled(ParProg *par);

int main(int argc, char **argv)
{
	if (argc != 4)
	{
		cout << "FILE: Execution Error! call: ./runMSCP <fileName> <#threads> <#threshold>" << endl;
		exit(EXIT_FAILURE);
	}

	char FILENAME[60];
	strcpy(FILENAME, argv[1]);

	if (PRINT)
		cout << "Reading file..." << FILENAME << endl;
	char prefix[80];
	strcpy(prefix, FILENAME);
	strcat(prefix, "_results_");
	strcat(prefix, argv[2]);
	strcat(prefix, "_");
	strcat(prefix, argv[3]);

	ifstream file(FILENAME);
	int numThreads = atoi(argv[2]);
	omp_set_num_threads(numThreads);
	int threshold = atoi(argv[3]);
	FILE *fp = fopen(prefix, "w");

	auto sT = high_resolution_clock::now();
		ulong durTotalVidca = 0;
		ulong durTotalGreedy = 0;
		int cardinalidadVidca = 0;
		int cardinalidadGreedy = 0;

		ulong durAnalyze = 0;
		ulong durGenerate = 0;
		ulong durVidca = 0;
		ulong durHeu = 0;
		ulong durGreedy = 0;

		if (PRINT)
		{
			cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
			cout << "                            Exp #1" << endl;
			cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
		}
		fprintf(fp, "Exp #Threads %d - #Threshold %d\n", numThreads, threshold);
		ParProg *par = new ParProg();
		par->numThreads = numThreads;
		par->threshold = threshold;
		readFsets(par, file);
		auto ss = high_resolution_clock::now();
		generateMatrixNM(par);
		auto tt = high_resolution_clock::now();
		durGenerate = duration_cast<microseconds>(tt - ss).count();
		std::cout << "Tiempo de generacion: " << durGenerate << std::endl;
		file.close();
		for (int rep = 0; rep < REP; rep++)
		{
			ss = high_resolution_clock::now();
			analyzeF(par);
			tt = high_resolution_clock::now();
			auto dur_analyze = duration_cast<microseconds>(tt - ss);
			durAnalyze += dur_analyze.count();
			std::cout << "Tiempo del analisis: " << durAnalyze << std::endl;
			// checkP(par);
			// getchar();
			if (0)
			{
				cout << "Checking elements in P_structure:\n"
					 << endl;
				checkP(par);
				cout << "Checking chi, Universe of elements:\n"
					 << endl;
				checkChi(par);
				cout << "Checking F_sets read:\n"
					 << endl;
				printVec(par->F_sets);
				cout << "Printing F_sets as succint space:\n"
					 << endl;
				// printSuccint(par->bF_sets,par->chiSize);
			}
			// checkP(par);
			// cout<<"Printing F_sets as succint space:\n"<<endl;
			// printSuccint(par->bF_sets,par->chiSize);
			initMemory(par);
            if (CHECK) checkMemory(par);

			// Call of Greedy alg.
			ss = high_resolution_clock::now();
			succintGreedyAlg(par);
			kGreedyAlg(par, 2);
			tt = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(tt - ss);
			durGreedy += duration.count();
			std::cout << "Greedy: " << durGreedy << " microseconds." << std::endl;
			cardinalidadGreedy = par->greedy_solution.size();

			ss = high_resolution_clock::now();
			preSetCover(par);
			tt = high_resolution_clock::now();
			duration = duration_cast<microseconds>(tt - ss);
			durVidca = duration.count();
			std::cout << "(h1) Presetcover: " << durVidca << " microseconds." << std::endl;
			// checkMemory(par);
			// checkP(par);
			ss = high_resolution_clock::now();
			if (par->I.filled != par->I.bits)
			{
				setCover(par);
			}
			tt = high_resolution_clock::now();
			duration = duration_cast<microseconds>(tt - ss);
			durVidca += duration.count();
			std::cout << "(h2) Setcover: " << duration.count() << " microseconds." << std::endl;
			// checking vidca solution completeness.
			if (TEST)
				testMemoryFilled(par);
			cardinalidadVidca = par->vidca_solution.size();

			durHeu += durVidca;
			durTotalVidca += durVidca + durAnalyze;
			durTotalGreedy += durGreedy;
			if (PRINT)
				cout << "Finalized." << endl;
		durAnalyze = durAnalyze / REP;
		durTotalGreedy = durTotalGreedy / REP;
		durTotalVidca = durTotalVidca / REP;
		durHeu = durHeu / REP;
		fprintf(fp, "Total: %d %lu\n", cardinalidadVidca, durTotalVidca);
		fprintf(fp, "Generate: %lu\n", durGenerate);
		fprintf(fp, "Greedy: %d %lu\n", cardinalidadGreedy, durTotalGreedy);
		fprintf(fp, "Analyze: %lu\n", durAnalyze);
		fprintf(fp, "SetCover: %lu\n", durHeu);

		cout << "-----------------------------------------------------------------------------------------------" << endl;
	}
	fclose(fp);
	auto ee = high_resolution_clock::now();
	auto dd = duration_cast<microseconds>(ee - sT);
	cout << "Total: " << dd.count() << endl;
	return 0;
}
// Generate vector F of subsets S_i to S_n to simulate a universe of elements.
void readFsets(ParProg *par, ifstream &file)
{
	int sets;
	string line, item;
	getline(file >> std::ws, line);
	istringstream iss(line);
	getline(iss, item, ' ');
	par->chiSize = stoi(item);
	getline(iss, item, ' ');
	sets = stoi(item);
	par->ns = sets;
	if (1)
	{
		cout << "Sets: " << sets << endl;
		cout << "chi: " << par->chiSize << endl;
	}

	for (int j = 0; j < sets; j++)
	{
		getline(file >> std::ws, line);
		istringstream iss(line);
		getline(iss, item, ' ');
		getline(iss, item, ' ');
		set<int> act;
		while (getline(iss, item, ' '))
		{
			if (CHECK)
				cout << item << " ";
			act.insert(stoi(item));
		}
		if (CHECK)
			cout << endl;
		par->F_sets.push_back(act);
	}
	getline(file, line);
	if (PRINT)
		cout << "Subsets loaded." << endl;
}
void generateMatrixNM(ParProg *par)
{
	int sizeWord = 1 + (par->chiSize / (8 * sizeof(ulong)));
	int elemPos = 0;
    //for each set to be mapped.
	for (const set<int> &set_i : par->F_sets)
	{
		ulong *current_set = new ulong[sizeWord];
        //init as a 0 vector
        for (int word_j = 0; word_j < sizeWord; word_j++)
        {
            current_set[word_j] = 0;
        }
        
        //for each element of that set
		for (const int &elem_j : set_i)
		{
            //Find if it exist.
			auto it = (par->chi_map).find(elem_j);
			if (it == (par->chi_map).end())
			{
				par->chi_map[elem_j] = elemPos;
				elemPos++;
			}
            //If it exist or is it new, then map it.
			setBit64(current_set, par->chi_map[elem_j]);
		}
		par->bF_sets.push_back(current_set);
	}
}

// Preprocess  of F to obtain universe, cardinality of each element and location.
void analyzeF(ParProg *par)
{
	// define values to iterate
	const int nsize = par->chiSize;
	const vector<ulong *> bfSets = par->bF_sets;
	const int msize = par->ns;
	const int numWords = 1 + (nsize / (8 * sizeof(ulong)));

	// initialize universe and P structure.
	vector<int> v(nsize);
	std::iota(begin(v), end(v), 0);
	par->chi = v;
	par->P = new item[nsize];
	item *P = par->P;

#pragma omp parallel for
	for (int i_word = 0; i_word < numWords; i_word++)
	{ // for each parallel word
		for (int j_set = 0; j_set < msize; j_set++)
		{ // for each secuencial set
			ulong *act_set = bfSets[j_set];
			if (act_set[i_word] > 0)
			{ // if word != 0
				for (int bit = 0; bit < 64; bit++)
				{									   // for each bit of the word
					int off_bit = bit + (i_word * 64); // offset to access the n-bit, according to the i_word visited.
					if (checkBit(act_set, off_bit))
					{
						P[off_bit].inSet.push_back(j_set);
					}
				}
			}
		}
	}
	// after all the j_set are inserted, we can compute the number of repetitions
	// and assign it value equal to the original position in this structure (before reordering)

#pragma omp parallel for
	for (int i = 0; i < nsize; i++)
	{
		P[i].value = i;
		P[i].repetitions = (P[i].inSet).size();
	}
}
// Check P struct of items, to validate each step, ordering and repetitions if needed.
void checkP(ParProg *par)
{
	for (int cont = 0; cont < 50; cont++)
	{
		item ind = par->P[cont];
		cout << "#" << setw(2) << cont << " >> ";
		cout << "Val: " << setw(2) << ind.value << " | ";
		cout << "Rep.: " << setw(2) << ind.repetitions << " | ";
		cout << "In Sets:";
		// for (int i : ind.inSet)
		for (int a = 0; a < ind.repetitions && a < 10; a++)
		{
			cout << " " << ind.inSet[a];
		}
		cout << "." << endl;
	}
	cout << "\n---------------------------------------------------------" << endl;
	for (int cont = 4820; cont < par->chiSize; cont++)
	{
		item ind = par->P[cont];
		cout << "#" << setw(2) << cont << " >> ";
		cout << "Val: " << setw(2) << ind.value << " | ";
		cout << "Rep.: " << setw(2) << ind.repetitions << " | ";
		cout << "In Sets:";
		// for (int i : ind.inSet)
		for (int a = 0; a < ind.repetitions && a < 10; a++)
		{
			cout << " " << ind.inSet[a];
		}
		cout << "." << endl;
	}
	cout << "\n---------------------------------------------------------" << endl;
}
// Check cardinality of Chi and its elements.
void checkChi(ParProg *par)
{
	cout << "|X| = " << par->chiSize << endl;
	cout << "Elements of the universe X:" << endl;
	for (int val : par->chi)
	{
		cout << val << " ";
	}
	cout << endl
		 << "---------------------------" << endl;
}
// Initialize the memory of bits needed.
void initMemory(ParProg *par)
{
	par->I.bits = par->chiSize;
	par->I.nW = 1 + (par->I.bits / (8 * sizeof(ulong)));
	par->I.mem = new ulong[par->I.nW];
	for (int i = 0; i < par->I.nW; i++)
	{
		par->I.mem[i] = 0;
	}
	par->M.bits = par->F_sets.size();
	par->M.nW = 1 + (par->M.bits / (8 * sizeof(ulong)));
	par->M.mem = new ulong[par->M.nW];
	for (int i = 0; i < par->M.nW; i++)
	{
		par->M.mem[i] = 0;
	}
	if (CHECK)
	{
		cout << "Bits for Elements: " << par->I.bits << endl;
		cout << "Bits for Sets: " << par->M.bits << endl;
		checkMemory(par);
	}
}
// Print bits of ulong pointers of elements and subsets covered.
void checkMemory(ParProg *par)
{
	cout << "I-Memory: ";
	int cant = par->I.bits;
	int size = W64;
	for (int i = 0; i < par->I.nW; i++)
	{
		if (cant > size)
		{
			printFirstBits(par->I.mem[i], size);
			cant -= W64;
		}
		else
		{
			printFirstBits(par->I.mem[i], cant);
		}
	}
	auto start = high_resolution_clock::now();
	auto c = countBits(par->I.mem, par->I.bits);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "\ncount: " << c << endl;
	cout << "M-Memory: ";
	cant = par->M.bits;
	for (int i = 0; i < par->M.nW; i++)
	{
		if (cant > size)
		{
			printFirstBits(par->M.mem[i], size);
			cant -= 8 * sizeof(ulong);
		}
		else
		{
			printFirstBits(par->M.mem[i], cant);
		}
	}
	start = high_resolution_clock::now();
	c = countBits(par->M.mem, par->M.bits);
	stop = high_resolution_clock::now();
	duration = duration_cast<microseconds>(stop - start);
	cout << "\ncount: " << c << endl;
	cout << "---------------------------------------------------------" << endl;
}

// new Greedy Algorithm
void succintGreedyAlg(ParProg *par)
{
	if (PRINT)
	{
		cout << "---------------------------------------------------------" << endl;
		cout << "        Succint Greedy Algorithm in progress." << endl;
		cout << "---------------------------------------------------------" << endl;
		cout << "Succint Greedy Search over F_sets." << endl;
	}
	// vector<int> chi(par->chi);
	int nsize = par->chiSize;
	set<int> local_chi(par->chi.begin(), par->chi.end());
	// vector<ulong *> b_sets = par->bf_sets;
	map<int, ulong *> candidates;
	int posicionAbs = 0;
	for (ulong *each : par->bF_sets)
	{
		candidates[posicionAbs] = each;
		posicionAbs++;
	}

	while (local_chi.size() > 0)
	{
		if (CHECK)
			cout << "Length: " << candidates.size() << endl;
		int conjOptimo = 0;
		set<int> lowest_chi;
		int mayorCover = local_chi.size();

		for (pair<int, ulong *> toCheck : candidates)
		{
			set<int> temp_chi = local_chi;
			for (int pos = 0; pos < nsize; pos++)
			{
				if (checkBit(toCheck.second, pos))
				{
					temp_chi.erase(pos);
				}
			}
			if (temp_chi.size() < mayorCover)
			{
				mayorCover = temp_chi.size();
				lowest_chi = temp_chi;
				conjOptimo = toCheck.first;
			}
		}
		local_chi = lowest_chi;
		if (CHECK)
		{
			cout << "--Best candidate set : " << conjOptimo + 1;
			cout << "  >Covering: " << mayorCover << " new elements." << endl;
		}
		candidates.erase(conjOptimo);
		par->greedy_solution.insert(conjOptimo);
	}
	// cout<<"chi OG size: "<<(par->chi).size()<<endl;
	if (PRINT)
	{
		cout << "Greedy Alg. solution: " << endl;
		cout << ">>Cardinality: " << par->greedy_solution.size() << ", Sets: ";
		for (int j : par->greedy_solution)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}
// Greedy Algorithm
void greedyAlg(ParProg *par)
{
	if (PRINT)
	{
		cout << "---------------------------------------------------------" << endl;
		cout << "            Greedy Algorithm in progress." << endl;
		cout << "---------------------------------------------------------" << endl;
		cout << "Greedy Search over F_sets." << endl;
	}
	vector<int> chi(par->chi);
	map<int, set<int>> candidates;
	int posicionAbs = 0;
	for (set<int> each : par->F_sets)
	{
		candidates[posicionAbs] = each;
		posicionAbs++;
	}

	while (chi.size() > 0)
	{
		if (CHECK)
			cout << "Length: " << candidates.size() << endl;
		int conjOptimo = 0, mayorCover = 0, iCover = 0;
		for (pair<int, set<int>> par : candidates)
		{
			iCover = intersect(par.second, chi);
			if (iCover > mayorCover)
			{
				conjOptimo = par.first;
				mayorCover = iCover;
			}
		}
		if (CHECK)
		{
			cout << "--Best candidate set : " << conjOptimo + 1;
			cout << "  >Covering: " << mayorCover << " new elements." << endl;
		}
		for (int val : candidates[conjOptimo])
		{
			// if(CHECK) cout<<val << " ";
			std::vector<int>::iterator it = std::find(chi.begin(), chi.end(), val);
			if (it != chi.end())
				chi.erase(it);
		}
		candidates.erase(conjOptimo);
		par->greedy_solution.insert(conjOptimo);
	}
	if (0)
	{
		cout << "Greedy Alg. solution: " << endl;
		cout << ">>Cardinality: " << par->greedy_solution.size() << ", Sets: ";
		for (int j : par->greedy_solution)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}

// Greedy Algorithm
void kGreedyAlg(ParProg *par, int k)
{
	if (PRINT)
	{
		cout << "---------------------------------------------------------" << endl;
		cout << "           k-Greedy Algorithm in progress." << endl;
		cout << "---------------------------------------------------------" << endl;
		cout << "Greedy Search over F_sets." << endl;
	}
	vector<int> chi(par->chi);
	map<int, set<int>> candidates;
	int posicionAbs = 0;
	for (set<int> each : par->F_sets)
	{
		candidates[posicionAbs] = each;
		posicionAbs++;
	}

	while (chi.size() > 0)
	{
		if (CHECK)
			cout << "Length: " << candidates.size() << endl;
		set<int> bestCandidates;
		int mayorCover = 0;
		bestCandidates = kCandidates(candidates, chi, k);
		for (int set : bestCandidates)
		{
			for (int val : candidates[set])
			{
				if (CHECK)
					cout << val << " ";
				std::vector<int>::iterator it = std::find(chi.begin(), chi.end(), val);
				if (it != chi.end())
					chi.erase(it);
			}
		}
		for (int set : bestCandidates)
		{
			par->k_greedy_solution.insert(set);
			candidates.erase(set);
		}
	}
	if (CHECK)
	{
		cout << "K Greedy Alg. solution: " << endl;
		cout << ">>Cardinality: " << par->k_greedy_solution.size() << ", Sets: ";
		for (int j : par->k_greedy_solution)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}
set<int> kCandidates(map<int, set<int>> candidates, vector<int> chi, int k)
{
	set<int> bestSol; // subconjuntos a utilizar
	int actualLeftToSearch = candidates.size();
	if (CHECK)
		cout << actualLeftToSearch << " left to search" << endl;
	if (actualLeftToSearch < k)
	{
		for (int i = 0; i < actualLeftToSearch; i++)
		{
			bestSol.insert(i);
		}

		return bestSol;
	}

	int mayorCover = 0; // cardinalidad de la mejor solucion de k subconjuntos.
	// We iter over the possible k-length combinations of sets
	vector<bool> v(actualLeftToSearch);
	fill(v.begin(), v.begin() + k, true);
	std::vector<int>::iterator it;
	// Test every possible combination of index-length sets.
	do
	{
		set<int> candidatos;
		for (int i = 0; i < actualLeftToSearch; ++i)
		{
			if (v[i])
			{
				candidatos.insert(i);
			}
		}
		// Makes the union of elements of every set involved.
		set<int> unionActual;
		for (int setActual : candidatos)
		{
			for (int elemActual : candidates[setActual])
			{
				it = std::find(chi.begin(), chi.end(), elemActual);
				if (it != chi.end())
					unionActual.insert(elemActual);
			}
		}
		if (CHECK)
		{
			cout << "Candidatos:" << endl;
			cout << "< ";
			for (int c : candidatos)
			{
				cout << c << " ";
			}
			cout << ">, Con: " << unionActual.size() << " elementos." << endl;
		}
		if (unionActual.size() > mayorCover)
		{
			mayorCover = unionActual.size();
			bestSol = candidatos;
		}

	} while (prev_permutation(v.begin(), v.end()));
	return bestSol;
}

// Preprocess for setCover heuristic.
void preSetCover(ParProg *par)
{
	auto start = high_resolution_clock::now();
	std::sort(std::execution::par_unseq, par->P, (par->P) + par->I.bits, compareRep);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	// cout << "Sort takes: " << duration.count() << " microseconds." << endl;
	// cout << "Max Size:: " << par->P[par->chiSize - 1].inSet.size() << endl;
	if (CHECK)
	{
		checkP(par);
	}
	if (PRINT)
	{
		cout << "---------------------------------------------------------" << endl;
		cout << "            (H1) PreSetcover in progress." << endl;
		cout << "---------------------------------------------------------" << endl;
	}
	int pos = 0;
	if (CHECK)
		cout << "Adding unique chi elements.\n"
			 << endl;
	while (par->P[pos].repetitions == 1)
	{
		int unique_set = par->P[pos].inSet[0]; // select the index of the unique set that contain the element {P[pos]}
		if (!checkBit(par->I.mem, par->P[pos].value))
		{
			ulong *actual = (par->bF_sets[unique_set]); // preload the set from bf_sets
			if (CHECK)
				cout << "adding set: " << unique_set << endl;
			int actFilled = 0;
			for (int pos = 0; pos < par->I.nW; pos++)
			{
				par->I.mem[pos] |= (actual)[pos];
				actFilled += __builtin_popcountl(par->I.mem[pos]);
			}
			par->I.filled = actFilled;
			setBit64(par->M.mem, unique_set);		   // set the bit of the subset to 1.
			par->vidca_solution.push_back(unique_set); // add the identificator of the set visited.(0 to m)
		}
		pos++;
		par->last_visited = pos; // update the last visited element.
	}

	if (PRINT)
	{
		cout << "---------------------------------------------------------" << endl;
	}
}

void setCover(ParProg *par)
{
	if (PRINT)
	{
		cout << "            (H2) Setcover in progress." << endl;
		cout << "---------------------------------------------------------" << endl;
	}
	int index = par->last_visited; // We scan over P from index i to n  (P[i..n])
	int nsize = par->chiSize;
	int gradeSet = par->P[index].inSet.size();
	double init_h2, init_par, end_par, end_h2;
	while (par->I.filled < nsize && index < nsize)
	{
		int value = par->P[index].value;
		if (!checkBit(par->I.mem, value))
		{
			pair<int, int> setSelected = findCandidate(par, index);

			if (CHECK)
			{
				cout << "Se elige: " << setSelected.first << setSelected.second << endl;
				cout << "     Partiendo con: " << par->I.filled << " elementos." << endl;
			}
			for (int i_word = 0; i_word < par->I.nW; i_word++)
			{
				par->I.mem[i_word] |= (par->bF_sets[setSelected.first])[i_word];
			}
			par->I.filled += setSelected.second;
			setBit64(par->M.mem, setSelected.first);
			if (CHECK)
				cout << "     Termina con: " << par->I.filled << " elementos." << endl;

			par->vidca_solution.push_back(setSelected.first);
		}
		if (checkBit(par->I.mem, value))
			index++;
	}

}

// Horizontal parallelism implemented
pair<int, int> findCandidate(ParProg *par, int index)
{
	const int nsize = par->chiSize;
	const int msize = par->ns;
	const int nW = par->I.nW;
	const int numThreads = par->numThreads;
	const int threshold = par->threshold;

	int gradeSet = (par->P[index].inSet).size();
	int bestCandidate = 0;
	int bestSuma = 0;
	int indice = index;
	item *P = par->P;
	ulong *mem = par->I.mem;
	vector<ulong *> bF = par->bF_sets;
	// new segmente code
	int end = index + 1;
	while (end < nsize && (par->P[end].inSet).size() == gradeSet)
	{
		end++;
	}
	int rango = end - index;
	int maxThreads = numThreads;

	pair<int, int> partialSolutions[maxThreads];
	int sizeChunk = (rango + maxThreads - 1) / maxThreads;
    //parallelism init with maxThreads number of threads
#pragma omp parallel num_threads(maxThreads)
	{
		int tid = omp_get_thread_num();
		int start = tid * sizeChunk + index;
		int tt = start + sizeChunk;
		int bestPartialCandidate = 0;
		int bestPartialSuma = 0;
		for (int ss = start; ss < tt && ss < end; ss++)
		{
			if (!checkBit(mem, P[ss].value))
			{
				for (int posibleSet : par->P[ss].inSet)
				{
					ulong *toCheck = par->bF_sets[posibleSet];
					int sumAct = 0;
					for (int pos = 0; pos < par->I.nW; pos++)
					{
						ulong e = toCheck[pos] & (~(mem[pos]));
						sumAct += countBits(&e, W64);
					}

					if (sumAct > bestPartialSuma)
					{
						bestPartialSuma = sumAct;
						bestPartialCandidate = posibleSet;
					}
				}
			}
		}
		
		partialSolutions[tid] = pair<int, int>{bestPartialCandidate, bestPartialSuma};
	}
		for (pair<int, int> solution : partialSolutions)
		{
		
			if (solution.second > bestSuma)
			{
				bestSuma = solution.second;
				bestCandidate = solution.first;
			}
		}
	return pair<int, int>(bestCandidate, bestSuma);
}

int intersect(set<int> conjunto, vector<int> chi)
{
	int cont = 0;
	std::vector<int>::iterator it;
	for (int val : conjunto)
	{
		it = std::find(chi.begin(), chi.end(), val);
		if (it != chi.end())
		{
			cont++;
		}
	}
	return cont;
}
// Used to sort P structure by number of repetitions.
bool compareRep(item a, item b)
{
	return (a.repetitions < b.repetitions);
}
bool compareRepMost(item a, item b)
{
	return (a.repetitions > b.repetitions);
}
void printVec(vector<set<int>> conjuntos)
{
	int contador = 0;
	for (set<int> a : conjuntos)
	{
		cout << "(" << contador << ")"
			 << " < ";
		for (int b : a)
		{
			cout << b << " ";
		}
		contador++;
		cout << ">" << endl;
	}
	cout << "---------------------------" << endl;
}
void printSuccint(vector<ulong *> bSets, int chiSize)
{
	int numWord = 1 + chiSize / 64;
	for (ulong *e : bSets)
	{
		for (int i = 0; i < numWord; i++)
		{
			printBitsUlong(e[i]);
			cout << "-";
		}
		cout << endl;
	}
}

int countBits(ulong *e, int bits)
{
	int pos = 0;
	int ret = 0;
	int len = (int)bits / (W64 + 1);
	while (pos <= len)
	{
		ret += __builtin_popcountl(e[pos]);
		pos++;
	}
	return ret;
}
int checkBit(ulong *e, int pos)
{
	return ((((e[pos / W64]) >> (63 - (pos % W64))) & 1));
}
void testMemoryFilled(ParProg *par)
{
	cout << "---------------------------------------------------------" << endl;
	cout << "Universe cardinality: " << par->chiSize << endl;
	cout << "Testing Vidca solution..." << endl;
	;
	int dif = par->chiSize - par->I.filled;
	if (dif)
	{
		if (CHECK)
		{
			cout << "I-Memory: ";
			int cant = par->I.bits;
			int size = 8 * sizeof(ulong);
			for (int i = 0; i < par->I.nW; i++)
			{
				if (cant > size)
				{
					printFirstBits(par->I.mem[i], size);
					cant -= 8 * sizeof(ulong);
				}
				else
				{
					printFirstBits(par->I.mem[i], cant);
				}
			}
		}
		cout << "\n[Not completed] --" << dif << "-- elements left to cover." << endl;
	}
	else
	{
		cout << "All elements were covered in (H1) and (H2)." << endl;
	}
}
