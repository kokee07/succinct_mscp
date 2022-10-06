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

using namespace std;
using namespace std::chrono; 
using namespace cds;

#define PRINT 0
#define CHECK 0
#define TEST 1

// Structure for each element of chi(universe).
typedef struct{
	int value;
	int repetitions=0;
	vector<int> inSet;
} item;
//Structure used for memory
typedef struct{
	int bits; //total number of bits needed.
	int nW; //number of cells/words.
	int filled=0;//summary of elements already covered
	ulong *mem; //bits already covered.
}memory;
// Structure with all globals parameters program.
typedef struct {
	vector<set<int>> F_sets; //Original sets of elements.
	vector<ulong*> bF_sets; //Binary representation of F
    set<int> chi; //Universe elements
	int chiSize;
	map<int,int> chi_map;//Map each element of chi from 0 to n. (index in P_structure to the real value)
	set<int> exha_solution;
    set<int> greedy_solution;
    set<int> vidca_solution;
	int last_visited;
    int ns;
	memory M;//bits for subsets already covered.
	memory I;//bits for elements of Chi already covered.
	item * P; //elements of chi, their cardinality and which subset included them.
    ifstream file;
} ParProg;

//Generators 
void genF(ParProg *par);
void read_Fset(ParProg *par, ifstream& file); //new: To read from file.
void analyzeF(ParProg *par);
void initMemory(ParProg *par);
//Checks
void checkChi(ParProg *par);
void checkP(ParProg *par);
void checkMemory(ParProg *par);
//Main methods
void preSetCover(ParProg *par);
void setCover(ParProg *par);
void h2(ParProg *par);
void greedyAlg(ParProg *par);
//Subrutines.
bool compareRep(item a,item b);
bool compareRepMost(item a,item b);
int dist(ParProg *par,int index);
void printVec(vector<set<int>> conjuntos);
void printSuccint(vector<ulong*> bSets, int chiSize);
int intersect(set<int> conjunto, set<int> chi);
int countBits(ulong* e, int bits);
int checkBit(ulong* e, int pos);
int bestCandidate(ParProg *par, int pos);
pair<int,int> bestC(ParProg *par, int pos);
pair<int,int> findCandidate(ParProg* par, int index);
//Test Methods
void testMemoryFilled(ParProg *par);


int main(int argc, char** argv){
	if(argc !=2){
		cout << "FILE: Execution Error! call: ./vidca <fileName>" << endl;
		exit(EXIT_FAILURE);
	}
	
	char FILENAME[40];
	strcpy(FILENAME,argv[1]);

	if(PRINT) cout << "Reading file..." <<FILENAME<< endl;
	char prefix[80];
	strcpy(prefix,FILENAME);
	strcat(prefix,"_results_HEU");
	FILE *fp = fopen(prefix, "w" );
	
    int experiment=1;
    ifstream file(FILENAME);
	auto ss = high_resolution_clock::now(); 
    while (!file.eof()) {
		ulong durTotalVidca=0;
		ulong durTotalGreedy=0;
		float cardinalidadVidca=0;
		float cardinalidadGreedy=0;
		string line,item;
        getline(file,line);
        if(file.eof()) break;
        int REP=stoi(line);
        //REP = 1;
        if(PRINT){
            cout<<"Initializing repetitions..."<<endl;
            cout<< "REP: "<<REP<<endl;
        }
        for (int k = 0; k < REP; k++)
        {	
            ulong durVidca=0;
            ulong durGreedy=0;
            if(PRINT){
				cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
				cout<< "                            Exp #"<< experiment<<", Iter #"<<k+1<<"."<<endl;
				cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
			}
			fprintf(fp,"Exp #%d, Iter #%d\n",experiment,k+1);
            ParProg *par = new ParProg();
            read_Fset(par,file);
            auto ss = high_resolution_clock::now(); 
            analyzeF(par);
            auto tt = high_resolution_clock::now(); 
            auto dur_analyze = duration_cast<microseconds>(tt-ss).count(); 
            if(CHECK){
                cout<<"Checking elements in P_structure:\n"<<endl;
                checkP(par);
                cout<<"Checking chi, Universe of elements:\n"<<endl;
                checkChi(par);
                cout<<"Checking F_sets read:\n"<<endl;
                printVec(par->F_sets);
				cout<<"Printing F_sets as succint space:\n"<<endl;
                printSuccint(par->bF_sets,par->chiSize);
            }
			//cout<<"Printing F_sets as succint space:\n"<<endl;
            //printSuccint(par->bF_sets,par->chiSize);
            initMemory(par);
            if(CHECK) checkMemory(par);
			
			//Call of Greedy alg.
			auto start = high_resolution_clock::now(); 
			greedyAlg(par);
			auto stop = high_resolution_clock::now(); 
			auto duration = duration_cast<microseconds>(stop - start); 
			durGreedy=duration.count();
			cardinalidadGreedy+=par->greedy_solution.size();
			start = high_resolution_clock::now(); 
			preSetCover(par);
			stop = high_resolution_clock::now(); 
			duration = duration_cast<microseconds>(stop - start);
			durVidca=duration.count();
			start = high_resolution_clock::now(); 
			if(par->I.filled!=par->I.bits){
				setCover(par);
			}
			stop = high_resolution_clock::now(); 
			duration = duration_cast<microseconds>(stop - start);
			durVidca+=duration.count();
			//checking vidca solution completeness.
			if(TEST) testMemoryFilled(par);
		
			cout<<"G:"<<par->greedy_solution.size()<<", V:"<<par->vidca_solution.size()<<endl;
            cardinalidadVidca+=par->vidca_solution.size();
			fprintf(fp,"%lu %lu\n",par->greedy_solution.size(),par->vidca_solution.size());
            fprintf(fp,"%lu %lu\n",durGreedy,durVidca+dur_analyze);
            durTotalVidca+=durVidca+dur_analyze;
            durTotalGreedy+=durGreedy;
            if(PRINT) cout<<"Finalized #"<<k<<endl;
			//cout<<"---------------------------------------------------------------------------"<<endl;
        }
        experiment++;
        cout << "################## " << endl;
        cout << "Time taken by "<<REP<<" reps.: "
			<<durTotalVidca+durTotalGreedy<< " microseconds" << endl; 
		cout << "Average time per cycle of Greedy: "<<durTotalGreedy/REP<< " microseconds"<<endl;
		cout << "Average cardinality per cycle Greedy: "<<cardinalidadGreedy/REP<< " subsets."<<endl;
		cout << "Average time per cycle of Vidca: "<<durTotalVidca/REP<< " microseconds"<<endl;
		cout << "Average cardinality per cycle Vidca: "<<cardinalidadVidca/REP<< " subsets."<<endl;
		cout<<"-----------------------------------------------------------------------------------------------"<<endl;

        fprintf(fp,"%f %f\n",cardinalidadGreedy/REP,cardinalidadVidca/REP);
    }
	file.close();
    fclose (fp);
	auto ee = high_resolution_clock::now(); 
	auto dd = duration_cast<microseconds>(ee - ss);
	cout<<"Total: "<<dd.count()<<endl;
	return 0;
}
//Generate vector F of subsets S_i to S_n to simulate a universe of elements. 
void read_Fset(ParProg *par, ifstream& file){
    int sets;
    string line,item;
	getline(file,line);
    sets=stoi(line);
    par->ns=sets;
    if (PRINT) cout<< "Sets: "<<sets<<endl;
    for (int j = 0; j < sets; j++)
    {
        getline(file,line);
        istringstream iss(line);
        set<int> act;
        while (getline(iss, item, ' ')) {
                if (CHECK) cout<<item<<" ";
                act.insert(stoi(item));
        }
        if (CHECK) cout<<endl;
        par->F_sets.push_back(act);
    }
    getline(file,line);
    if(PRINT) cout<< "Subsets loaded."<<endl;
}
//Preprocess  of F to obtain universe, cardinality of each element and location.
void analyzeF(ParProg *par){

	unordered_map<int,int> valor_cant; //unordered because chi is not necessarily filled from 0 to n.
	unordered_map<int,vector<int>> inSub;
	auto start = high_resolution_clock::now(); 
	int cont=0; //Needed to mark the index of subset visited.
	
	for(set<int> set_i:par->F_sets){
		for(int elem_j:set_i){
			valor_cant[elem_j]++;
			inSub[elem_j].push_back(cont); //the actual index (in which subset we found elem_j) is added to the map.
			par->chi.insert(elem_j);
		}
		cont++;
	}

	par->chiSize=par->chi.size();
	par->P= new item[par->chiSize];
	int pos=0;
	for(pair<int,int> values:valor_cant){ //iter through pair elements
        (par->P)[pos].value=values.first;
		par->chi_map[values.first]=pos; 
		(par->P)[pos].repetitions=values.second;
		for(int indexSet:inSub[values.first]){
			(par->P)[pos].inSet.push_back(indexSet);
		}
		pos++;
    }    
	int sizeWord=1+(par->chiSize/(8*sizeof(ulong)));
	for (set<int> set_i:par->F_sets)
	{
		ulong* act =new ulong[sizeWord];
		for (int j = 0; j < sizeWord; j++)
		{
			act[j]=0;
		}
		
		for(int elem:set_i){
			setBit64(act,par->chi_map[elem]);
		}
		par->bF_sets.push_back(act);
	}
	auto stop = high_resolution_clock::now();
	if(1) cout<< "Subsets analyzed in "<<(duration_cast<microseconds>(stop-start)).count()<<" microseconds."<<endl;
}
//Check P struct of items, to validate the preprocess.
void checkP(ParProg *par){
	
	for(int cont=0;cont<par->chiSize;cont++){
		item ind=par->P[cont];
		cout <<"#"<<setw(2)<<cont<<" >> ";
		cout<<"Val: "<<setw(2)<<ind.value<< " | ";
		cout<<"Rep.: "<<setw(2)<<ind.repetitions<<" | ";
		cout<<"In Sets:";
		for(int i:ind.inSet){
			cout<<" "<<i;
		}
		cout<<"."<<endl;
	}
	cout<<"\n---------------------------------------------------------"<<endl;
}
//Check cardinality of Chi and his elements.
void checkChi(ParProg *par){
	cout << "|X| = "<<par->chiSize<<endl;
	cout << "Elements of the universe X:"<<endl;
	for(int val:par->chi){
		cout<< val << " ";
	}
	cout<<endl<<"---------------------------"<<endl;
}
//Initialize the memory of bits needed.
void initMemory(ParProg *par){
	par->I.bits=par->chiSize;
	par->I.nW=1+(par->I.bits/(8*sizeof(ulong)));
	par->I.mem=new ulong[par->I.nW];
	for (int i = 0; i < par->I.nW; i++)
	{
		par->I.mem[i]=0;
	}
	par->M.bits=par->F_sets.size();
	par->M.nW=1+(par->M.bits/(8*sizeof(ulong)));
	par->M.mem=new ulong[par->M.nW];
	for (int i = 0; i < par->M.nW; i++)
	{
		par->M.mem[i]=0;
	}
	if(CHECK){
		cout<<"Bits for Elements: "<<par->I.bits<<endl;
		cout<<"Bits for Sets: "<<par->M.bits<<endl;
		checkMemory(par);
	} 
	
}
//Print bits of ulong pointers of elements and subsets covered.
void checkMemory(ParProg *par){
	cout<<"I-Memory: ";
	int cant=par->I.bits;
	int size=W64;
	for (int i = 0; i < par->I.nW; i++)
	{
		if(cant>size){
			printFirstBits(par->I.mem[i],size);
			cant-=W64;
		}else{
			printFirstBits(par->I.mem[i],cant);
		}
	}
    auto start = high_resolution_clock::now(); 
    auto c =countBits(par->I.mem,par->I.bits);
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start);
    cout<<"\ncount: "<<c<<endl;
	cout<<"M-Memory: ";
	cant=par->M.bits;
	for (int i = 0; i < par->M.nW; i++)
	{
		if(cant>size){
			printFirstBits(par->M.mem[i], size);
			cant-=8*sizeof(ulong);
		}else{
			printFirstBits(par->M.mem[i], cant);
		}
	}
    start = high_resolution_clock::now(); 
    c =countBits(par->M.mem,par->M.bits);
    stop = high_resolution_clock::now(); 
    duration = duration_cast<microseconds>(stop - start);
    cout<<"\ncount: "<<c<<endl;
	cout<<"---------------------------------------------------------"<<endl;
}

//Greedy Algorithm
void greedyAlg(ParProg *par){
    if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"            Greedy Algorithm in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"Greedy Search over F_sets."<<endl;
	} 
    set<int> chi(par->chi); 
	map<int,set<int>> candidates;
	int posicionAbs=0;
	for (set<int> each: par->F_sets){
		candidates[posicionAbs]=each;
		posicionAbs++;
	}
	
    while(chi.size()>0){
		if(CHECK) cout<<"Length: "<<candidates.size()<<endl;
        int conjOptimo=0,mayorCover=0,iCover=0;  
		for(pair<int,set<int>> par:candidates){
			iCover=intersect(par.second,chi);
			if(iCover>mayorCover){
				conjOptimo=par.first;
				mayorCover=iCover;
			}
		}
        if(CHECK){
            cout<<"--Best candidate set : "<< conjOptimo+1;
            cout<<"  >Covering: "<<mayorCover<<" new elements."<<endl;
        }
        for (int val:candidates[conjOptimo])
        {
           // if(CHECK) cout<<val << " ";
            chi.erase(val);
		}
		candidates.erase(conjOptimo);
        par->greedy_solution.insert(conjOptimo);
    }
    if(TEST){
		cout<<"Greedy Alg. solution: "<<endl;
        cout<<">>Cardinality: "<<par->greedy_solution.size()<<", Sets: ";
		for(int j:par->greedy_solution){
			cout<<j<<" ";
		}
		cout<<endl;
    }
}
//Preprocess for setCover heuristic.
void preSetCover(ParProg *par){
	sort(par->P,(par->P)+par->I.bits,compareRep); //Sort the P structure of items.
	if(CHECK){
			checkP(par);
		}	
	if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"            (H1) PreSetcover in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
	}
	int pos=0;
	if(CHECK) cout<<"Adding unique chi elements.\n"<<endl;
	while(par->P[pos].repetitions==1){
		int actual_set=par->P[pos].inSet[0]; //select the index of the unique set that contain the element {P[pos]}
		// ***Puede ser reemplazdo por checkBit del elemento en cuestion***
		if(!checkBit(par->I.mem,par->chi_map[par->P[pos].value])){
			set<int> actual(par->F_sets[actual_set]); //load the set from f_sets
			if(CHECK) cout<<"adding set: "<<actual_set<<endl;
			for(int otherVal:actual){
				if(!checkBit(par->I.mem,par->chi_map[otherVal])){ //If the element isn't already included, include it
					if(CHECK) cout<<"    value: "<<otherVal<<endl; 
					setBit64(par->I.mem,par->chi_map[otherVal]); //set the bit of the element to 1.
					par->I.filled++;
				}
			}
			setBit64(par->M.mem,actual_set);//set the bit of the subset to 1.
			par->vidca_solution.insert(actual_set); //add the identificator of the set visited.(0 to k)
		}
		pos++;
		par->last_visited=pos; //update the last visited element.
	}
	if (!pos) cout<<"Ningun elemento agregado en H1"<<endl;
	if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
	}
}

void setCover(ParProg *par){
	int index = par->last_visited;//We scan over P from index m to n  (P[m..n])
	if(PRINT){
		cout<<"            (H2) Setcover in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
	} 
	while(index<par->chiSize && par->I.filled<par->I.bits){
		int value = par->P[index].value;
		//cout<<"REVISANDO: "<<index<<", "<<checkBit(par->I.mem,par->chi_map[value])<<" - "<<par->chi_map[value]<<endl;
		if(!checkBit(par->I.mem,par->chi_map[value])){
            pair<int,int> setSelected = findCandidate(par,index);
            if(CHECK){
                cout<<"Se elige: "<<setSelected.first<<endl;
                cout<<"     Partiendo con: "<<par->I.filled<<" elementos."<<endl;
            }
            
            for (int pos = 0; pos < par->I.nW; pos++)
            {
               par->I.mem[pos]|=(par->bF_sets[setSelected.first])[pos];
            }
            par->I.filled+=setSelected.second;
            if(CHECK) cout<<"     Termina con: "<<par->I.filled<<" elementos."<<endl;

			par->vidca_solution.insert(setSelected.first);
		}
		if(checkBit(par->I.mem,par->chi_map[value])) index++;
	}

}


pair<int,int> findCandidate(ParProg* par, int index){
    int bestCandidate=0;
    int bestSuma=0;
	int value;
	int numSets= (par->P[index].inSet).size();
	int indice=index;
	//cout<<"-------INDEX: "<<index<<"-------- NUMS: "<<numSets<<endl;
	while(indice<par->chiSize && (par->P[indice].inSet).size()==numSets){
		value = par->P[indice].value;
		if(!checkBit(par->I.mem,par->chi_map[value])){
			for(int posibleSet:par->P[indice].inSet){
				ulong* toCheck = par->bF_sets[posibleSet];
				int sumAct=0;
				for (int pos = 0; pos < par->I.nW; pos++)
				{
					ulong e = toCheck[pos]&(~(par->I.mem[pos]));
					sumAct+=countBits(&e,W64);
				}
				if (CHECK) cout<<"> "<<posibleSet<<", suma parcial: "<<sumAct<<endl;
				if(sumAct>bestSuma){
					bestSuma=sumAct;
					bestCandidate=posibleSet;
				}
			}
		}
		indice++;
	}
    
    return pair<int,int> (bestCandidate,bestSuma);
}
void h2(ParProg *par){
	//sort(par->P,(par->P)+par->I.bits,compareRepMost);
	if(CHECK){//CHECK){
			checkP(par);
		}	
	
	if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"            (H2*) in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
	}
	vector<int> leftToCover;
	int pos= par->last_visited;
	while(pos<par->I.bits && par->I.filled<par->I.bits){
		//cout<<">>>VALOR: "<<par->P[pos].value<<"<<<"<<endl;
		if(!checkBit(par->I.mem,par->chi_map[par->P[pos].value])){
			int best=bestCandidate(par,pos);
			if(best){
				int actual_set=par->P[pos].inSet[best]; //select the index of the unique set that contain the element {P[pos]}	
				set<int> actual(par->F_sets[actual_set]); //load the set from f_sets
				if(CHECK) cout<<"adding set: "<<actual_set<<endl;
				for(int otherVal:actual){
					if(!checkBit(par->I.mem,par->chi_map[otherVal])){ //If the element isn't already included, include it
						if(CHECK) cout<<"    value: "<<otherVal<<endl; 
						setBit64(par->I.mem,par->chi_map[otherVal]); //set the bit of the element to 1.
						par->I.filled++;
					}
				}
				setBit64(par->M.mem,actual_set);//set the bit of the subset to 1.
				par->vidca_solution.insert(actual_set); //add the identificator of the set visited.(0 to k)
			}else{
				leftToCover.push_back(pos);
			}
		}
		pos++;
		pos%=par->I.bits;

	}
}

int dist(ParProg *par,int index){
	int mejorDist=0;
	int mejorOpcion=0;
	vector<int> sets(par->P[index].inSet);
	for(int set:sets){
		int dist=0;
		for(int toCheck:par->F_sets[set]){
			if(!checkBit(par->I.mem,par->chi_map[toCheck])){
			//if(!getNum64(par->I.mem,par->chi_map[toCheck],1)){
				dist++;
			}
		}
		if(dist>mejorDist){
			mejorDist=dist;
			mejorOpcion=set;
		}
	}
	return mejorOpcion;
}
int intersect(set<int> conjunto, set<int> chi){
    int cont =0;
    for(int val:conjunto){
            if(chi.count(val)>0){
                cont++;
            }
    }
    return cont;
}
//Used to sort P structure by number of repetitions.
bool compareRep(item a,item b){
	return(a.repetitions<b.repetitions);
}
bool compareRepMost(item a,item b){
	return(a.repetitions>b.repetitions);
}
void printVec(vector<set<int>> conjuntos){
    int contador=0;
    for(set<int> a:conjuntos){
        cout<<"("<<contador<<")"<<" < ";
        for(int b:a){
            cout << b << " ";
        }
        contador++;
        cout<<">"<<endl;
    }
    cout<<"---------------------------"<<endl;
}
void printSuccint(vector<ulong*> bSets, int chiSize){
    int numWord=1+chiSize/64;
	for(ulong* e: bSets){
		for (int i = 0; i < numWord; i++)
		{
			printBitsUlong(e[i]);
			cout<<"-";
		}
		cout<<endl;
	}
}

int countBits(ulong* e, int bits){
	int pos=0;
	int ret=0;
    int len = (int) bits/(W64+1);
	while(pos<=len){
		ret+=__builtin_popcountl(e[pos]);
		pos++;
	}
	return ret;
}
int checkBit(ulong *e, int pos){
	return((((e[pos/W64])>>(W64-1-(pos%W64)))&1));
}
void testMemoryFilled(ParProg *par){
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Universe cardinality: "<<par->chiSize<<endl;
	cout<<"Testing Vidca solution..."<<endl;;
	int dif = par->chiSize-par->I.filled;
	if(dif){
		if(CHECK){
			cout<<"I-Memory: ";
			int cant=par->I.bits;
			int size=8*sizeof(ulong);
			for (int i = 0; i < par->I.nW; i++)
			{
				if(cant>size){
					printFirstBits(par->I.mem[i],size);
					cant-=8*sizeof(ulong);
				}else{
					printFirstBits(par->I.mem[i],cant);
				}
			}
		}
		cout<<"\n[Not completed] --"<<dif<<"-- elements left to cover."<<endl;
	}else{
		cout<<"All elements were covered in (H1) and (H2)."<<endl;
	}
}