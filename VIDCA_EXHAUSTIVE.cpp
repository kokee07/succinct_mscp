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

using namespace std;
using namespace std::chrono; 
using namespace cds;

#define PRINT 0
#define CHECK 0
#define TEST 0
#define REP 1
#define maxCantElem 15
#define minCantElem 4
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
	vector<set<int>> F_sets; //sets of elements.
    set<int> chi; //Universe elements
	map<int,int> chi_map;//Map each element of chi from 0 to n. (index in P_structure to the real value)
	set<int> exha_solution;
    set<int> greedy_solution;
    set<int> vidca_solution;
	int maxLength;
	int last_visited;
	memory M;//bits for subsets already covered.
	memory I;//bits for elements of Chi already covered.
	item * P; //elements of chi, their cardinality and which subset included them.
	ulong ns;
	ulong MAX_val;
} ParProg;

//Generators 
void genF(ParProg *par);
void analyzeF(ParProg *par);
void initMemory(ParProg *par);
//Check
void checkChi(ParProg *par);
void checkP(ParProg *par);
void checkMemory(ParProg *par);
//Main methods
void exhaustiveSearch(ParProg *par);
void preSetCover(ParProg *par);
void setCover(ParProg *par);
void greedyAlg(ParProg *par);
//Subrutines.
bool compareRep(item a,item b);
int dist(ParProg *par,int index);
void printVec(vector<set<int>> conjuntos);
int intersect(set<int> conjunto, set<int> chi);
//Test Methods
void testMemoryFilled(ParProg *par);
void testSolutionSets(ParProg *par);


int main(int argc, char** argv){
	if(argc != 3){
		cout << "Execution Error! call: ./vidca <n_sets> <MAX_val>" << endl;
		exit(EXIT_FAILURE);
	}
	
	int n = atoi(argv[1]);
	int Mval = atoi(argv[2]);
	default_random_engine generator(time(NULL));
	if(PRINT){
		cout << "Parameters..." << endl;
		cout << "# subsets to generate = " << n << endl;
		cout << "Max value of elements = " << Mval << endl;
	}
	ulong durPrep=0;
	ulong durTotalVidca=0;
	ulong durTotalGreedy=0;
    ulong durTotalExhaustive=0;
	ulong cardinalidadVidca=0;
	ulong cardinalidadGreedy=0;
    ulong cardinalidadExhaustive=0;
	if(PRINT) cout<<"Initializing repetitions..."<<endl;
	for (int k = 0; k < REP; k++)
	{	
		if(PRINT) cout<<"Iter #"<<k<<"."<<endl;
		ParProg *par = new ParProg();
		par->ns=n;
		par->MAX_val=Mval;
		genF(par);
		auto start = high_resolution_clock::now(); 
		analyzeF(par);
		auto stop = high_resolution_clock::now(); 
		auto duration = duration_cast<microseconds>(stop - start); 
		durPrep+=duration.count();
		if(PRINT){
			cout<<"Checking elements in P_structure:\n"<<endl;
			checkP(par);
			cout<<"Checking chi, Universe of elements:\n"<<endl;
			checkChi(par);
			cout<<"Checking F_sets generated:\n"<<endl;
			printVec(par->F_sets);
		}
		initMemory(par);
		//Exhaustive search for the cardinality best solution (minimun set cover)
		//Or at least one of the best solutions if there are more than one.
		start = high_resolution_clock::now(); 
		//exhaustiveSearch(par);
		stop = high_resolution_clock::now(); 
		duration = duration_cast<microseconds>(stop - start); 
		cardinalidadExhaustive+=par->exha_solution.size();
		durTotalExhaustive+=duration.count();
        //Call of Greedy alg.
        start = high_resolution_clock::now(); 
		greedyAlg(par);
		stop = high_resolution_clock::now(); 
		duration = duration_cast<microseconds>(stop - start); 
		durTotalGreedy+=duration.count();
		cardinalidadGreedy+=par->greedy_solution.size();
		//Call of SetCover methods (Heuristics H1 and H2)
		if(CHECK){
			checkP(par);
		}	
		start = high_resolution_clock::now(); 
		preSetCover(par);
		stop = high_resolution_clock::now(); 
		duration = duration_cast<microseconds>(stop - start);
		durTotalVidca+=duration.count();
		if(CHECK){
			checkMemory(par);
		}
		start = high_resolution_clock::now(); 
		if(par->I.filled!=par->I.bits){
			setCover(par);
		}
		stop = high_resolution_clock::now(); 
		duration = duration_cast<microseconds>(stop - start);
		if(CHECK){
			checkMemory(par);
		} 
		durTotalVidca+=duration.count();
		cardinalidadVidca+=par->vidca_solution.size();

		if(TEST){
			testSolutionSets(par);
		}
		if(PRINT) cout<<"Finalized #"<<k<<endl;
	}
	cout << "################## " << endl;
	cout << "Cycles"<<endl;
    cout << "Time taken by "<<REP<<" reps.: "
         << durTotalExhaustive+durTotalVidca+durTotalGreedy<< " microseconds" << endl; 
	cout << "Average time per cycle of Exhaustive: "<<durTotalExhaustive/REP<< " microseconds"<<endl;
	cout << "Average cardinality per cycle Exhaustive: "<<cardinalidadExhaustive/REP<< " subsets."<<endl;
	cout << "Average time per cycle of Greedy: "<<durTotalGreedy/REP<< " microseconds"<<endl;
	cout << "Average cardinality per cycle Greedy: "<<cardinalidadGreedy/REP<< " subsets."<<endl;
    cout << "Average time per cycle of Vidca: "<<durTotalVidca/REP<< " microseconds"<<endl;
	cout << "Average cardinality per cycle Vidca: "<<cardinalidadVidca/REP<< " subsets."<<endl;
	cout << "################## " << endl;
	cout << "Time taken by Preprocess: "<<durPrep/REP<<" microseconds."<<endl;
	cout << "################## " << endl;
	return 0;
}

//Generate vector F of subsets S_i to S_n to simulate a universe of elements. 
void genF(ParProg *par){
    int cantElem,val;
    for (ulong i = 0; i < par->ns; i++)
    {   
		default_random_engine generator(time(NULL));
        set<int> actual;
        cantElem=rand()%(maxCantElem-minCantElem + 1) + minCantElem;
		int mean=rand()%par->MAX_val; //random mean for each set
		normal_distribution<double> distribution(mean,(par->MAX_val)/8); //normal distribution for each set
		while((int)actual.size()<cantElem){
			val= ((int)distribution(generator))%(par->MAX_val);
			actual.insert(abs(val));
		}
		par->F_sets.push_back(actual);
    }
    if(PRINT) cout<< "Subsets Generated."<<endl;
}
//Preprocess  of F to obtain universe, cardinality of each element and location.
void analyzeF(ParProg *par){

	unordered_map<int,int> valor_cant;
	unordered_map<int,vector<int>> inSub;
	int maxSetSize=par->F_sets[0].size();
	int cont=0; //Needed to mark the index of subset visited.
	for(set<int> set_i:par->F_sets){
		if((int)set_i.size()>maxSetSize){
			maxSetSize=(int)set_i.size();
		}
		for(int elem_j:set_i){
			valor_cant[elem_j]++;
			inSub[elem_j].push_back(cont); //the actual index (in which subset we found elem_j) is added to the map.
			par->chi.insert(elem_j);
		}
		cont++;
	}
	par->maxLength=maxSetSize;
	par->P= new item[valor_cant.size()];
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
	if(PRINT) cout<< "Subsets analyzed."<<endl;
}
//Check P struct of items, to validate the preprocess.
void checkP(ParProg *par){
	
	for(ulong cont=0;cont<par->chi.size();cont++){
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
	cout << "|X| = "<<par->chi.size()<<endl;
	cout << "Elements of the universe X:"<<endl;
	for(int val:par->chi){
		cout<< val << " ";
	}
	cout<<endl<<"---------------------------"<<endl;
}
//Initialize the memory of bits needed.
void initMemory(ParProg *par){
	par->I.bits=par->chi.size();
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
	cout<<"\nM-Memory: ";
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
	cout<<"\n---------------------------------------------------------"<<endl;
}
//Exhaustive search over the F_sets to find the best solution
void exhaustiveSearch(ParProg *par){
	set<int> bestSol;
	int numBestSol=par->ns;
	bool solved=false;
	//Index is calculated as the minimun numbers of set needed to cover all the elements of chi.
	int index=par->chi.size()/par->maxLength + (par->chi.size() % par->maxLength != 0);
    if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"            Exhaustive Search in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"Exhaustive Search over F_sets."<<endl;
	} 
	if(CHECK){
		cout<<"Minimun number of sets needed: "<<index<<endl;
		cout<<"To complete: "<<par->chi.size()<<endl;
	}
	//We iter over the possible length of sets combination
	while(index <=(int)par->ns && !solved){
        vector<bool> v(par->ns);
        fill(v.begin(), v.begin() + index, true);
		//Gen every possible combination of index-length sets.
        do {
            set<int> candidatos;
            for (int i = 0; i < (int)par->ns; ++i) {
                if (v[i]) {
                    candidatos.insert(i);
                }
			}
			//Makes the union of elements of every set involved.
			set<int> unionActual;
			for (int setActual:candidatos){
				for(int elemActual:par->F_sets[setActual]){
					unionActual.insert(elemActual);
				}
			}
			if(CHECK){
				cout<<"Candidatos:"<<endl;
				cout<<"< ";
				for(int c:candidatos){
					cout<<c<<" ";
				}
				cout<<">, Con: "<<unionActual.size()<<" elementos."<<endl;
			}
			//if all elements of chi are covered.
			if(unionActual.size()==par->chi.size()){
				//if is better than the actual bestSol.
		
					bestSol=candidatos;
					numBestSol=candidatos.size();
					
					solved=true;
			}
        } while (!solved && prev_permutation(v.begin(), v.end()));
		index++;
    }
	//Update exhaustive solution in par struct.
	par->exha_solution=bestSol;
	if(TEST){
		cout<<"Exhaustive Solution: "<<endl;
		cout<<">>Cardinality: "<<numBestSol<<", Sets: ";
		for(int j:bestSol){
			cout<<j<<" ";
		}
		cout<<endl;
	}
}
void greedyAlg(ParProg *par){
    if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"            Greedy Algorithm in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"Greedy Search over F_sets."<<endl;
	} 
    set<int> chi(par->chi); 
    vector<set<int>> candidates(par->F_sets);
    while(chi.size()>0){
		if(CHECK) cout<<"Length: "<<candidates.size()<<endl;
        int conjOptimo=0,mayorCover=0,iCover=0;  
		//ToSolve: Iteraciones demás al releer sets seleccionados.
        for (int i = 0; i < (int)candidates.size(); i++)
        {
            iCover=intersect(candidates[i],chi);
            if(iCover>mayorCover){
                conjOptimo=i;
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
		set<int> actual(par->F_sets[actual_set]); //load the set from f_sets
		if(CHECK) cout<<"adding set: "<<actual_set<<endl;
		for(int otherVal:actual){
			if(!getNum64(par->I.mem,par->chi_map[otherVal],1)){ //If the element isn't already included, include it
				if(CHECK) cout<<"    value: "<<otherVal<<endl; 
				setBit64(par->I.mem,par->chi_map[otherVal]); //set the bit of the element to 1.
				par->I.filled++;
			}
		}
		setBit64(par->M.mem,actual_set);//set the bit of the subset to 1.
		par->vidca_solution.insert(actual_set); //add the identificator of the set visited.(0 to k)
		
		
		pos++;
		par->last_visited=pos; //update the last visited element.
	}
	if(TEST){
		if(par->last_visited){
			cout<< "(H1) Elements of P already visited: "<<par->last_visited<<endl;
			cout<< "preSetCover partial solution sets: "<<endl;
			cout<<">>Cardinality: "<<par->vidca_solution.size()<<", Set(s): ";
			for(int k:par->vidca_solution){
				cout<<k<<" ";
			}
			cout<<endl;
		}else{
			cout<< "(H1) No element was added."<<endl;
		}
		testMemoryFilled(par);
	}
	if(PRINT){
		cout<<"---------------------------------------------------------"<<endl;
	}
}

void setCover(ParProg *par){
	int index = par->last_visited;//We scan from P from index m to n  (P[m..n])
	if(PRINT){
		cout<<"            (H2) Setcover in progress."<<endl;
		cout<<"---------------------------------------------------------"<<endl;
	} 
	while(par->I.filled<par->I.bits){
		int value = par->P[index].value;
		if(!getNum64(par->I.mem,par->chi_map[value],1)){
			int setSelected=dist(par,index);
			if(CHECK) cout<<"### Best set option: "<<setSelected<<endl;
			for(int otherVal:par->F_sets[setSelected]){
				if(!getNum64(par->I.mem,par->chi_map[otherVal],1)){
					if(CHECK) cout<<"    value: "<<otherVal<<endl; 
					setBit64(par->I.mem,par->chi_map[otherVal]); //set the bit of each element in the same subset.
					par->I.filled++;
				}
			}
			setBit64(par->M.mem,setSelected);
			par->vidca_solution.insert(setSelected);
		}
		index++;
	}
	
	if(TEST){
		cout<<"(H2) setCover Solution sets: "<<endl;
		cout<<">>Cardinality: "<<par->vidca_solution.size()<<", Set(s): ";
		for(int k:par->vidca_solution){
			cout<<k<<" ";
		}
		cout<<endl;
		testMemoryFilled(par);
		cout<<"---------------------------------------------------------"<<endl;
	}

}
int dist(ParProg *par,int index){
	int mejorDist=0;
	int mejorOpcion=0;
	vector<int> sets(par->P[index].inSet);
	for(int set:sets){
		int dist=0;
		for(int toCheck:par->F_sets[set]){
			if(!getNum64(par->I.mem,par->chi_map[toCheck],1)){
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
void testMemoryFilled(ParProg *par){
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Testing memory..."<<endl;;
	int dif = par->I.bits-par->I.filled;
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
		cout<<"\n--"<<dif<<"-- elements left to cover."<<endl;
	}else{
		cout<<"All elements are covered in (H1) and (H2)."<<endl;
	}
}
void testSolutionSets(ParProg *par){
	bool equal=true;
	cout<<"Testing solutions..."<<endl;
    cout<<"\nComparing Exhaustive and Greedy solutions."<<endl;
    if(par->exha_solution.size()==par->greedy_solution.size()){
        for (int j:par->exha_solution){
            if(!(find(par->greedy_solution.begin(),par->greedy_solution.end(),j)!=par->greedy_solution.end())){
                cout<<"Both solutions sets are different but of equal cardinality."<<endl;
				equal=false;
                break;
            }
        }
		if(equal){
			cout<<"Both solutions sets are equal."<<endl;
		}
    }else{
		cout<<"The cardinality of the solutions set is different."<<endl;
	}
    cout<<"\nComparing Exhaustive and Vidca solutions."<<endl;
    if(par->exha_solution.size()==par->vidca_solution.size()){
        for (int j:par->exha_solution){
            if(!(find(par->vidca_solution.begin(),par->vidca_solution.end(),j)!=par->vidca_solution.end())){
                cout<<"Both solutions sets are different but of equal cardinality."<<endl;
				equal=false;
                break;
            }
        }
		if(equal){
			cout<<"Both solutions sets are equal."<<endl;
		}
    }else{
		cout<<"The cardinality of the solutions set is different."<<endl;
	}
    cout<<"\nComparing Greedy and Vidca solutions."<<endl;
    if(par->greedy_solution.size()==par->vidca_solution.size()){
        for (int j:par->greedy_solution){
            if(!(find(par->vidca_solution.begin(),par->vidca_solution.end(),j)!=par->vidca_solution.end())){
                cout<<"Both solutions sets are different but of equal cardinality."<<endl;
				equal=false;
                break;
            }
        }
		if(equal){
			cout<<"Both solutions sets are equal."<<endl;
		}
    }else{
		cout<<"The cardinality of the solutions set is different."<<endl;
	}
	cout<<"---------------------------------------------------------"<<endl;
}