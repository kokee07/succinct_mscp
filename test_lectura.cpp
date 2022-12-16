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
int main(int argc, char** argv){
	if(argc !=2){
		cout << "FILE: Execution Error! call: ./test <fileName>" << endl;
		exit(EXIT_FAILURE);
	}
	
	char FILENAME[40];
	strcpy(FILENAME,argv[1]);

	cout << "Reading file..." <<FILENAME<< endl;
	char prefix[80];
	strcpy(prefix,FILENAME);
	strcat(prefix,"_results_HEU");
	FILE *fp = fopen(prefix, "w" );
    ifstream file(FILENAME);
	auto ss = high_resolution_clock::now(); 
    string line,item;
    int sets, col, row;
    // getline(file>>std::ws,line);
    // cout<<"L: "<<line<<endl;
    // sets = stoi(line);
    cout << "SETS: "<<sets<<endl;
    int count = 0;
    while (!file.eof()) {
        if (count > 10) break;
        getline(file,line);
        cout<<"init: "<<line<<endl;
        istringstream iss(line);
        while (iss>>item){
            cout<<"item: "<<item<<" "<<endl;
        }
        cout<<endl;
        // istringstream iss(line);

        // while (iss>>item) {
        //     cout<<item<<" ";
        // }
        // cout<<endl;
        count++;
        // if(file.eof()) break;
        // int REP=stoi(line);
        // REP = 1;
        // int sets;
        // getline(file,line);
        // sets=stoi(line);
        // cout<< "Sets: "<<sets<<endl;
        // for (int j = 0; j < sets; j++)
        // {
        //     getline(file,line);
        //     istringstream iss(line);
        //     set<int> act;
        //     while (getline(iss, item, ' ')) {
        //             cout<<item<<" ";
        //             act.insert(stoi(item));
        //     }
        //     cout<<endl;
        // }
    }
    cout<<"????? "<<count<<endl;
	file.close();
    fclose (fp);
	auto ee = high_resolution_clock::now(); 
	auto dd = duration_cast<microseconds>(ee - ss);
	cout<<"Total: "<<dd.count()<<endl;
	return 0;
}