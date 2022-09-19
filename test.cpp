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
#include <map>
using namespace std;
using namespace cds;

#define PRINT 1
#define TEST 1
#define REP 1
#define maxCantElem 15
#define minCantElem 5
#define maxVal 40

int countBits(ulong* e, int bits);
int checkBit(ulong *e, int pos);
void printVec(vector<set<int>> vec);

int main(int argc, char** argv){ // maximum number of stars to distribute
    ulong* arr= new ulong[2];
    cout<<"arr Bits: "<<countBits(arr,128)<<endl;
    arr[0]=0;
    arr[1]=3;
    cout<<"arr Bits: "<<countBits(arr,128)<<endl;
    vector<ulong*> data;
    data.push_back(arr);
    for(ulong* a:data){
        cout<<"arr: "<<countBits(a,128)<<endl;
    }
    cout<<"65/64= "<<65/64<<endl;
	return 0;
}
void printVec(vector<set<int>> vec){
    for(set<int> a:vec){
        for(int b:a){
            cout<<b<<" ";
        }
        cout<<endl;
    }
}
int countBits(ulong* e, int bits){
	int pos=0;
	int ret=0;
	while(pos<=(int)(bits/(W64+1))){
		ret+=__builtin_popcountl(e[pos]);
		pos++;
	}
	return ret;
}
int checkBit(ulong *e, int pos){
	return((((e[pos/W64])>>(W64-1-(pos%W64)))&1));
}
