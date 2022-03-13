#include <bits/stdc++.h>
#include "hgraph.h"
#include<windows.h>



using namespace std;

//class type global declaration

// global parameters
int SEED =3;
double RATIO =0.5;
bool MODE;
int initcut;




int main()
{
    //set for time counter
    LARGE_INTEGER nFreq;
	LARGE_INTEGER nBeginTime;
	LARGE_INTEGER nEndTime;
	QueryPerformanceFrequency(&nFreq);
    int opt;
    cout << "1 for speedtest mode, 2 for visualization mode:  " << endl;
    cin >> opt;
    if (opt == 1){
        MODE = true;
    }
    else if (opt == 2){
        MODE = false;
    }
    else{
        cout << "Invalid input, program exit." << endl;
        exit(0);
    }
//    fstream file;
//    char buffer[256];
//    file.open(inputname, ios::in);
//    if(!file){
//        cout << "Failed to open the file." << endl;
//        exit(0);
//    }
    char *inputname = "ibm09.net";
    char name[]="partition_result.txt";
    double time;
    //init the graph parameters
    parthgraph hg;
    QueryPerformanceCounter(&nBeginTime);  // start timer
    hg.getgraph(inputname);
    hg.initgains();
    hg.part();
    QueryPerformanceCounter(&nEndTime);    // stop timer
    time=(double)(nEndTime.QuadPart-nBeginTime.QuadPart)/(double)nFreq.QuadPart;

    //output file setting
    char *fname;
    fname = name;

    ofstream myresult(fname,ios::app);
    if (!myresult)
    {
        cout<<"error opening file ";
        cin>>name;
        exit(0);
    }
    else{
        char string[100];
        hg.print(string);
        myresult << "Benchmark: " << inputname << "  ";
        myresult<<"Seed: "<< hg.currTime() << "(System Time in seconds)" <<"  ";
        myresult << "Init Cutset: " << initcut << "  " ;
        myresult<<"Final Cutset: "<<hg.cutset() << "  ";
        myresult << "Improvement percentage: " <<  ((double)initcut - (double)hg.cutset()) / (double)initcut * 100 << "%   ";
        myresult<<"Time: "<< time*1000 << "ms" <<"  ";
        myresult<<string;
        myresult<<'\n';
    }


    return 0;
}
