/*
 * C++ Program to Implement Adjacency List
 */
#include <iostream>
#include<fstream>
#include<cstdlib>
#include<cstring>
#include<string>
#include<vector>
#include"graph_es.h"
#include<time.h>
#include<set>
#include <algorithm>
#include <map>
#include <random>
#include <sys/time.h>
using namespace std;




//
//====================================================get_wall_time====================================================
//
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

//
//====================================================main====================================================
//
int main(int argc,char **argv)
{
	cout.setf(ios::fixed,ios::floatfield);//Non-scientific notation
	if (argc!=6)
	{cout<<"please input five parameters: infile, motiffile, outfile, delta, edge sampling P!"<<endl; exit(1);}
	string infile=argv[1];
	string motiffile=argv[2];
	string outfile=argv[3];
	int delta=atoi(argv[4]);
	float P=atof(argv[5]);

	cout<<"input file: "<<infile<<endl;
	cout<<"motif file: "<<motiffile<<endl;
	cout<<"output file: "<<outfile<<endl;
	cout<<"delta: "<<delta<<endl;
	cout<<"P: "<<P<<endl;

	ifstream in(infile.c_str());
	if(!in){cout<<"open infile failed\n";exit(1);}
	in.close();
	in.open(motiffile);
	if(!in){cout<<"open motiffile failed\n";exit(1);}
	in.close();
	ofstream out(outfile.c_str(),ios::app);	
	out.setf(ios::fixed,ios::floatfield);//Non-scientific notation
	if(!out){cout<<"open outfile failed\n";exit(1);}
	out.close();
	int found=infile.find_last_of("/");
	string s(infile.substr(found+1,3));
	cout<<"name of prefile"<<s<<endl;

//
//====================================================read the graph====================================================
//
	Graph G;
	//Graph *g=&G;
	char buf[100];
	const char* d=" ";
	long a[3];
	double start,time0,time1;
	in.open(infile.c_str());
	while (in.getline(buf,100)) // 
	{ 
		char* p=strtok(buf,d);
		int i=0;
		while (p){			
			a[i]=atol(p); // 
			p=strtok(NULL,d);
			i=i+1;	
		}
		//cout<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
		delete p;
		if (a[0] != a[1]) { 
			//G.edges.push_back({a[0], a[1], a[2], 0L});
			G.addEdge(a[0],a[1],a[2],0L);
		}
	}
	G.Initialize();
	/*for(int i=0;i<G.edges_.size();i++)
	{
		cout<<G.edges_[i].src<<" "<<G.edges_[i].dst<<" "<<G.edges_[i].tim<<" "<<G.edges_[i].id<<endl;
	}*/
	cout<<"Finish reading "<<G.getEn()<<" edges and the last edge is "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
	/*G.Printgraph();
	for (int i=0;i<G.getVn();i++)
	{
		sum_d=sum_d+G.nodeDegree(i);
	}
	cout<<"Sum of in and out degree: "<<sum_d<<endl;
	cout<<"the number of edges in hash table: "<<G.getHashn()<<endl;
	G.printGraph();*/
	in.close();

//
//====================================================sampling====================================================
//
	int seeds[10]={5,10,15,20,25,30,35,40,45,50};
	for (int jj=0;jj<10;jj++)
	{
		int seed_p=seeds[jj];
		double rn;
		vector<long> ids;
		default_random_engine e(seed_p);
		uniform_real_distribution<double> u(0.0,1.0);
		for(long i=0;i<G.getEn();i++)
		{
			rn=u(e);
			//cout<<rn<<endl;
			if(rn<=P){
			/*if(G.edges_[i].tim==1197014240)
			{*/
				//cout<<i<<" ";
				ids.push_back(i);
			}
		}
		cout<<"number of sampling:"<<ids.size()<<endl;
		//cout<<G.getHashn()<<endl;
	
//
//====================================================read a motif====================================================
//
		/*int mid[5]={66,12,35,24,7};
		for (int ii=0;ii<5;ii++)
		{*/
			/*string motiffile;
			motiffile="/home/jingjing/Cmotif/basic/motifs/motif"+to_string(mid[ii])+".txt";
			string outfile="/home/jingjing/Cmotif/basic/results/ES/Table1/"+s+to_string(mid[ii])+".txt";*/
			//cout<<motiffile<<endl;
			
			G.cle_motif();
			in.open(motiffile);
			int b[3];
			while (in.getline(buf,100)) // 
			{ 
				char* p=strtok(buf,d);
				int i=0;
				while (p){			
					b[i]=atol(p); // 
					p=strtok(NULL,d);
					i=i+1;	
				}
				delete p;
				if (b[0] != b[1]) { 
					cout<<b[0]<<" "<<b[1]<<endl;
					G.motif(b[0], b[1]);
				}
				else{cout<<"the motif has self-loops"<<b[0]<<" "<<b[1]<<endl; exit(1);}
			}
			in.close();
			cout<<"finishing reading the motif"<<endl;
			G.setVm_();
			start=get_wall_time(); 
			//cout<<"begin"<<endl;
			double T=G.ExactCountMotifs(delta,ids)/P;
			time0=get_wall_time()-start;

			out.open(outfile.c_str(),ios::app);
			//cout<<"T="<<T<<" T/"<<G.edgesM_.size()<<"="<<T/G.edgesM_.size()<<" "<<time0<<endl;
			out<<T/G.edgesM_.size()<<" "<<time0<<endl;
			out.close();
		//}
	}
}

