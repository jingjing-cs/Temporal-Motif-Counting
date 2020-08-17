/*
 * C++ Program to Implement Adjacency List
 */
#include <iostream>
#include<fstream>
#include<cstdlib>
#include<cstring>
#include<string>
#include<vector>
#include"graph_ews.h"
#include<time.h>
#include<set>
#include <algorithm>
#include <map>
#include <random>
#include <sys/time.h>
using namespace std;

//declear



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
//====================================================Function declaration====================================================
//
//for triangle
long Cmotif_T2(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T3(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T0(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T1(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T10(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T11(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T8(Graph *G,long id, int delta,float Q,int seed,Graph *M);
long Cmotif_T9(Graph *G,long id, int delta,float Q,int seed,Graph *M);

//for wedge
long Cmotif_0_2(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_1_3(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_4_6(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_5_7(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_8_9_16_17(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_10_11(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_12_13(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_14_15_22_23(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_18_19(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);
long Cmotif_20_21(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M);


//long Cpair(Graph* G, long id,int delta, float Q, int seed_q, Graph* M)
//{
//	long T;
//	T=G->Cmotif_Pair(id,delta,Q,seed_q,M);
//	return T;		//T=qp(&G,ids[i],delta,Q,seed_q,&M);
//}


//
//====================================================main====================================================
//
int main(int argc,char **argv)
{
	cout.setf(ios::fixed,ios::floatfield);//Non-scientific notation
	if (argc!=8)
	{cout<<"please input seven parameters: infile, motiffile, outfile, delta, edge sampling P, wedge sampling Q,seed for P!"<<endl; exit(1);}
	string infile=argv[1];
	string motiffile=argv[2];
	string outfile=argv[3];
	int delta=atoi(argv[4]);
	float P=atof(argv[5]);
	float Q=atof(argv[6]);
	int seed_p=atoi(argv[7]);
	int seed_q=seed_p+2;

	cout<<"input file: "<<infile<<endl;
	cout<<"motif file: "<<motiffile<<endl;
	cout<<"output file: "<<outfile<<endl;
	cout<<"delta: "<<delta<<endl;
	cout<<"P: "<<P<<". Q: "<<Q<<endl;
	cout<<"seed_p: "<<seed_p<<". seed_q: "<<seed_q<<endl;

	ifstream in(infile.c_str());
	if(!in){cout<<"open infile failed\n";exit(1);}
	in.close();
	in.open(motiffile.c_str());
	if(!in){cout<<"open motiffile failed\n";exit(1);}
	in.close();
	ofstream out(outfile.c_str(),ios::app);	
	out.setf(ios::fixed,ios::floatfield);//Non-scientific notation
	if(!out){cout<<"open outfile failed\n";exit(1);}
	out.close();

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
	double rn;
	vector<long> ids;
	default_random_engine e(seed_p);
	uniform_real_distribution<double> u(0.0,1.0);
	for(long i=0;i<G.getEn();i++)
	{
		rn=u(e);
		//cout<<rn<<endl;
		if(rn<=P)
		{
			ids.push_back(i);
		}
	}
	cout<<"number of sampling:"<<ids.size()<<endl;
	
//
//====================================================read a motif====================================================
//
	int id, classMotif,center;
	long T=0, sumT=0;
	Graph M;
	in.open(motiffile.c_str());
	while (in.getline(buf,100)) // 
	{ 
		char* p=strtok(buf,d);
		int i=0;
		while (p){			
			a[i]=atol(p); // 
			p=strtok(NULL,d);
			i=i+1;	
		}
		delete p;
		if (a[0] != a[1]) { 
			
			M.addEdge(a[0],a[1],a[2],0L);
		}
		else{cout<<"the motif has self-loops"<<a[0]<<" "<<a[1]<<endl; exit(1);}
	}
	M.Initialize();

//
//====================================================counting====================================================
//
	start=get_wall_time(); 
	string MM;
	int Vm=M.getVn();
	switch(Vm)
	{
	case 2:
		classMotif=1;
		time0=get_wall_time()-start;
		cout<<"classMotif="<<classMotif<<endl;
		start=get_wall_time();
		for (int i=0;i<ids.size();i++)
		{
			T=G.Cmotif_Pair(ids[i],delta,Q,seed_q,&M);			
			sumT=sumT+T;
		}
		time1=get_wall_time()-start;
		break;
	case 3:
		classMotif=3;
		for (int i=0;i<Vm;i++)
		{
			if(M.nodeDegree(i)==3)
			{
				center=i;
				classMotif=2;
				break;
			}	
		}
		long (*qq)(Graph *,long, int,float,int,Graph *);
		long (*pp) (Graph *,long, int, int,int order[][3],float,int, Graph *);
		if(classMotif==3)
		{
			id=M.cal_Tid();
			time0=get_wall_time()-start;
			switch(id)
			{
			case 0:
				MM="M23";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T0;
				break;
			case 1:
				MM="M24";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T1;
				break;
			case 2:
				MM="M13";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T2;
				break;
			case 3:
				MM="M14";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T3;
				break;
			case 8:
				MM="M45";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T8;
				break;
			case 9:
				MM="M46";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T9;
				break;
			case 10:
				MM="M35";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T10;
				break;
			case 11:
				MM="M36";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				qq=Cmotif_T11;
				break;
			}
			start=get_wall_time();
			for (int i=0;i<ids.size();i++)
			{
				T=qq(&G,ids[i],delta,Q,seed_q,&M);
				sumT=sumT+T;
			}
			time1=get_wall_time()-start;
			
		}
		else
		{
			int order[3][3];
			id=M.cal_Wid(center,order);
			time0=get_wall_time()-start;
			//cout<<"id"<<id<<endl;
			time0=get_wall_time()-start;
			switch(id)
			{
			case 0:
				MM="M63";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_0_2;
				break;
			case 2:
				MM="M53";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_0_2;
				break;
			case 1:
				MM="M64";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_1_3;
				break;
			case 3:
				MM="M54";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_1_3;
				break;
			case 4:
				MM="M55";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_4_6;
				break;
			case 6:
				MM="M65";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_4_6;
				break;
			case 5:
				MM="M56";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_5_7;
				break;
			case 7:
				MM="M66";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_5_7;
				break;
			case 8:
				MM="M41";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_8_9_16_17;
				break;
			case 9:
				MM="M42";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_8_9_16_17;
				break;
			case 16:	
				MM="M43";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_8_9_16_17;
				break;
			case 17:
				MM="M44";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_8_9_16_17;
				break;
			case 10:
				MM="M31";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_10_11;
				break;
			case 11:
				MM="M23";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_10_11;
				break;
			case 12:
				MM="M22";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_12_13;
				break;
			case 13:
				MM="M21";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_12_13;
				break;
			case 14:
				MM="M12";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_14_15_22_23;
				break;
			case 15:
				MM="M11";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_14_15_22_23;
				break;
			case 22:	
				MM="M15";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_14_15_22_23;
				break;
			case 23:
				MM="M16";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_14_15_22_23;
				break;
			case 18:
				MM="M33";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_18_19;
				break;
			case 19:
				MM="M34";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_18_19;
				break;
			case 20:	
				MM="M25";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_20_21;
				break;		
			case 21:
				MM="M26";
				cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
				pp=Cmotif_20_21;
				break;		
			}
			start=get_wall_time();
			for (int i=0;i<ids.size();i++)
			{
				T=pp(&G,ids[i],delta,center,order,Q,seed_q,&M);
				sumT=sumT+T;
			}
			time1=get_wall_time()-start;
		}
		break;
	case 4:
		classMotif=4;
		if(M.edges[1].src==M.edges[0].dest)
		{
			MM="M7";
			time0=get_wall_time()-start;
			cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
			start=get_wall_time();
			for (int i=0;i<ids.size();i++)
			{
				T=G.Cmotif_F7(ids[i], delta, &M);
				sumT=sumT+T;
			}
			time1=get_wall_time()-start;
		}
		else
		{
			MM="M8";
			time0=get_wall_time()-start;
			cout<<"classMotif="<<classMotif<<" Mid="<<MM<<endl;
			start=get_wall_time();
			for (int i=0;i<ids.size();i++)
			{
				T=G.Cmotif_F8(ids[i], delta, &M);
				sumT=sumT+T;
			}
			time1=get_wall_time()-start;
		}
		break;
	}
	/*cout<<"finish reading the motif, classMotif="<<classMotif<<endl;
	cout<<"id"<<id<<endl;*/

	out.open(outfile.c_str(),ios::app);
	double sum=double(sumT)/double(M.getEn())/(P*Q);
	double time_=time0+time1;
	cout<<"sumT: "<<sumT<<" sum:"<<sum<<endl;
	cout<<"time0:"<<time0<<" time1:"<<time1<<" time:"<<time_<<endl;
	out<<sum<<" "<<time_<<endl;
	out.close();




	/*	for(vector<pair<int,int>>::iterator it=edgesM.begin();it!=edgesM.end();it++)
	{
		cout<<(*it).first<<" "<<(*it).second<<endl;
	}*/

	/*double time_c=get_wall_time();
	for (int j=0;j<ids.size();j++)
	{
		
	}
	*/	

}


long Cmotif_T2(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	//M13: i->j,k->j,i->k;
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;

	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M13: i->j,k->j,i->k; find the third edge.check the second edge;O+   (0,2,1)
	//M13: i->j,k->j,i->k; find the third edge. check the first edge. I+ (1,2,0)
	//M13: i->j,k->j,i->k; find the first edge. O- check the second edge.(2 0 1)
		int order0[3]={0,2,1};
		int order1[3]={1,2,0};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri01(id,delta,src,dest,order0,Q,seed,M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri11(id,delta,src,dest,order1,Q,seed,M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id,delta,src,dest,order2,Q,seed,M);
		sum=sum+t;
	}
	else
	{
		//M13: i->j,k->j,i->k; find the second edge.I+ check the third edge.(0,1,2)
		//M13: i->j,k->j,i->k;  find the first edge. I- check the third edge.(1,0,2)
		//M13: i->j,k->j,i->k; find the second edge O-. check the first edge.(2,1,0)
		int order0[3]={0,1,2};
		int order1[3]={1,0,2};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri11(id,delta,dest,src,order0,Q,seed,M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri10(id,delta,dest,src,order1,Q,seed,M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id,delta,dest,src,order2,Q,seed,M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T3(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M14: i->j,k->j,k->i;  find the third edge.check the second edge; I+  (0,2,1)
	//M14: i->j,k->j,k->i; find the third edge. check the first edge. O+(1,2,0)
	//M14: i->j,k->j,k->i; find the second edge. O- check the first edge. (2,1,0)
		int order0[3]={0,2,1};
		int order1[3]={1,2,0};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri01(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M14: i->j,k->j,k->i; find the second edge. I+ check the third edge.(0,1,2)
		//M14: i->j,k->j,k->i;  find the first edge. I- check the third edge.(1,0,2)
		//M14: i->j,k->j,k->i; find the first edge O-. check the second edge.(2,0,1)
		int order0[3]={0,1,2};
		int order1[3]={1,0,2};
		int order2[3]={2,0,1};

		//as first edge
		t=G->Cmotif_Tri11(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri10(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T0(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	//M23: i->j,j->k,i->k;  find the third edge.check the second edge;O+   (0,2,1)
	//M23: i->j,j->k,i->k; find the first edge. check the third edge. I- (1,0,2)
	//M23: i->j,j->k,i->k; find the first edge. O- check the second edge.

	//M23: i->j,j->k,i->k; find the second edge. O+ check the third edge.
	//M23: i->j,j->k,i->k;  find the third edge. I+ check the first edge.
	//M23: i->j,j->k,i->k; find the second edge I-. check the first edge.
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M23: i->j,j->k,i->k;  find the third edge.check the second edge;O+   (0,2,1)
	//M23: i->j,j->k,i->k; find the first edge. check the third edge. I- (1,0,2)
	//M23: i->j,j->k,i->k; find the first edge. O- check the second edge.(2,0,1)
		int order0[3]={0,2,1};
		int order1[3]={1,0,2};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//2)i->j,j->k,i->k; find the second edge. O+ check the third edge.(0,1,2)
		//2)i->j,j->k,i->k;  find the third edge. I+ check the first edge.(1,2,0)
		//2)i->j,j->k,i->k; find the second edge I-. check the first edge.(2,1,0)
		int order0[3]={0,1,2};
		int order1[3]={1,2,0};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri11(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T1(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M24: i->j,j->k,k->i;  find the third edge.check the second edge; I+  (0,2,1)
	//M24: i->j,j->k,k->i; find the first edge. check the third edge. I- (1,0,2)
	//M24: i->j,j->k,k->i; find the second edge. I- check the first edge.(2,1,0)
		int order0[3]={0,2,1};
		int order1[3]={1,0,2};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M24: i->j,j->k,k->i; find the second edge. O+ check the third edge.(0,1,2)
		//M24: i->j,j->k,k->i;  find the third edge. O+ check the first edge.(1,2,0)
		//M24: i->j,j->k,k->i; find the first edge O-. check the second edge.(2,0,1)
		int order0[3]={0,1,2};
		int order1[3]={1,2,0};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T10(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M35: i->j,k->i,j->k;  find the second edge. check the third edge; I+  (0,1,2)
	//M35: i->j,k->i,j->k; find the third edge. check the first edge. I+ (1,2,0)
	//M35: i->j,k->i,j->k; find the first edge. I- check the second edge.(2,0,1)
		int order0[3]={0,1,2};
		int order1[3]={1,2,0};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri11(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M35: i->j,k->i,j->k; find the third edge,O+ check the second edge.	(0,2,1)
		//M35: i->j,k->i,j->k;  find the first edge. O- check the third edge. (1,0,2)
		//M35: i->j,k->i,j->k; find the second edge O-. check the first edge.(2,1,0)
		int order0[3]={0,2,1};
		int order1[3]={1,0,2};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri00(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;

}
long Cmotif_T11(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M36: i->j,k->i,k->j;  find the second edge. check the third edge; I+  (0,1,2)
	//M36: i->j,k->i,k->j; find the third edge. check the first edge. O+
	//M36: i->j,k->i,k->j; find the second edge. O- check the first edge.
		int order0[3]={0,1,2};
		int order1[3]={1,2,0};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri01(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri00(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M36: i->j,k->i,k->j; find the third edge,I+ check the second edge.
		//M36: i->j,k->i,k->j;  find the first edge. O- check the third edge.
		//M36: i->j,k->i,k->j; find the first edge I-. check the second edge.
		int order0[3]={0,2,1};
		int order1[3]={1,0,2};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri00(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T8(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M45: i->j,i->k,j->k; find the second edge.check the third edge;O+ 	(0,1,2)
	//M45: i->j,i->k,j->k; find the first edge. check the third edge.O- (1,0,2)
	//M45: i->j,i->k,j->k; find the first edge. I- check the second edge.(2,0,1)
		int order0[3]={0,1,2};
		int order1[3]={1,0,2};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri00(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M45: i->j,i->k,j->k; find the third edge O+ check the second edge..
		//M45: i->j,i->k,j->k;  find the third edge. I+ check the first edge.
		//M45: i->j,i->k,j->k; find the second edge I-. check the first edge.
		int order0[3]={0,2,1};
		int order1[3]={1,2,0};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri11(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}

long Cmotif_T9(Graph *G,long id, int delta,float Q,int seed,Graph *M)
{
	long sum=0,t=0;
	long src=G->edges[id].src;
	long dest=G->edges[id].dest;
	if (G->nodeDegree(src)<=G->nodeDegree(dest))
	{
	//M46: i->j,i->k,k->j;  find the second edge, check the third edge;O+  (0,1,2)
	//M46: i->j,i->k,k->j; find the first edge. check the third edge.O- (1,0,2)
	//M46: i->j,i->k,k->j; find the second edge. I- check the first edge.
		int order0[3]={0,1,2};
		int order1[3]={1,0,2};
		int order2[3]={2,1,0};
		//as first edge
		t=G->Cmotif_Tri01(id, delta, src, dest, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri00(id, delta, src, dest, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, src, dest, order2, Q, seed, M);
		sum=sum+t;
	}
	else
	{
		//M46: i->j,i->k,k->j; find the third edge I+ check the second edge.
		//M46: i->j,i->k,k->j;  find the third edge. O+ check the first edge.
		//M46: i->j,i->k,k->j; find the first edge I-. check the second edge.
		int order0[3]={0,2,1};
		int order1[3]={1,2,0};
		int order2[3]={2,0,1};
		//as first edge
		t=G->Cmotif_Tri11(id, delta, dest, src, order0, Q, seed, M);
		sum=sum+t;
		//as second edge
		t=G->Cmotif_Tri01(id, delta, dest, src, order1, Q, seed, M);
		sum=sum+t;
		//as third edge
		t=G->Cmotif_Tri10(id, delta, dest, src, order2, Q, seed, M);
		sum=sum+t;	
	}
	return sum;
}


long Cmotif_0_2(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//000(M63) i->j;i->j;i->k; center is i;
	//010(M53) i->j;j->i;i->k; center=i
	//as first edge, find third edge O+, check the second edge.
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//000(M63) i->j;i->j;i->k; center is i;
	//as second edge, find third edge O+, check the first edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//000(M63) i->j;i->j;i->k; center is i;
	//as third edge, find first edge O-, check the second edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;

	return sum;
}

long Cmotif_1_3(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//001(M64) i->j;i->j;k->i;centter is i
	//as first edge,find third edge I+,check the second edge.
	//as second edge, find third edge I+, check the first edge.
	//as third edge, find first edge O-, check the second edge.

	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge O-, check the second edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;

	return sum;
}


long Cmotif_4_6(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//100(M55) i->j;j->i;j->k; center=j;
	//110(M65) i->j;i->j;j->k;center=j;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge, find the third edge O+, check the second edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find the third edge O+, check the first edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find the first edge I-, check the second edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;

	return sum;
}

long Cmotif_5_7(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//101(M56) i->j;j->i;k->j; center=j; 
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;

	return sum;
}

long Cmotif_8_9_16_17(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//O+,O-,O-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;
	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;
	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}

long Cmotif_10_11(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//I+,O-,I-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}

long Cmotif_12_13(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//O+,I-,O-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}

long Cmotif_14_15_22_23(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//I+,I-,I-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}

long Cmotif_18_19(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//I+,O-,O-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed11(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed00(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}

long Cmotif_20_21(Graph *G,long id, int delta, int center,int order[][3],float Q,int seed, Graph *M)
{
	//O+,I-,I-;
	long es,ed,et,sum=0,t=0;
	es=G->edges[id].src;
	ed=G->edges[id].dest;
	et=G->edges[id].tim;
	//as first edge,find third edge I+,check the second edge.
	t=G->Cmotif_Wed01(es,ed,et,delta,center,order[0],Q,seed,M);
	sum=sum+t;

	//as second edge, find third edge I+, check the first edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[1],Q,seed,M);
	sum=sum+t;

	//as third edge, find first edge I-, check the second edge.
	t=G->Cmotif_Wed10(es,ed,et,delta,center,order[2],Q,seed,M);
	sum=sum+t;
	return sum;
}