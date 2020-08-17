/*
 * C++ Program to Implement Adjacency List
 */
#include <iostream>
#include<fstream>
#include<cstdlib>
#include<cstring>
#include<string>
#include<vector>
#include<time.h>
#include<set>
#include <algorithm>
#include <map>
#include <random>
#include <sys/time.h>
#include<unordered_map>
#include <climits>
#include<algorithm>
using namespace std;


struct TEdge
{
public:
	long src,dest,tim,id;
	TEdge(long src,long dest,long tim,long id)
	{
		this->src=src;
		this->dest=dest;
		this->tim=tim;
		this->id=id;
	}
	const bool operator<(const TEdge& e) const{
		if(e.tim!=tim) return tim<e.tim;
		if(e.src!=src) return src<e.src;
		if(e.dest!=dest) return dest<e.dest;
		if(e.id!=id) return id<e.id;
		return false;
	}
};

class Etim
{
public:
	long node,tim,id; 
	Etim(long node,long tim,long id)
	{
		this->node=node;
		this->tim=tim;
		this->id=id;
	}
	const bool operator<(const Etim& e) const{
		return tim<e.tim;
	}
};

class Graph
{
//private:
public:	
	long Vn;
	long En;
	vector<vector<Etim>> Oadjlist;
	vector<vector<Etim>> Iadjlist;
	vector< unordered_map<long, vector<Etim> > > Hash;
	vector<TEdge> edges;

	
public:
	//Graph();
	inline void addEdge(long src,long dest,long tim,long id)
	{
		TEdge e(src, dest, tim, id);
		edges.push_back(e);
	}

	void Initialize()
	{
		unordered_map<long,long> vmap;
		int id=0;
		for(vector<TEdge>::iterator it=edges.begin();it!=edges.end();it++)
		{
			if(!vmap.count((*it).src)){
				vmap[(*it).src]=id++;
			}
			if(!vmap.count((*it).dest)){
				vmap[(*it).dest]=id++;
			}
		}
		Vn=id;
		En=edges.size();
		for(vector<TEdge>::iterator it=edges.begin();it!=edges.end();it++)
		{
			(*it).src=vmap[(*it).src];
			(*it).dest=vmap[(*it).dest];
		}
		vmap.clear();
		sort(edges.begin(),edges.end());

		for (int i = 0; i < (int) edges.size(); i++) {
			edges[i].id = i;
		}

		Oadjlist.clear();
		Iadjlist.clear();
		Hash.clear();

		Oadjlist.resize(id);
		Iadjlist.resize(id);
		Hash.resize(id);

		
		for(vector<TEdge>::iterator it=edges.begin();it!=edges.end();it++)
		{
			Etim _es((*it).src,(*it).tim,(*it).id);
			Etim _ed((*it).dest,(*it).tim,(*it).id);
			Etim _eh(0L,(*it).tim,(*it).id);
			Oadjlist[(*it).src].push_back(_ed);
			Iadjlist[(*it).dest].push_back(_es);
			Hash[(*it).src][(*it).dest].push_back(_eh);
		}
	}

	long nodeDegree(long node)
	{
		return Oadjlist[node].size()+Iadjlist[node].size();
	}

	long getVn()
	{
		return Vn;
	}

	long getEn()
	{
		return En;
	}

	long getHashn()
	{
		long sum=0;
		for (int i=0;i<Vn;i++)
		{
			for(int j=0;j<Vn;j++)
				if(Hash[i].find(j)!=Hash[i].end())
				{sum=sum+Hash[i][j].size();}
		}
		return sum;
	}

	void Printgraph()
	{
		for (int i=0;i<En;i++)
			cout<<edges[i].src<<" "<<edges[i].dest<<" "<<edges[i].tim<<" "<<edges[i].id<<endl;
	}
	int Cmotif_F7(long id, int delta, Graph*M);
	long Check_F(vector<Etim>*X,long left,long right);

};
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
	if (argc!=5)
	{cout<<"please input four parameter: infile, outfile, delta, edge sampling P!"<<endl; exit(1);}
	string infile=argv[1];
	string outfile=argv[2];
	int delta=atoi(argv[3]);
	float P=atof(argv[4]);
	
	string motiffile="../motifs/motif7.txt";
	cout<<"input file: "<<infile<<endl;
	cout<<"motif file: "<<motiffile<<endl;
	cout<<"output file: "<<outfile<<endl;
	cout<<"delta: "<<delta<<endl;
	cout<<"P: "<<P<<endl;
	

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
	char buf[100];
	const char* d=" ";
	long a[3];
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
		delete p;
		if (a[0] != a[1]) { 
			G.addEdge(a[0],a[1],a[2],0L);
		}
	}
	G.Initialize();
	cout<<"Finish reading "<<G.getEn()<<" edges and the last edge is "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
	in.close();


//
//====================================================read a motif====================================================
//
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
	in.close();

//
//====================================================sampling====================================================
//
	int seeds[10]={5,10,15,20,25,30,35,40,45,50};
	for (int ii=0;ii<10;ii++)
	{
		int seed_p=seeds[ii];
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
//====================================================counting====================================================
//
		double start,time1;
		long T=0, sumT=0;
		start=get_wall_time(); 		
		for (int i=0;i<ids.size();i++)
		{
			T=G.Cmotif_F7(ids[i], delta, &M);
			sumT=sumT+T;
		}
		time1=get_wall_time()-start;

		out.open(outfile.c_str(),ios::app);
		double sum=double(sumT)/(M.getEn()*P);
		cout<<"sumT: "<<sumT<<" sum:"<<sum<<endl;
		cout<<" time1:"<<time1<<endl;
		out<<sum<<" "<<time1<<endl;
		out.close();
	}
}

int Graph::Cmotif_F7(long id, int delta, Graph*M)
{
	//find the outdegree of dest.
	long es,ed,et,left0,right0,left,right,left1,right1,n0,n1,pos0,pos1,pos2,pos3,it,t=0,sum=0;
	double rn;
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;
	
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	TEdge* s3;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<Etim>* Z;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);

	int order[4][4]={{0,1,3,2},{1,2,0,3},{2,3,1,0},{3,0,2,1}};
	for(int i=0;i<4;i++)
	{
		long head0=0,head1=0;
		s0=&(M->edges[order[i][0]]);
		s1=&(M->edges[order[i][1]]);
		s2=&(M->edges[order[i][2]]);
		s3=&(M->edges[order[i][3]]);
		long ts[4]={-1,-1,-1,-1};
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;
		ts[order[i][0]]=et;
		X=&Oadjlist[ed];
		if(X->size())
		{
			//cout<<"size of X="<<X->size()<<endl;
			if(order[i][0]<order[i][1])
			{
				left0=ts[order[i][0]]+1;
				right0=ts[order[i][0]]+delta;
			}
			else
			{
				left0=ts[order[i][0]]-delta;
				right0=ts[order[i][0]]-1;
			}
			if(left0<=right0)
			{
				bool small=(X->size()<16);
				if(!small)
				{
					head0=lower_bound(X->begin(),X->end(),Etim(0,left0,0))-X->begin();				
				}	
				/*pos0=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
				pos1=upper_bound(X->begin(),X->end(),Etim(0,right,0))-X->begin()-1;*/
				//cout<<"pos0="<<pos0<<" pos1="<<pos1<<endl;
				//if(pos0<=X->size()&&pos1>=0&&pos0<=pos1)
				//{
					for (long j=head0;j<X->size();j++)
					{
						long _t=X->at(j).tim;
						if(_t<left0||_t>right0)
						{
							if(small) continue;
							else break;
						}
						n0=X->at(j).node;
						if(n0!=es)
						{
							mapMtoG[s1->dest]=n0;
							ts[order[i][1]]=_t;
							//find the third edge
							Y=&Iadjlist[es];
							if(Y->size())
							{
								//cout<<"size of Y="<<Y->size()<<endl;
								if(order[i][2]<order[i][0]&&order[i][2]<order[i][1]&&order[i][0]<order[i][1])
								{
									left=ts[order[i][1]]-delta;
									right=ts[order[i][0]]-1;
								}
								else if (order[i][2]<order[i][0]&&order[i][2]<order[i][1]&&order[i][0]>order[i][1])
								{
									left=ts[order[i][0]]-delta;
									right=ts[order[i][1]]-1;
								}
								else if (order[i][2]>order[i][0]&&order[i][2]>order[i][1]&&order[i][0]<order[i][1])
								{
									left=ts[order[i][1]]+1;
									right=ts[order[i][0]]+delta;
								}
								else if (order[i][2]>order[i][0]&&order[i][2]>order[i][1]&&order[i][0]>order[i][1])
								{
									left=ts[order[i][0]]+1;
									right=ts[order[i][1]]+delta;
								}
								else if (order[i][2]<order[i][0]&&order[i][2]>order[i][1])
								{
									left=ts[order[i][1]]+1;
									right=ts[order[i][0]]-1;
								}
								else{ std::cout<<"the order is wrong! first three orders are"<<order[i][0]<<order[i][0]<<order[i][0]; exit(1); }
								if(left<=right)
								{
									bool small=(Y->size()<16);
									if(!small)
									{
										head1=lower_bound(Y->begin(),Y->end(),Etim(0,left,0))-Y->begin();				
									}	
									//pos2=lower_bound(Y->begin(),Y->end(),Etim(0,left,0))-Y->begin();
									//pos3=upper_bound(Y->begin(),Y->end(),Etim(0,right,0))-Y->begin()-1;	
									//cout<<"pos2="<<pos2<<" pos3="<<pos3<<endl;
									//if (pos2<=Y->size()&&pos3>=0&&pos2<=pos3)
									{
										for(long k=head1;k<Y->size();k++)
										{
											long _tt=Y->at(k).tim;
											if(_tt<left||_tt>right)
											{
												if(small) continue;
												else break;
											}
											n1=Y->at(k).node;
											if(n1!=ed&&n1!=n0)
											{
												mapMtoG[s2->src]=n1;
												ts[order[i][2]]=_tt;
												//check
												Z=&Hash[mapMtoG[s3->src]][mapMtoG[s3->dest]];
												if(Z->size())
												{
													switch(order[i][3]){
													case 0:
														{
															left1=ts[3]-delta;
															right1=ts[1]-1;
															break;
														}
													case 1:
														{
															left1=ts[0]+1;
															right1=ts[2]-1;
															break;
														}
													case 2:
														{
															left1=ts[1]+1;
															right1=ts[3]-1;
															break;
														}
													case 3:
														{
															left1=ts[2]+1;
															right1=ts[0]+delta;
															break;
														}
													default:
														{
															std::cout<<"wrond order"<<order[i][3];
															exit(1);
														}
													}
													if(left<=right)
													{
														t=Check_F(Z,left1,right1);
														sum=sum+t;
													}
												}
											}//endif
										}//endfor
									}//endif
								}//endif
							}//endif
						}//endif
					}//endfor
				//}//endif
			}//end if
		}//end if
	}//end for 
	return sum;
}

long Graph::Check_F(vector<Etim>* X, long left,long right)
{
	long pos0,pos1,tr=0;
	pos0=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
	pos1=upper_bound(X->begin(),X->end(),Etim(0,right,0))-X->begin()-1;
	if (pos0>X->size()||pos1<0||(pos0>pos1))
	{
		tr=0;
	}
	else 
	{
		tr=pos1-pos0+1;
	}		
	return tr;
}
