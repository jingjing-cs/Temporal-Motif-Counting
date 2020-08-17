#ifndef graph_ews_H
#define graph_ews_H
#include<iostream>
#include<vector>
#include<unordered_map>
#include<string>
#include <climits>
#include<random>
#include<algorithm>
using namespace std;

//#define N 2147483647
//#define bit 32


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

	int cal_Wid( int center, int order[][3]);
	int cal_Tid();
	int cal_Rid();
	long Cmotif_Pair(long id, int delta, float Q,int seed, Graph* M);
	long Cmotif_Tri00(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M);
	long Cmotif_Tri01(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M);
	long Cmotif_Tri10(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M);
	long Cmotif_Tri11(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M);

	long Cmotif_Wed00(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M);
	long Cmotif_Wed01(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M);
	long Cmotif_Wed10(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M);
	long Cmotif_Wed11(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M);

	int Cmotif_F7(long id, int delta, Graph*M);
	int Cmotif_F8(long id, int delta, Graph*M);

	long Check(vector<Etim>* X, long tim0,long tim1,int order, int delta);
	long Check_F(vector<Etim>*X,long left,long right);
	double Rnum(int seed);
};


int Graph::cal_Wid(int center, int order[][3])
{
	int ids[4],e[3];
	int bit_=16;
	for (int i=0;i<3;i++)
	{
		if(edges[i].src<edges[i].dest)
			{e[i]=(edges[i].src<<bit_|edges[i].dest);}
		else {e[i]=(edges[i].dest<<bit_|edges[i].src);}
	}
	//****e0=e1(in undirected graph)****;
	//8 000(M63) i->j;i->j;i->k;
	//9 001(M64) i->j;i->j;k->i;
	//10 010(M53) i->j;j->i;i->k;
	//11 011(M54) i->j;j->i;k->i;
	//12 100(M55) i->j;j->i;j->k;
	//13 101(M56) i->j;j->i;k->j;
	//14 110(M65) i->j;i->j;j->k;
	//15 111(M66) i->j;i->j;k->j;

	//****e0=e2(in undirected graph)****;
	//16 000(M41) i->j;i->k;i->j;
	//17 001(M42) i->j;i->k;j->i;
	//18 010(M31) i->j;k->i;i->j;
	//19 011(M32) i->j;k->i;j->i;
	//20 100(M22) i->j;j->k;j->i;
	//21 101(M21) i->j;j->k;i->j;
	//22 110(M12) i->j;k->j;j->i;
	//23 111(M11) i->j;k->j;i->j;
	
	//****e1=e2(in undirected graph)****;
	//24 000(M43) i->j;i->k;i->k;
	//25 001(M44) i->j;i->k;k->i;
	//26 010(M33) i->j;k->i;i->k;
	//27 011(M34) i->j;k->i;k->i;
	//28 100(M25) i->j;j->k;j->k;
	//29 101(M26) i->j;j->k;k->j;
	//30 110(M15) i->j;k->j;j->k;
	//31 111(M16) i->j;k->j;k->j;
	if(e[0]==e[1])
	{
		ids[3]=0;
		order[0][0]=0;
		order[0][1]=2;
		order[0][2]=1;
		order[1][0]=1;
		order[1][1]=2;
		order[1][2]=0;
		order[2][0]=2;
		order[2][1]=0;
		order[2][2]=1;

	}
	else if (e[0]==e[2])
	{
		ids[3]=8;
		order[0][0]=0;
		order[0][1]=1;
		order[0][2]=2;
		order[1][0]=1;
		order[1][1]=0;
		order[1][2]=2;
		order[2][0]=2;
		order[2][1]=1;
		order[2][2]=0;
	}
	else if (e[1]==e[2])
	{
		ids[3]=16;
		order[0][0]=0;
		order[0][1]=1;
		order[0][2]=2;
		order[1][0]=1;
		order[1][1]=0;
		order[1][2]=2;
		order[2][0]=2;
		order[2][1]=0;
		order[2][2]=1;
	}
	for (int i=0;i<3;i++)
	{
		if (edges[i].src==center)
		{
			ids[i]=0;
		}
		else if(edges[i].dest==center)
		{
			ids[i]=1;
		}
		else 
		{
			cout<<"wrong wedge motif."<<endl;
			exit(1);
		}
	}
	int id=(ids[0]<<2)|(ids[1]<<1)|(ids[2]);
	id=id+ids[3];
	return id;
}


int Graph::cal_Tid()
{	
	int ids[4];
	if (edges[1].src==edges[0].dest||edges[1].dest==edges[0].dest)
	{
		ids[0]=0;
	}
	else
	{
		ids[0]=1;
	}

	for (int i=0;i<Vn;i++)
	{
		if(edges[i].src<edges[i].dest)
		{	ids[i+1]=0; }
		else
		{	ids[i+1]=1;}
	}
	int ID=(ids[0]<<3)|(ids[1]<<2)|(ids[2]<<1)|(ids[3]);
	return ID;
}

long Graph::Cmotif_Pair(long id, int delta, float Q,int seed, Graph *M)
{
	int order[3][3]={{0,1,2},{1,0,2},{2,0,1}};
	long es,ed,et,left,right,it,t=0,T=0;
	double rn;
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;

	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);
	for (int i=0;i<3;i++)
	{
		s0=&(M->edges[order[i][0]]);
		s1=&(M->edges[order[i][1]]);
		s2=&(M->edges[order[i][2]]);
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;
		long head=0;
		
		X=&Hash[mapMtoG[s1->src]][mapMtoG[s1->dest]];
		bool small=(X->size()<16);
		if(X->size())
		{
			if(order[i][0]<order[i][1])
			{
				left=et+1;
				right=et+delta;
			}
			else 
			{
				left=et-delta;
				right=et-1;
			}
			if(left<=right)
			{
				if(!small)
				{
					head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
				}
					
				for(it=head;it<X->size();it++)
				{
						
					long _t=X->at(it).tim;
					if(_t<left||_t>right)
					{
						if(small) continue;
						else break;
					}
					rn=Rnum(seed);
					if(rn<=Q)
					{
							
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(et<_t)
							{
								t=Check(Y,et,_t,order[i][2],delta);
							}
							else
							{
								t=Check(Y,_t,et,order[i][2],delta);
							}
							//cout<<"t"<<t<<endl;
							T=T+t;
						}
						//std::cout<<"order[i][2]"<<order[i][2]<<" "<<t<<" "<<T<<"\n";
					}
				}
					
			}
		}
		
	}	
	return T;
}


long Graph::Check(vector<Etim>* X, long tim0,long tim1,int order, int delta)
{
	long left,right,tr=0;
	int pos0,pos1;
	switch(order)
	{
	case 0:
		{
			left=tim1-delta;
			right=tim0-1;
			break;
		}
	case 1:
		{
			left=tim0+1;
			right=tim1-1;
			break;
		}
	case 2:
		{
			left=tim1+1;
			right=tim0+delta;
			break;
		}
	default:
		{
			cout<<"invalid order: "<<order<<"\n";
			exit(1);
		}
	}
	if(left<=right)
	{
		pos0=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
		pos1=upper_bound(X->begin(),X->end(),Etim(0,right,0))-X->begin()-1;
		//cout<<"check pos0 pos1"<<pos0<<" "<<pos1<<endl;
		if (pos0>X->size()||pos1<0||(pos0>pos1))
			{
				tr=0;
			}
			else 
			{
				tr=pos1-pos0+1;
			}	
	}
	return tr;
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


long Graph::Cmotif_Tri00(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M)
{
	long no,es,ed,et,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;
	X=&Oadjlist[src];
	if(X->size())
	{
		mapMtoG.resize(M->getVn(),-1);
		s0=&(M->edges[order[0]]);
		s1=&(M->edges[order[1]]);
		s2=&(M->edges[order[2]]);
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;
		left=et-delta;
		right=et-1;
		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
				
			}		
			for(it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].dest;
				if(no!=dest)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->dest]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
							}
							//cout<<"t"<<t<<endl;
							sum=sum+t;
						}
					}
				}
			}
			
		}
	}
	return sum;
}
long Graph::Cmotif_Tri01(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M)
{
	long no,es,ed,et,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;
	X=&Oadjlist[src];
	if(X->size())
	{
		mapMtoG.resize(M->getVn(),-1);
		s0=&(M->edges[order[0]]);
		s1=&(M->edges[order[1]]);
		s2=&(M->edges[order[2]]);
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;

		left=et+1;;
		right=et+delta;

		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
				
			}
		
			for(it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
			
				no=X->at(it).node;//edges[X->at(it).id].dest;
				if(no!=dest)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->dest]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
							}
							//cout<<"t"<<t<<endl;
							sum=sum+t;
						}
					}
				}
			}
		}
	}
	return sum;
}
long Graph::Cmotif_Tri10(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M)
{
	long no,es,ed,et,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;
	X=&Iadjlist[src];
	if(X->size())
	{
		mapMtoG.resize(M->getVn(),-1);
		s0=&(M->edges[order[0]]);
		s1=&(M->edges[order[1]]);
		s2=&(M->edges[order[2]]);
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;

		left=et-delta;
		right=et-1;

		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
				
			}
			
			for(it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].src;
				if(no!=dest)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->src]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
							}
							//cout<<"t"<<t<<endl;
							sum=sum+t;
						}
					}
				}
			}
		}
	}
	return sum;
}
long Graph::Cmotif_Tri11(long id, int delta, long src,long dest, int order[],float Q,int seed,Graph *M)
{
	long no,es,ed,et,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	
	es=edges[id].src;
	ed=edges[id].dest;
	et=edges[id].tim;
	X=&Iadjlist[src];
	if(X->size())
	{
		mapMtoG.resize(M->getVn(),-1);
		s0=&(M->edges[order[0]]);
		s1=&(M->edges[order[1]]);
		s2=&(M->edges[order[2]]);
		mapMtoG[s0->src]=es;
		mapMtoG[s0->dest]=ed;

		left=et+1;
		right=et+delta;

		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
			}
			for(it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].src;
				if(no!=dest)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->src]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
							}
							//cout<<"t"<<t<<endl;
							sum=sum+t;
						}
					}
				}
			}
		}
	}
	return sum;
}



//for star motif: outdegree-vector<TEdge> edges
long Graph::Cmotif_Wed00(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M)
{
	long no,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);

	s0=&(M->edges[order[0]]);
	s1=&(M->edges[order[1]]);
	s2=&(M->edges[order[2]]);

	mapMtoG[s0->src]=es;
	mapMtoG[s0->dest]=ed;

	X=&Oadjlist[mapMtoG[center]];
	if(X->size())
	{
		left=et-delta;
		right=et-1;
		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();				
			}
			for (it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].dest;
				if(no!=es&&no!=ed)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->dest]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
								sum=sum+t;
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
								sum=sum+t;
							}
						}
					}
				}		
			}
			
		}
	}
	return sum;
}

//for star motif outdegree+
long Graph::Cmotif_Wed01(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M)
{
	long no,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);

	s0=&(M->edges[order[0]]);
	s1=&(M->edges[order[1]]);
	s2=&(M->edges[order[2]]);

	mapMtoG[s0->src]=es;
	mapMtoG[s0->dest]=ed;

	X=&Oadjlist[mapMtoG[center]];
	if(X->size())
	{
		left=et+1;
		right=et+delta;
		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
			}						
			for (it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].dest;
				if(no!=es&&no!=ed)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->dest]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
								sum=sum+t;
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
								sum=sum+t;
							}
						}
					}
				}		
			}
		}
	}
	return sum;
}
//for star motif: indegree-
long Graph::Cmotif_Wed10(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M)
{
	long no,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);

	s0=&(M->edges[order[0]]);
	s1=&(M->edges[order[1]]);
	s2=&(M->edges[order[2]]);

	mapMtoG[s0->src]=es;
	mapMtoG[s0->dest]=ed;

	X=&Iadjlist[mapMtoG[center]];
	if(X->size())
	{
		left=et-delta;
		right=et-1;
		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();
			}			
			for (it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].src;
				if(no!=es&&no!=ed)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->src]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
								sum=sum+t;
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
								sum=sum+t;
							}
						}
					}
				}		
			}
		}
	}
	return sum;
}

//for star motif: indegree+
long Graph::Cmotif_Wed11(long es,long ed, long et, int delta, int center, int order[],float Q,int seed,Graph *M)
{
	long no,left,right,it,t=0,sum=0,head=0;
	double rn;
	TEdge* s0;
	TEdge* s1;
	TEdge* s2;
	vector<Etim>* X;
	vector<Etim>* Y;
	vector<long> mapMtoG;
	mapMtoG.resize(M->getVn(),-1);

	s0=&(M->edges[order[0]]);
	s1=&(M->edges[order[1]]);
	s2=&(M->edges[order[2]]);

	mapMtoG[s0->src]=es;
	mapMtoG[s0->dest]=ed;

	X=&Iadjlist[mapMtoG[center]];
	if(X->size())
	{
		left=et+1;
		right=et+delta;
		if(left<=right)
		{
			bool small=(X->size()<16);
			if(!small)
			{
				head=lower_bound(X->begin(),X->end(),Etim(0,left,0))-X->begin();				
			}						
			for (it=head;it<X->size();it++)
			{
				long _t=X->at(it).tim;
				if(_t<left||_t>right)
				{
					if(small) continue;
					else break;
				}
				no=X->at(it).node;//edges[X->at(it).id].src;
				if(no!=es&&no!=ed)
				{
					rn=Rnum(seed);
					if(rn<=Q)
					{
						mapMtoG[s1->src]=no;
						Y=&Hash[mapMtoG[s2->src]][mapMtoG[s2->dest]];
						if(Y->size())
						{
							if(order[0]<order[1])
							{
								t=Check(Y,et,_t,order[2],delta);
								sum=sum+t;
							}
							else
							{
								t=Check(Y,_t,et,order[2],delta);
								sum=sum+t;
							}
						}
					}
				}		
			}		
		}
	}
	return sum;
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



int Graph::Cmotif_F8(long id, int delta, Graph*M)
{
	//find the outdegree of dest.
	long es,ed,et,left,right,left0,right0,left1,right1,n0,n1,pos0,pos1,pos2,pos3,it,t=0,sum=0;
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

	int order[4][4]={{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
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
		X=&Oadjlist[es];
		if(X->size())
		{
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
				pos1=upper_bound(X->begin(),X->end(),Etim(0,right,0))-X->begin()-1;
				if(pos0<=X->size()&&pos1>=0&&pos0<=pos1)*/
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
						if(n0!=ed)
						{
							mapMtoG[s1->dest]=n0;
							ts[order[i][1]]=_t;
							//find the third edge
							Y=&Iadjlist[ed];
							if(Y->size())
							{
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
									/*pos2=lower_bound(Y->begin(),Y->end(),Etim(0,left,0))-Y->begin();
									pos3=upper_bound(Y->begin(),Y->end(),Etim(0,right,0))-Y->begin()-1;	
									if (pos2<=Y->size()&&pos3>=0&&pos2<=pos3)*/
									//{
										for(long k=head1;k<Y->size();k++)
										{
											long _tt=Y->at(k).tim;
											if(_tt<left||_tt>right)
											{
												if(small) continue;
												else break;
											}
											n1=Y->at(k).node;
											if(n1!=es&&n1!=n0)
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
														//if(t>0)
														/*{
															std::cout<<"results="<<t<<"nodes: "<<src<<" "<<dest<<" "<<n0<<" "<<n1<<"\n";
															std::cout<<"edges: "<<ws[0]<<" "<<ws[1]<<" "<<ws[2]<<ws[3]<<"\n";
															std::cout<<"order: "<<order[i][0]<<" "<<order[i][1]<<" "<<order[i][2]<<" "<<order[i][3]<<"\n";
														}*/
														sum=sum+t;
													}
												}
											}
										}
									//}
								}
							}//end if
						}//end if
					}//end if
				//}//end for
			}//end if
		}//end if
	}//end for 
	return sum;
}

double Graph::Rnum(int seed)
{
	static std::default_random_engine e(seed);
	static std::uniform_real_distribution<double> u(0.0,1.0);
	return u(e);
}



#endif