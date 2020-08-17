#ifndef graph_es_H
#define graph_es_H
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
	long src,dst,tim,id;
	TEdge(long src,long dst,long tim,long id)
	{
		this->src=src;
		this->dst=dst;
		this->tim=tim;
		this->id=id;
	}
	const bool operator<(const TEdge& e) const{
		if(e.tim!=tim) return tim<e.tim;
		if(e.src!=src) return src<e.src;
		if(e.dst!=dst) return dst<e.dst;
		if(e.id!=id) return id<e.id;
		return false;
	}
};



class Graph
{
//private:
public:	
	long Vn;
	long En;

	vector<TEdge> edges_;
	vector< vector< TEdge > > adj_list_;
	vector< vector< TEdge > > revadj_list_;
	vector< unordered_map<int, vector<TEdge> > > adjMap_;

	//Structure for algo
	vector<int> edgeCount_;
	vector<int> mapGM_;
	vector<int> mapMG_;
	//Motifs
	int Vm_;
	vector< pair<int, int> > edgesM_;
	
public:
	//Graph();
	void addEdge(long src,long dst,long tim,long id)
	{
		TEdge e(src, dst, tim, id);
		edges_.push_back(e);
	}

	void Initialize()
	{
		unordered_map<long,long> vmap;
		int id=0;
		for(vector<TEdge>::iterator it=edges_.begin();it!=edges_.end();it++)
		{
			if(!vmap.count((*it).src)){
				vmap[(*it).src]=id++;
			}
			if(!vmap.count((*it).dst)){
				vmap[(*it).dst]=id++;
			}
		}
		Vn=id;
		En=edges_.size();
		for(vector<TEdge>::iterator it=edges_.begin();it!=edges_.end();it++)
		{
			(*it).src=vmap[(*it).src];
			(*it).dst=vmap[(*it).dst];
		}
		vmap.clear();
		sort(edges_.begin(),edges_.end());

		for (int i = 0; i < (int) edges_.size(); i++) {
			edges_[i].id = i;
		}

		adj_list_.clear();
		revadj_list_.clear();
		adjMap_.clear();

		adj_list_.resize(id);
		revadj_list_.resize(id);
		adjMap_.resize(id);

		
		for(vector<TEdge>::iterator it=edges_.begin();it!=edges_.end();it++)
		{
			
			adj_list_[(*it).src].push_back(*it);
			revadj_list_[(*it).dst].push_back(*it);
			adjMap_[(*it).src][(*it).dst].push_back(*it);
		}
		

	}
	void cle_motif()
	{
		edgesM_.clear();
		Vm_=-1;
	}
	void motif(int src, int dst)
	{
		edgesM_.push_back(make_pair(src,dst));
	}

	void setVm_()
	{
		int maxnode=0;
		for (int i=0;i<edgesM_.size();i++)
		{
			maxnode=max(maxnode,edgesM_[i].first+1);
			maxnode=max(maxnode,edgesM_[i].second+1);
		}
		Vm_=maxnode;
	}

	long nodeDegree(long node)
	{
		return adj_list_[node].size()+revadj_list_[node].size();
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
				if(adjMap_[i].find(j)!=adjMap_[i].end())
				{sum=sum+adjMap_[i][j].size();}
		}
		return sum;
	}
	void Printgraph()
	{
		for (int i=0;i<En;i++)
		{cout<<edges_[i].src<<" "<<edges_[i].dst<<" "<<edges_[i].tim<<" "<<edges_[i].id<<endl;}
	}

	double ExactCountMotifs(int delta, vector<long>&ids) {
		//cout<<"test ExactCountMotifs"<<endl;
		vector<TEdge> fs;
		for(int i=0;i<ids.size();i++)
		{
			fs.push_back(edges_[ids[i]]);
		}
		
		if (edgesM_.size()==3)
		{
			int S_order[4][4]={{0,2,1,-1},{1,0,2,-1},{2,0,1,-1},{-1,-1,-1,-1}};
			return Cmotifs(delta,-1,edges_,fs, S_order);
		}
		if (edgesM_.size()==4)
		{
			int S_order[4][4]={{0,3,1,2},{1,0,3,2},{2,3,0,1},{3,0,1,2}};
			return Cmotifs(delta,-1,edges_,fs, S_order);
		}
    }

	double Cmotifs(int delta, int p,vector<TEdge>&edges,vector<TEdge>&fs,int S_order[][4])
	{
		
		double sum=0.0;
		int l=edgesM_.size();
		//cout<<"l="<<l<<" "<<edges.size()<<" "<<fs[0].id<<endl;
		for (int i=0;i<l;i++)
		{
			
			double res = 0.0;
			int flag=0;
			//cout<<"i="<<i<<endl;
			edgeCount_.clear();
    		edgeCount_.resize(En, 0);
    		mapGM_.clear();
    		mapGM_.resize(En, -1);
    		mapMG_.clear();
    		mapMG_.resize(Vm_, -1);


			//vector<int> order;
			vector<int> ts;
			vector<int>MM;
			MM.resize(l,-1);
    		vector<int> eStack;
			//vector<int> t_max;
			//t_max.resize(l,-1);
			

    		int j=0, eG = 0, eM = S_order[i][j], uG = -1, vG = -1, uM = -1, vM = -1;
    		const int INF = 2147483647;//numeric_limits<int>::max();
    		int _t = INF;
			//t_max[eM]=_t;
    		while (flag==0){
				//cout<<"eM="<<eM<<" before eG="<<eG<<" ";
    			int last = eStack.empty() ? -1 : eStack.back();
    			eG = FindNextMatch(eM, eG, mapMG_, mapGM_, _t, edges, fs, last,MM,i);
				//cout<<" eG="<<eG<<" edges[eG].tim="<<edges[eG].tim<<" ";
				//cout<<"size of order="<<order.size()<<" i="<<i<<" eM="<<eM<<" eG="<<eG<<" time="<<edges[eG].tim<<" _t="<<_t<<endl;
    			TEdge edge = {0, 0, 0, INF};
    			if (eG < (int) edges.size()) {
    				if (eStack.size()== int(edgesM_.size() - 1)) {
    					// Apply the weight function
    					if (p == -1) res += 1;
    					else res += 1.0/p;
						//cout<<"\n res="<<res<<" "<<MM[0]<<" "<<MM[1]<<" "<<MM[2]<<" "<<edges[eStack[0]].tim<<" "<<edges[eStack[1]].tim<<" "<<edges[eG].tim<<endl;
						eG+=1;//last edge is not 1.
    				} else {
    					edge = edges[eG];
    					uG = edge.src, vG = edge.dst; 
    					uM = edgesM_[eM].first, vM = edgesM_[eM].second;

    					mapGM_[uG] = uM;
    					mapGM_[vG] = vM;
    					mapMG_[uM] = uG;
    					mapMG_[vM] = vG;

    					edgeCount_[uG] += 1;
    					edgeCount_[vG] += 1;
						ts.push_back(_t);
						/*if(eM==0&&t_max[2]==-1)
						{
							t_max[2] = edge.tim + delta;
						}
						else if(eM==1&&t_max[0]==-1)
						{
							int tem=edge.tim-delta;
							t_max[0]=max(tem,0);
						}
						else if (eM==2)
						{
							if(t_max[0]==-1)
							{
								int tem=edge.tim-delta;
								t_max[0]=max(tem,0);
							}
							if(t_max[1]==-1)
							{
								t_max[1]=edge.tim-1;
							}
						}*/
    					eStack.push_back(eG);
						MM[eM]=eG;
						j++;
						eM=S_order[i][j];
						
						if(eM==0)
						{
							if(l==4&&j==2)
							{ eG=MM[2]-1;}
							else
							{ eG=eG-1;}
							if(l==4&&MM[3]!=-1)
							{
								_t=edges[MM[3]].tim-delta;
							}
							else if(MM[2]!=-1)
							{
								_t=edges[MM[2]].tim-delta;
							}
							else if(MM[1]!=-1)
							{
								_t=edges[MM[1]].tim-delta;
							}
						}
						else  if (eM==1)
						{
							if(MM[0]!=-1)
							{
								eG=MM[0]+1;
							}
							if(MM[2]!=-1)
							{
								_t=edges[MM[2]].tim-1;
							}
							else if(l==4&&MM[3]!=-1)
							{
								_t=edges[MM[3]].tim-1;
							}
							
						}
						else if(eM==2)
						{
							if(MM[1]!=-1)
							{
								eG=MM[1]+1;
							}
							else if(MM[0]!=-1)
							{
								eG=MM[0]+1;
							}
							if(l==3)
							{
								_t=edges[MM[0]].tim+delta;
							}
							else
							{
								_t=edges[MM[3]].tim-1;
							}
						}
						else if(eM==3)
						{
							if(j==2)
							{
								eG=MM[1]+1;
							}
							else
							{
								eG=eG+1;
							}
							if(MM[0]!=-1)
							{
								_t=edges[MM[0]].tim+delta;
							}
							else
							{
								_t=edges[MM[2]].tim+delta;
							}
						}
						//cout<<"next  eM="<<eM<<" _t="<<_t<<" eG="<<eG<<endl;
						/*if(t_max[eM]!=-1)
						{
							_t=t_max[eM];
						}
						else
						{
							cout<<"t_max["<<eM<<"]"=t_max[eM]<<endl;
							exit(1);
						}*/
    				}
    			}
				//cout<<"eM="<<eM<<" edge: "<<edge.id<<" "<<edge.dst<<" "<<edge.src<<" "<<edge.tim<<" _t"<<_t<<endl;
								
				while (eG >= edges.size() || ((eM!=0)&&(edges[eG].tim > _t))||(edge.id!=INF&&eM==0&&j!=0&&edge.tim<_t)) {
					//cout<<"_t"<<_t<<" ";
    				if (!eStack.empty()) {
						j--;
						if(S_order[i][j]==0&&j!=0)
						{
							eG = eStack.back() - 1;
							edge = edges[eG + 1];
						}
						else
						{
							eG = eStack.back() + 1;
							edge = edges[eG - 1];
						}
    				_t=ts.back();
					ts.pop_back();
    				eStack.pop_back();
					MM[S_order[i][j]]=-1;

    				uG = edge.src, vG = edge.dst;
    				uM = edgesM_[eM].first, vM = edgesM_[eM].second;

    				if (eStack.empty()) {
    					_t = INF;
    				}
    				edgeCount_[uG] -= 1;
    				edgeCount_[vG] -= 1;

    				if (edgeCount_[uG] == 0) {
    					uM = mapGM_[uG];
    					mapMG_[uM] = -1;
    					mapGM_[uG] = -1;
    				}

    				if (edgeCount_[vG] == 0) {
    					vM = mapGM_[vG];
    					mapMG_[vM] = -1;
    					mapGM_[vG] = -1;
    				}

    				eM =S_order[i][j];
						//cout<<"eM="<<eM<<" size of order="<<order.size()<<" eG="<<eG<<" time="<<edge.tim<<" max="<<_t<<endl;
    				} else {//eM=i
						//cout<<res<<endl;
    					sum=sum+res;
						//cout<<"sum="<<sum<<" eM="<<eM<<" eG="<<eG<<endl;
						flag=1;
						break;
    				}
    			}
			}
			//cout<<"\n"<<i<<"round"<<" res="<<res<<" sum="<<sum<<endl;
		}
		return sum;
	}
	
	inline int FindNextMatch(
    	int eM, 
    	int eG, 
    	vector<int>& mapMG, 
    	vector<int>& mapGM, 
    	int _t,
    	vector<TEdge>& edges,
		vector<TEdge>& fs,
    	int last,
		vector<int>&MM,
		int round) {

    	int uM, vM, uG, vG;
    	uM = edgesM_[eM].first;
    	vM = edgesM_[eM].second;
		
    	uG = mapMG_[uM];
    	vG = mapMG_[vM];
		
    	vector<TEdge>* S;
    	int head = -1;
    	if (uG >= 0 && vG >= 0) {
    		S = &adjMap_[uG][vG];
    	} else if (uG >= 0) {
    		S = &adj_list_[uG];
    	} else if (vG >= 0) {
    		S = &revadj_list_[vG];
    	} else {
    		S = &edges;
    	}
		if(last==-1)
		{
			S=&fs;
		}
		/*cout<<"\nthe size of S="<<S->size()<<endl;
		for (int k=0;k<S->size();k++)
		{
			cout<<(*S)[k].tim<<" ";
		}
		cout<<endl;*/
    	bool small = (S->size() < 16);
    	if (!small) {
			if(eM==0&&last!=-1)
	    	{
				head = upper_bound(S->begin(), S->end(), edges[eG]) - S->begin()-1;
				//返回数组中第一个大于该元素的下标 -1即返回最后一个小于等于该元素的下标。
			}
			else
			{
				head = lower_bound(S->begin(), S->end(), edges[eG]) - S->begin();
			}//返回数组中第一个大于等于该元素的下标。}
	    }
		
		if(eM==0&&last!=-1)
		{
			if(head==-1)
			{
				head=S->size()-1;
			}
			//cout<<"S"<<" ";
			for (int i = head; i >=0; i--) {
    			auto edge = (*S)[i];
				//cout<<edge.tim<<" ";
    			if (edge.id > eG || edge.tim < _t) {
    				if (small) continue;
    				else break;
    			}
				int tem;
				if(round==2){tem=MM[2];}
				else{tem=last;}
    			if (last != -1 && edge.tim>=edges[tem].tim) {
    				continue;
    			}

    			int _eG = edge.id;
    			int _uG = edge.src, _vG = edge.dst;
    			if (uG == _uG || (uG < 0 && mapGM_[_uG] < 0)) {
    				if (vG == _vG || (vG < 0 && mapGM_[_vG] < 0)) {
    					return _eG;
    				}
    			} 
    		}
		}
		else
		{
			if (head==-1)
			{
				head=0;
			}
    		for (int i = head; i < (int) S->size(); i++) {
    			auto edge = (*S)[i];
    			if (edge.id < eG || edge.tim > _t) {
    				if (small) continue;
    				else break;
    			}
				int tem;
				if(eM-1>=0&&MM[eM-1]!=-1)
				{
					tem=MM[eM-1];
				}
				else if(eM-2>=0&&MM[eM-2]!=-1)
				{
					tem=MM[eM-2];
				}
				else if(eM-3>=0&&MM[eM-3]!=-1)
				{
					tem=MM[eM-3];
				}
    			if (last != -1 && edge.tim <= edges[tem].tim) {
    				continue;
    			}
				
    			int _eG = edge.id;
				
    			int _uG = edge.src, _vG = edge.dst;
    			if (uG == _uG || (uG < 0 && mapGM_[_uG] < 0)) {
    				if (vG == _vG || (vG < 0 && mapGM_[_vG] < 0)) {
    					return _eG;
    				}
    			} 
    		}
		}
		
    	return edges.size();
    }


};

#endif