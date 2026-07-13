#include <cstdio>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
using namespace std;
typedef uint64_t LK;
static const int AOFF=200,BOFF=200,COFF=20000,ASPAN=400,BSPAN=400,CSPAN=40000;
struct Pt{int x,y;};
Pt c4(int x,int y,int r,int n){ if(r==0)return{x,y}; if(r==1)return{n-1-y,x}; if(r==2)return{n-1-x,n-1-y}; return{y,n-1-x}; }
static inline int igcd(int a,int b){a=abs(a);b=abs(b);while(b){int t=a%b;a=b;b=t;}return a;}
static inline LK line_of(Pt p,Pt q){
    int a=-(q.y-p.y),b=(q.x-p.x),c=(q.y-p.y)*p.x-(q.x-p.x)*p.y;
    int g=igcd(abs(a),abs(b));g=igcd(g,abs(c));
    if(g!=0){a/=g;b/=g;c/=g;}
    if(a<0||(a==0&&b<0)||(a==0&&b==0&&c<0)){a=-a;b=-b;c=-c;}
    LK K=(LK)(a+AOFF); K=K*BSPAN+(b+BOFF); K=K*CSPAN+(c+COFF);
    return K;
}
int main(){
    int cells[10][2]={{0,7},{1,5},{3,6},{4,2},{4,9},{5,3},{6,1},{8,2},{8,7},{9,0}};
    int m=10,n=2*m;
    vector<Pt> pts; for(int i=0;i<m;i++)for(int r=0;r<4;r++)pts.push_back(c4(cells[i][0],cells[i][1],r,n));
    printf("npts=%d\n",(int)pts.size());
    unordered_map<LK,unordered_set<int>> line_pts;
    for(int i=0;i<(int)pts.size();i++)for(int j=0;j<i;j++){
        if(pts[i].x==pts[j].x&&pts[i].y==pts[j].y)continue;
        LK k=line_of(pts[i],pts[j]); line_pts[k].insert(i); line_pts[k].insert(j);
    }
    int bad=0;
    for(auto&kv:line_pts){ int t=(int)kv.second.size(); if(t>=3){ bad+=t*(t-1)*(t-2)/6;
        int aa,bb,cc; uint64_t Kk=kv.first; cc=(int)(Kk%CSPAN)-COFF; Kk/=CSPAN; bb=(int)(Kk%BSPAN)-BOFF; Kk/=BSPAN; aa=(int)Kk-AOFF;
        bool truly=false;
        auto it=kv.second.begin(); int i1=*it; ++it; int i2=*it; ++it; int i3=*it;
        int dx=pts[i2].x-pts[i1].x, dy=pts[i2].y-pts[i1].y;
        if(dx*(pts[i3].y-pts[i1].y)==dy*(pts[i3].x-pts[i1].x)) truly=true;
        printf("LINE t=%d (a,b,c)=(%d,%d,%d) truly_collinear=%d: pts",t,aa,bb,cc,truly?1:0);
        for(int idx:kv.second)printf(" (%d,%d)",pts[idx].x,pts[idx].y);
        printf("\n");
    }}
    printf("total_bad(packed)=%d\n",bad);
    return 0;
}
