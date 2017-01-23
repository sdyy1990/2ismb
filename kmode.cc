#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <bitset>
#include <cstring>
#include <algorithm>
#include <map>
using namespace std;
const static int NNL = 160;
typedef bitset<NNL> BS;
vector<BS> Vmode;
vector< vector<BS> > Vcluster;
vector<vector<int> > Vcount1;
vector<int> Vsize;
double tmp1, tmp2;
void printgroupinfo(vector<BS> & grp, bool shorter) {
    vector<int> count(NNL,0);
    for (auto point: grp) {
        for (int i = 0 ; i < NNL; i++)
            if (point[i]) count[i]++;
    }
    map<int,int> mmp;
    for (int i = 0 ; i < NNL; i++) mmp[count[i]] ++;
    printf("siz=%d, vis=%.3lf ", grp.size(), (mmp[0]+mmp[grp.size()])*1.0/NNL);
    tmp1 += grp.size(); tmp2 += grp.size()*(mmp[0]+mmp[grp.size()])*1.0/NNL;
    if (shorter) printf("\t");
    else {
        for (auto a:mmp)
            printf("%d:%d\t", a.first, a.second);
        printf("\n");
    }
}

BS mean(vector<BS> & grp) {
    if (grp.size() ==0) {
        BS a; return a;
    }
    return grp[rand()% grp.size()];
    vector<int> count(NNL,0);
    for (auto point: grp) {
        for (int i = 0 ; i < NNL; i++)
            if (point[i]) count[i]++;
    }
    BS ret;
    for (int i = 0 ; i < NNL; i++)
        if (count[i] >= grp.size()/2) ret.set(i);
    return ret;
}
void printclusterinfo( vector< vector< BS> > & vcluster, bool shorter = false) {
    printf("clustercount = %d\n", vcluster.size());
    tmp1 = tmp2 = 0.0;
    for (int i = 0 ; i < vcluster.size(); i++)
        printgroupinfo(vcluster[i], shorter);
    printf("\n");
    printf("%.3lf\n", tmp2/tmp1);
}
int BCCNT;
void printall(char * fname, vector< vector< BS> > & vcluster ) {
    char buf[1024]; memset(buf,0,sizeof(buf));
    sprintf(buf,"%s.%d.res",fname,BCCNT);
    string sfname (buf);
    FILE * fout; fout = fopen(sfname.c_str(),"w");
    for (int i = 0 ; i < vcluster.size(); i++) {
        fprintf(fout,"\n");
        for (auto a: vcluster[i]) 
            fprintf(fout,"%s\n", a.to_string().c_str());
    }
    fclose(fout);
}
double dist(BS &p, vector<int> cnt, int siz) {
    if (siz ==0) {
        double ans = 0.0;
        for (int i = 0 ; i < NNL; i++)
            if (((int) p[i])!=cnt[i]) ans +=1.0;
        return ans;
    }
    double ans = 0.0;
    for (int i = 0; i< NNL; i++) {
        int cmp = cnt[i];
        if (!p[i]) cmp = siz - cmp; //for 0 digits, compare 0 with 0.
        if (cmp == 0) ans += 0.5; // which means all 0, and here comes a 1, hence we need more digits;
        ans += (siz-cmp)*1.0/siz;
    }
    return ans;
}
int main(int argc, char * argv[]) {
    if (argc<2) {
        printf("kmer-bitarry-fname  clustercount \n 1count,000,abundance NNL=%d\n",NNL);
        return 0;
    }
    int CCNT; sscanf(argv[2],"%d",&CCNT);
    BCCNT = CCNT;
    FILE *fin; fin = fopen(argv[1],"r");
    char buf[1024];
    char s[1024];
    vector< bitset<NNL> > VBS;
    while (fgets(buf,1024,fin)) {
        if (strlen(buf)<NNL) continue;
        int cnt, cnt2;
        sscanf(buf,"%d,%[01],%d",&cnt,s,&cnt2);
        bitset<NNL> BS;
        for (int i = 0; i  < NNL; i++)
            if (s[i]=='1') BS.set(NNL-i-1);
        VBS.push_back(BS);
    }
    CCNT = VBS.size()/CCNT; 
    random_shuffle(VBS.begin(),VBS.end());
    cout << VBS.size() << endl;
    
    for (int i = 0 ; i < CCNT; i++) {
        Vmode.push_back(VBS[i]);
        vector<int> count(NNL,0);
        for (int j = 0 ; j < NNL; j++)
            if (Vmode[i][j]) count[j]++;
        Vcount1.push_back(count);
        vector<BS> cluster; //cluster.push_back(VBS[i]);
        Vcluster.push_back(cluster);
        Vsize.push_back(0);
    }
    //round
     {
    for (auto point : VBS) {
        double min = 1e9;
        int choice = -1;
        for (int i = 0; i < CCNT; i++) {
            double dd = dist(point, Vcount1[i], Vsize[i]);
            if (min > dd) {
                min = dd;
                choice = i;
            }
        }
        Vsize[choice]++;
        Vcluster[choice].push_back(point);
        for (int i = 0 ; i < NNL; i++)
            if (point[i]) Vcount1[choice][i]++;
    }
    printclusterinfo(Vcluster,false);
    printall(argv[1], Vcluster);

    for (int i = 0 ; i < Vcluster.size(); i++) {
        Vmode[i] = mean(Vcluster[i]);
        vector<int> count(NNL,0);
        for (int j = 0 ; j < NNL; j++)
            if (Vmode[i][j]) count[j]++;
        Vcount1[i] = count; //push_back(count);
        vector<BS> cluster; //cluster.push_back(VBS[i]);
        Vcluster[i] = (cluster);
        Vsize[i] = 0; //.push_back(0);
        
    }
    random_shuffle(VBS.begin(),VBS.end());

    }
}
