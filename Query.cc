#include <iostream>

#include "filegrouper.h"
#include "olt.h"
#include "kmercompress.h"

int main(int argc, char * argv[]) {
    if (argc < 3) {
        printf("args: NodeFile rawQueryFile Kmer_length nodeIsCompressed[Y/N]=N UseBatchQuery[Y/N]=N DetailedMap[Y/N]=N\n");
    }
    int Kmer_length = atoi(argv[3]);
    FILE *fin;
    fin = fopen64(argv[2],"rb");
    bool nodeIsCompressed = ((argc >4) && argv[4][0] == 'Y');
    bool useBatchQuery = ((argc >5) && argv[5][0] == 'Y');
    bool detailedMap = ((argc >6) && argv[6][0] == 'Y');
    char buf[1048576];
    vector<string> vSeq;
    while ( fgets(buf,sizeof(buf),fin)!= NULL) {
        char * p; p = &buf[0];
        while (*p == 'A' || *p == 'T' || *p == 'G' || *p == 'C' || *p == 'N') *p++;
        *p = '\0';
        if (strlen(buf)>=Kmer_length) 
            vSeq.push_back(string(buf));
    }
    fclose(fin);
    //raw no batch query
    memset(buf,0,sizeof buf); 
    char buf2[32];
    char buf3[32];
    memset(buf2 , 0 , sizeof buf2);
    memset(buf3 , 0 , sizeof buf3);
    ConstantLengthKmerHelper<uint64_t, uint16_t> helper(Kmer_length,0);
    fin = fopen(argv[1],"rb");
    const static int NNNL = 152;
    OthelloLargeNode<uint64_t, NNNL>  OLN(fin);
    fclose(fin);
    vector<vector<uint64_t> > requests;

    for (auto &str: vSeq)  {
        vector<uint64_t> vu64(str.size()-Kmer_length + 1);
        for (int i = 0 ; i < str.size() - Kmer_length + 1; i++) {
            memcpy(buf,str.data()+i,Kmer_length);
            uint64_t key;
            if (nodeIsCompressed) {
                compress(buf,buf2,Kmer_length);
                helper.convert(buf2,&key);
                reversecomp(buf2,buf,Kmer_length);
                compress(buf2,buf3,Kmer_length);
                uint64_t key2;
                helper.convert(buf3,&key2);
                vu64[i] = (key2<key)?key2:key;
            }
            else {
                helper.convert(buf,&key);
                vu64[i] = helper.minSelfAndRevcomp(key);
            }
        }
        requests.push_back(vu64);
    }
    //uint64_t kkk = 0x89fc04f33; 
    //uint64_t kkk = 0x002271333; 
    //bitset<NNNL> VV = OLN.query(kkk);
    //printf("%s\n",VV.to_string().c_str());
    //return 0;
    vector< vector<bitset<NNNL> > > res = OLN.query(requests,useBatchQuery);
    for (int t = 0 ; t < vSeq.size(); t++) {
        printf("%s\n",vSeq[t].c_str());
        vector<int> ans (NNNL);
        for (auto &bs : res[t]) 
            for (int i = 0 ; i < NNNL; i++)
                if (bs[i])
                    ans[i] ++;
        for (int i = 0 ; i < NNNL; i++)
            printf("%d\n",ans[i]);
        if (detailedMap) 
            for (int i = 0 ; i < res[t].size(); i++) {
                string sss = res[t][i].to_string();
                reverse(sss.begin(),sss.end());
                printf("%s:%s:%012llx:%d\n",sss.c_str(), vSeq[t].substr(i,Kmer_length).c_str(),requests[t][i],res[t][i].count());
            }
    }
    
}

