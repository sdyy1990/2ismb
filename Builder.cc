#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <bitset>
#include "filegrouper.h"
#include <chrono>
#include <inttypes.h>
#include <string>
#include <map>
#include <unordered_map>
#include "othello.h"
#include "olt.h"
using namespace std;
typedef uint64_t keyT;
typedef uint16_t valueT;


int main(int argc, char * argv[]) { 
    if (argc < 5) {
        printf("args: descriptiveFilename KmerFnamePrefix Kmer_length OutputFile workingDirectory Kmer[R/C/U/H](for Raw/Compressed/Universial/Hybrid\n");
        printf("when KmerLength <0, use SortedBinaryFiles\n");
        return 0;
    }    
    string Kmer_length_str = argv[3];
    int Kmer_length = atoi(Kmer_length_str.c_str());
    bool useBinaryKmerFile = (Kmer_length < 0);
    if (Kmer_length<0) Kmer_length = -Kmer_length; 
    IOHelper<keyT,uint16_t> *helper;
    const static int NNNL = 152;
    bool isKmerCompressed = (argv[6][0]=='C' || argv[6][0]=='U');
    bool isKmerUniversial = (argv[6][0]=='U');
    bool isKmerHybrid = (argv[6][0]=='H');
    KmerGroupedReader<uint64_t, NNNL> reader(argv[1],argv[2],argv[5],Kmer_length, useBinaryKmerFile);
   
   OthelloLargeNode<uint64_t, NNNL>  OLN(&reader, 148,isKmerCompressed, isKmerUniversial, isKmerHybrid);
   FILE *fout; fout = fopen64(argv[4],"wb");
   OLN.exportToFile(fout);
   fclose(fout);
   /*
   FILE *fin;
   fin = fopen(argv[4],"rb");
   OthelloLargeNode<uint64_t, NNNL>  OLN2(fin);
   fclose(fin);
    KmerGroupedReader<uint64_t, NNNL> reader2(argv[1],argv[2],argv[5],Kmer_length, useBinaryKmerFile);

    BinaryBitSet<uint64_t, NNNL> buf;
    while (reader2.getNext(&buf)) {
        uint64_t k = buf.k;
        bitset<NNNL> ans1 = OLN.query(k);
        bitset<NNNL> ans2 = OLN2.query(k);
        printf("%s\n%s\n\n",ans1.to_string().c_str(), ans2.to_string().c_str());
    }
    
*/
    /*
    BinaryKmerWriter< BinaryBitSet<uint64_t, NNNL> > writer(argv[4]);
    BinaryBitSet<uint64_t, NNNL> buf;
    map<int,int> count;

   
    vector<keyT> kV;
    vector<uint8_t> vV;
    vector< vNode<uint64_t, NNNL> > vNodes;
    for (int i = 0; i<= NNNL; i++) {
        vNodes.push_back(vNode<uint64_t,NNNL>(i,148));
    }
    

    while (builder.getNext(&buf)) {
        //for (int i = 0 ; i < NNNL; i++) if (buf.v[i]) printf("%d ",i); printf("\n");
	    writer.write(&buf);
        uint8_t cnt; 
        count[cnt = buf.count()]++;
        kV.push_back(buf.k);
        vV.push_back(cnt);
        vNodes[cnt].addKey(buf);
    }
    Othello<uint64_t> *oth; 
    oth = new Othello<uint64_t> (8, kV,vV, true, 0);

    for (int i = 0; i <= NNNL; i++) {
        vNodes[i].construct();
    }

    BinaryBitSet<uint64_t, NNNL> aal;
    printf("%d\n",sizeof(aal));
    for (auto a:count) {
        printf("%d : %d \t", a.first, a.second);
        if (a.first % 5 == 0) printf("\n");
    }
    
    writer.finish();
    printf("\n");
    BinaryKmerReader< BinaryBitSet<uint64_t, NNNL> > reader(argv[4]);
    int nn = 0;
    while (reader.getNext(&buf)) {
        printf("keys: %llx %llx\n", kV[nn], buf.k);
        int count = buf.count();
        for (int i = 0 ; i < NNNL; i++) if (buf[i]) printf("%d ",i); printf(":::: %d\n",count);
        uint32_t othquery = oth->queryInt(buf.k);        
        bitset<NNNL> ans = (othquery<=maxn)?(vNodes[othquery].query(&buf.k)):bitset<NNNL>();
        count = 0;
        for (int i = 0 ; i < NNNL; i++) {
            if (ans[i]) {               printf("%d ",i); count ++; }
        }
        printf(":::: %d \n", count);
        nn++;
    }

    */

}






