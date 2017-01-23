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
using namespace std;
typedef uint64_t keyT;
typedef uint16_t valueT;


int main(int argc, char * argv[]) { 
    if (argc < 5) {
        printf("args: descriptiveFilename KmerFnamePrefix Kmer_length OutputFile workingDirectory\n");
        return 0;
    }    
    string Kmer_length_str = argv[3];
    int Kmer_length = atoi(Kmer_length_str.c_str());
    
    IOHelper<keyT,uint16_t> *helper;
    const static int NNNL = 160;
    KmerGroupedReader<uint64_t, uint16_t , 160> builder(argv[1],argv[2],argv[5],Kmer_length, false);
    bitset<NNNL> v;
    uint64_t k;
    BinaryKmerWriter< BinaryBitSet<uint64_t, NNNL> > writer(argv[4]);
    BinaryBitSet<uint64_t, NNNL> buf;
    map<int,int> count;
    while (builder.getNext(&buf.k,(uint16_t *) (&buf.v))) {
        //for (int i = 0 ; i < NNNL; i++) if (buf.v[i]) printf("%d ",i); printf("\n");
	writer.write(&buf);
        count[buf.v.count()]++;
	buf.v.reset();
    }
    for (auto a:count) {
        printf("%d : %d \t", a.first, a.second);
    }
    writer.finish();
}






