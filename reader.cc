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

template <size_t NNNL> 
struct bitset_comp{
   bool operator () (const bitset<NNNL> &b1, const bitset<NNNL> &b2) {
      //if (b1.to_ulong() != b2.to_ulong())
      //  return b1.to_ulong() < b2.to_ulong();
       for (int i = NNNL-1; i >= 0; i--) 
          if (b1[i] ^ b2[i]) return b1[i];
       return false;
   }
};

int main(int argc, char * argv[]) { 
    if (argc < 1) {
        printf("args: filename \n");
        return 0;
    }
      
	const static int NNNL=160;
    BinaryKmerReader< BinaryBitSet<uint64_t, NNNL> > reader(argv[1]);
    BinaryBitSet<uint64_t, NNNL> buf;
    map< bitset<NNNL> ,int, bitset_comp<NNNL> > count;
    vector< map< bitset<NNNL> ,int, bitset_comp<NNNL> > > Vcount;
    for (int i = 0; i<30; i++) {
        Vcount.push_back(count);
    }
    while (reader.getNext(&buf)){
        // if (buf.v.count()>aa || buf.v.count()<bb) continue;
        Vcount[buf.v.count()/5][buf.v]++;
	    buf.v.reset();
    }
    for (int i = 0 ; i <  30; i++) {
        char buf[30]; memset(buf,0,sizeof(buf));
        sprintf(buf,"out%d-%d.txt",i*5,i*5+5);
        FILE * fout; fout = fopen (buf,"w");
        for (auto a:Vcount[i]) {
            fprintf(fout,"%d,%s,%d\n", a.first.count(), a.first.to_string().c_str() , a.second);
        }
        fclose(fout);
    }
    reader.finish();
}






