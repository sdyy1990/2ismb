#include <cstdio>
#include "io_helper.h"
#include <vector>
#include <algorithm>
using namespace std;
int main(int argc, char * argv[]) {
    int kmerlength;
    if (argc<2) {
        printf("SortAndCompress [text kmer filename] [binary kmer filename] [kmerlength]\n");
    }
    sscanf(argv[3],"%d", &kmerlength);
    ConstantLengthKmerHelper<uint64_t, uint16_t> iohelper(kmerlength,0);
    KmerFileReader< uint64_t,uint16_t > freader(argv[1], &iohelper,false);
    KVpair<uint64_t, uint16_t> pp;
    vector<uint64_t> VKmer;
    while (freader.getNext(&pp.k, &pp.v)) {
        VKmer.push_back(pp.k);
    }
    printf("Sorting\n");
    sort(VKmer.begin(),VKmer.end());
    BinaryKmerWriter< uint64_t > fwriter(argv[2]);
    for (int i = 0; i < VKmer.size(); i++ ) {
        if (i>0)
            if (VKmer[i] == VKmer[i-1]) continue;
        fwriter.write(&VKmer[i]);
    }
    fwriter.finish();
   return 0; 

}
