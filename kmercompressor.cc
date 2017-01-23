#include <cstdio>
#include <cstdint>
#include "kmercompress.h"
int main(int argc, char * argv[]) {
    int kmerlength;
    sscanf(argv[1],"%d",&kmerlength);
    size_t size;
    char *line = NULL;
    char *buf = 0;
    char cmp[256];
    char rev[256];
    char cmprev[256];
    char kmer[256];
    memset(rev,0,sizeof(rev));
    memset(cmp,0,sizeof(rev));
    memset(cmprev,0,sizeof(rev));   
    rev[kmerlength] = '\n';
    while (getline(&buf,&size, stdin)!= -1) {
        memcpy(cmp,buf,32);
        memcpy(cmprev,buf,32);
        compress(buf,cmp, kmerlength);
        reversecomp(rev,buf,kmerlength);
        compress(rev,cmprev, kmerlength);
        if (strcmp(cmp,cmprev)<0)
            printf("%s", cmp);
        else 
            printf("%s", cmprev);
    }
    return 0;
}
