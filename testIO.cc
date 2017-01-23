#include "io_helper.h"
#include <string>
#include "kmercompress.h"
using namespace std;
int main(int argc, char * argv[]) {
    printf("%016llx\n",compressed20merToID( 0x11113b7264ULL));
    ConstantLengthKmerHelper<uint64_t, uint32_t> helper(20,0);
    string aa(argv[1]);
    aa.resize(20);
    char buf1[40];
    char buf2[40];
    uint64_t k;
    strcpy(buf2,aa.c_str());
    helper.convert(&(buf2[0]),&k);
    memset(buf1,0,sizeof(buf1));
    memset(buf2,0,sizeof(buf2));
    helper.convertstring(buf1,&k);
    uint64_t k2 = helper.reverseComplement(k);
    helper.convertstring(buf2,&k2);
    printf("%s+\n%s+\n%016llx \n%016llx\n", buf1, buf2,k,k2);
    return 0;



}
