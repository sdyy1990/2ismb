#pragma once
#include <cstring>
#include <cstdlib>
void compress(char * ori, char * ret, int L) {
    *ret = *ori;
    char * ret0; ret0 = ret;
    for (int i = 1 ; i < L; i++) 
        if (ori[i] != ori[i-1]) 
        {
        ret++;
        *ret = ori[i];
        }
    while (ret< ret0 + L -1) {
        *(ret+1) = *ret;
        ret++;
    }
}
void reversecomp(char *dst, char *ori, int L) {
    char *pp; pp = ori + L;
    dst += (L-1);
    while (ori < pp) {
        if (*ori == 'A') *dst = 'T';
        if (*ori == 'C') *dst = 'G';
        if (*ori == 'G') *dst = 'C';
        if (*ori == 'T') *dst = 'A';
        ori ++;
        dst--;
    }
}
uint64_t compressed20merToID(uint64_t key) {
   uint64_t ans;
   ans = key >> (19*2);
   uint64_t s[21] = {0, 4, 16, 52, 160, 484, 1456, 4372, 13120, 39364, 118096, 354292, 
       1062880, 3188644, 9565936, 28697812, 86093440ULL, 258280324ULL, 774840976ULL, 
           2324522932ULL,6973568800ULL};
   uint8_t v[16] = {4,0,1,2,
                    0,4,1,2,
                    0,1,4,2,
                    0,1,2,4};
   //[39][38][37][36]
   int shift = 0;
   for (int i = 36; i>=0; i-=2) {
        uint8_t vadd =  v[ (key >> i) & 0xF];
        if (vadd ==4) break;
        ans *= 3;
        ans += vadd;
        shift ++;
   }

    return s[shift] + ans;

}
