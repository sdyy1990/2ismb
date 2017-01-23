#pragma once
#include "othello.h"
#include <unordered_map>
#include "kmercompress.h"
using namespace std;
template <typename keyType, unsigned int NNNL> 
class vNode{
    vector<uint64_t> values;
    vector<keyType> keys;
    public:
    uint64_t keycount() {
        return keys.size();
    }
    class arrayHash {
        public:
              std::size_t operator()(std::array<unsigned char, NNNL/8> const &arr) const {
                 std::size_t sum(0);
                 for (int i = 0; i <NNNL/8; i+=8) {
                     uint64_t *p; p = (uint64_t *)  &arr[i];
                     sum ^= std::hash<uint64_t>()(*p);
                 }
                 return sum;
              }
    };

    unordered_map< array<unsigned char, NNNL/8> , int, arrayHash> linemap;
    vector< array<unsigned char, NNNL/8> > lines;
    
    static const int EMPTY = 1;
    static const int OCTVALUES = 2;
    static const int MAPP = 4;
    public:
    bool isCompliment = false;
    int type= 0;
    int valuecnt;
    int maxn;
    public:
    void definetypes() {
        lines.clear();
        lines.push_back(array<unsigned char, NNNL/8>());
        if (valuecnt == 0) {
            type = EMPTY;
            return;
        }
        if (valuecnt == maxn) {
            type = EMPTY;
            isCompliment = true;
            valuecnt = 0;
            return;
        }
        if (valuecnt <= 8 || valuecnt>=maxn-8) {
            type = OCTVALUES;
            if (isCompliment = (valuecnt>8)) 
                valuecnt = maxn-valuecnt;
            return;
        }
        type = MAPP;
    }
    vNode(int _valuecnt,int _maxn) {
        maxn = _maxn;
        valuecnt = _valuecnt;
        definetypes();
    }
    Othello<keyType> *oth;
    uint64_t getvalue(keyType *k) {
        uint64_t ans = oth->queryInt(*k);
        return ans;
    }
    bitset<NNNL> query(keyType *k) {
        bitset<NNNL> ans;
        if (type == EMPTY) {
            if (isCompliment) ans.flip();
            return ans;
        }
        uint64_t value = getvalue(k);
        if (value == 0) return ans;
        if (type == OCTVALUES) {
            if (value) {
                while (value) {
                    if ( (value & 0xFF) < maxn) 
                        ans.set(value & 0xFF);
                    value >>= 8;
                }
            }
            if (isCompliment) ans.flip();
            return ans;
        }
        if (type == MAPP) {
            array<unsigned char, NNNL/8> pp; pp.fill(0);
            if (value >= lines.size()) return ans;
            memcpy(pp.data(),&lines[value],NNNL/8);
            for (int i = 0 ; i < maxn; i++)
                if ((pp[i>>3]>>(i&7)) & 1) {
                    if (i>=maxn) continue;
                    ans.set(i);
                }
            return ans;
        }
    }
    void addKey(BinaryBitSet<keyType,NNNL> &buf) {
        if (type == OCTVALUES) {
            uint64_t value = 0ULL;
            for (int i = 0 ; i < maxn; i++) 
                if (!isCompliment){
                   if (buf[i]) {
                      value <<= 8;   value += (i & 0xFF);
                    } 
                }
                else{
                   if (!buf[i]) {
                      value <<= 8;   value += (i & 0xFF);
                    } 
                }
            values.push_back(value);
            keys.push_back(buf.k);
            return;
        }
        if (type == MAPP) {
            array<unsigned char, NNNL/8> mm; 
            memcpy(mm.data(), buf.m, NNNL/8);
            if (linemap.find(mm) == linemap.end())  {
                lines.push_back(mm);
                values.push_back(linemap[mm] = lines.size());
             }
            else 
                values.push_back(linemap[mm]);
            keys.push_back(buf.k);
        }
    }
    void construct() {
        int L = 64;
        if (type == OCTVALUES) L = valuecnt * 8;
        if (type == MAPP) {
            L = 8;
            while ((1<<L)<values.size()) L++;
        }
        oth = new Othello<uint64_t> (L, keys,values, true, 0);
        
    }
    void writeDataToBinaryFile(FILE *fout) {
        unsigned char buf[0x20];
        memset(buf,0,sizeof(buf));

        uint32_t linescount = lines.size();
        memcpy(buf, &linescount, 4);
        if (isCompliment) valuecnt = maxn-valuecnt;
        memcpy(buf+4,&valuecnt, 4);
        if (isCompliment) valuecnt = maxn-valuecnt;
        memcpy(buf+8,&maxn, 4);
        fwrite(buf,sizeof(buf),1,fout);
        oth->exportInfo(buf);
        fwrite(buf,sizeof(buf),1,fout);
        oth->writeDataToBinaryFile(fout);
        if (type == MAPP) {
            for (int i = 0 ; i < linescount; i++)
                fwrite(lines[i].data(),NNNL/8,1,fout);
        }
        
    }
    uint32_t linescount;
    vNode(unsigned char * p) {
        memcpy(&linescount, p, 4);
        memcpy(&valuecnt, p+4, 4);
        memcpy(&maxn, p+0x8, 4);
        definetypes();
    }
    void loadDataFromBinaryFile(FILE *fin) {
        unsigned char buf[0x20];
        fread(buf,sizeof(buf),1,fin);
        oth = new Othello<uint64_t> (buf);
        oth->loadDataFromBinaryFile(fin);
        if (type == MAPP) {
            for (int i = 0 ; i < linescount; i++) {
                array<unsigned char, NNNL/8> aa;
                fread(aa.data(), 1, NNNL/8, fin);
                lines.push_back(aa);
            }
        }
    }
};

template <typename keyType, unsigned int NNNL>
class OthelloLargeNode {
public:
    int maxn;
    Othello<keyType> *freqOth = NULL;
    bool kmerCompressed;
    bool kmerUniversial;
    bool kmerHybrid;
    vector< vNode<keyType, NNNL> *> vNodes;
    Othello<keyType> *oth;
    vector<unsigned char> *freqmap;
    static const uint64_t FREQMAPSIZE = 6973568800ULL;
    static const uint64_t MAXHASH = 2145390523ULL;
    unsigned char get8bit(keyType k) {
        uint8_t ans = 0;
        uint8_t * pk = (uint8_t *) &k;
        for (int q = 0; q <8; q++) {
            ans ^= *pk;
            pk++;
        }
        if (ans ==0) return k | 0x1;
        return ans;
    }
    uint64_t gethash(keyType k) {
        __uint128_t vv = k;
        vv = vv * vv;
        return (vv % MAXHASH);
    }
    OthelloLargeNode(KmerGroupedReader<keyType,NNNL> *builder,int _maxn,bool _kmerCompressed, bool _kmerUniversial, bool _isHybrid) {
        kmerCompressed = _kmerCompressed;
        kmerUniversial = _kmerUniversial;
        kmerHybrid = _isHybrid;
        maxn = _maxn;
        vector<uint32_t> vNodeID(256);
        if (kmerUniversial) {
            freqmap = new vector<unsigned char> (FREQMAPSIZE);
            for (int i = 2 ; i <= 8; i++) {
                vNodeID[i] = vNodes.size();
                vNodes.push_back(new vNode<keyType,NNNL>(i,maxn));
            }
            for (int i = maxn-8 ; i <= maxn; i++) {
                vNodeID[i] = vNodes.size();
                vNodes.push_back(new vNode<keyType,NNNL>(i,maxn));
            }
            int offset = vNodes.size();
            int cnt = 255 - maxn - offset;
            while (vNodes.size()<255-maxn) 
                vNodes.push_back(new vNode<keyType,NNNL>(maxn/2,maxn));
            //exactly 128 vNodes.
            for (int i = 9; i<maxn-8; i++)
                vNodeID[i] = (i%cnt) + offset; //<128

        }
        else if (kmerHybrid) {
            freqmap = new vector<unsigned char> (MAXHASH);
            while (vNodes.size()<256) 
                vNodes.push_back(new vNode<keyType,NNNL>(maxn/2,maxn));

        }
        else {
            for (int i = 0; i<= NNNL; i++) 
                vNodes.push_back(new vNode<keyType,NNNL>(i,maxn));
        }
        vector<keyType> kV;
        vector<uint16_t> vV;
       
        BinaryBitSet<keyType, NNNL> buf;
        while (builder->getNext(&buf)) {
            int cnt = buf.count();
            if (kmerUniversial) {
                uint64_t id = compressed20merToID(buf.k);
                if (cnt>1) {
                    (*freqmap)[id] = vNodeID[cnt];
                    vNodes[vNodeID[cnt]] -> addKey(buf);
                }
                else  {
                    int pp = 0; while (buf[pp] == 0) pp++;
                    (*freqmap)[id] = vNodes.size()+pp; //which one;
                }
            }
            else if (kmerHybrid) {
                uint64_t id = gethash(buf.k);
                if ((*freqmap)[id] == 0) {
                    (*freqmap)[id] = get8bit(buf.k);
                }
                vNodes[(*freqmap)[id]] -> addKey(buf);
            }
            else {
                kV.push_back(buf.k);
                vV.push_back(cnt);
                vNodes[cnt]->addKey(buf);
            }
        }

        if ((!kmerUniversial) && (!kmerHybrid)) {
            int LLfreq = 8;
            while ((1<<LLfreq)<NNNL) LLfreq++;
            freqOth = new Othello<keyType> (LLfreq, kV, vV, true,10);
        }
        for (int i = 0 ; i < vNodes.size(); i++)
            printf("%d \t",vNodes[i]->keycount());
        printf("\n");
        for (int i = 0 ; i < vNodes.size(); i++)
            vNodes[i]->construct();

    }

    OthelloLargeNode(FILE *fin) {
        unsigned char buf[0x20];
        memset(buf,0,sizeof(buf));
        fread(buf,sizeof(buf),1,fin);
        memcpy(&maxn,buf,4);
        uint32_t ioflags;
        memcpy(&ioflags,buf+4,4);
        kmerCompressed = (ioflags & 1);
        kmerUniversial = (ioflags & 2);
        kmerHybrid = (ioflags & 4);
        uint32_t nodecount;
        memcpy(&nodecount,buf+8,4);
        if (kmerUniversial || kmerHybrid) {
            uint64_t siz = kmerUniversial?FREQMAPSIZE:MAXHASH;
            freqmap = new vector<unsigned char> (siz);
            fread(&(*freqmap)[0],sizeof (unsigned char),siz,fin);
        }
        else {
            fread(buf,sizeof(buf),1,fin);
            freqOth = new Othello<keyType> ( buf);
            freqOth -> loadDataFromBinaryFile(fin);
        }
        for (int i = 0; i<nodecount; i++) {
            fread(buf,sizeof(buf),1,fin);
            vNodes.push_back(new vNode<keyType,NNNL>(buf));
            vNodes[i]->loadDataFromBinaryFile(fin);
        }
    }

    void exportToFile(FILE *fout) {
        unsigned char buf[0x20];
        memset(buf,0,sizeof(buf));
        memcpy(buf,&maxn,4);
        uint32_t ioflags = 0;
        if (kmerCompressed) ioflags+=1;
        if (kmerUniversial) ioflags+=2;
        if (kmerHybrid) ioflags+=4;
        memcpy(buf+4,&ioflags,4);
        uint32_t nodecount = vNodes.size();
        memcpy(buf+8,&nodecount,4);
        fwrite(buf,1,sizeof(buf),fout);

        if (kmerUniversial || kmerHybrid) {
            fwrite(&(*freqmap)[0],sizeof (unsigned char),freqmap->size(),fout);
        }
        else {
            freqOth -> exportInfo(buf);
            fwrite(buf,1,sizeof(buf),fout);
            freqOth -> writeDataToBinaryFile(fout);
        }
        for (int i = 0; i<vNodes.size(); i++) {
            if (vNodes[i] == NULL) 
                vNodes[i] = new vNode<keyType,NNNL>(0,maxn);
            vNodes[i]->writeDataToBinaryFile(fout);
        }
        
    }
    bitset<NNNL> query(keyType &k) {
        if (k & 0x8000000000000000ULL) return bitset<NNNL>();
        uint32_t othquery;
        if (kmerUniversial) {
            othquery = (*freqmap)[compressed20merToID(k)];
            if (othquery >= vNodes.size()) {
                bitset<NNNL> ans;
                ans.set(othquery - vNodes.size());
                return ans;
            }
            bitset<NNNL> ans = vNodes[othquery]->query(&k);
            return ans;
        }
        else if (kmerHybrid) {
            othquery = (*freqmap)[gethash(k)];
            bitset<NNNL> ans = (othquery>0)?(vNodes[othquery]->query(&k)):bitset<NNNL>();
            return ans;

        }
        else
        othquery = freqOth->queryInt(k);        
  //      printf("oth-> %d \n",othquery);
  //      if (othquery <=maxn) {
  //          printf("iscomplement %d", vNodes[othquery]->isCompliment);
  //      }
        bitset<NNNL> ans = (othquery<=maxn)?(vNodes[othquery]->query(&k)):bitset<NNNL>();
        return ans;
    }
    vector<vector<bitset<NNNL>>> query(vector<vector<keyType>> &req, bool batch) {
        vector<vector<bitset<NNNL>>> ans;
        if (batch) {
            vector< vector< pair<int,int> > > querylist(maxn+1);
            for (int i = 0 ; i < req.size(); i++) 
                for (int j = 0 ; j < req[i].size(); j++) if (req[i][j] < 0x8000000000000000ULL) {
                    uint32_t othquery = freqOth->queryInt(req[i][j]);        
                    if (othquery <=maxn && othquery>0) {
                        querylist[othquery].push_back(make_pair(i,j));
                    }
//                    bitset<NNNL> ans = (othquery<=maxn)?(vNodes[othquery]->query(&k)):bitset<NNNL>();

                }
            for (int i = 0; i < req.size(); i++)
                ans.push_back(vector<bitset<NNNL>>(req[i].size()));
            for (int i = 1; i<=maxn; i++) {
                for (auto &v: querylist[i]) {
                    ans[v.first][v.second] = vNodes[i]->query(&req[v.first][v.second]);
                }
            }
            return ans;
        }
        else {
            for (auto &vk: req) {
                vector<bitset<NNNL>> vbs;
                for (auto &k : vk)
                    vbs.push_back(query(k));
                ans.push_back(vbs);
            }
            return ans;
        }
    }
};

