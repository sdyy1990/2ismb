/*!
 * \file filecombi.h
 * Contains IO utilities.
 */
#pragma once
#include "io_helper.h"
#include <bitset>
/*
template <typename keyType, typename valueType>
class IOHelper {
public:
    virtual bool convert(char *s, keyType *T, valueType *V)  = 0;
    virtual bool convert(char *s, keyType *T)  = 0;
    virtual void splitgrp(const keyType &key, uint32_t &grp, keyType &keyInGroup) = 0;
    virtual void combgrp(keyType &key, uint32_t &grp, keyType &keyInGroup) = 0;
};

template<typename keyType, typename valueType>
class FileReader {
public:
    IOHelper<keyType, valueType> *helper;
    virtual bool getFileIsSorted() = 0;
    virtual bool getNext(keyType *T, valueType *V) = 0;
    virtual void finish() =0;
    virtual void reset() = 0;
    virtual ~FileReader() {
    }
};
*/

template <typename keyType, int NNL>
struct BinaryBitSet {
public:
    unsigned char m[NNL/8];
    keyType k;
    static_assert(NNL % 8 == 0, "multiples of 8");
    unsigned int count() {
        uint64_t u64[NNL/64+1];
        memset(u64,0,sizeof(u64));
        memcpy(u64,m,NNL/8);
        int cnt = 0;
        for (int i = 0 ; i < NNL; i+= 64) {
            cnt += __builtin_popcountll(u64[i>>6]);
        }
        return cnt;
    } 
    void setvalue(unsigned int t) {
        m[t>>3] |= (1<<(t&7));
    }
    void reset() {
        memset(m,0,NNL/8);
    }
    bool operator [] (const int t) {
        return (m[t>>3] >> (t&7) & 1);
    } 
} __attribute__((packed));

template <typename keyType, int NNL>
class KmerGroupedReader{ //: public FileReader <keyType, unsigned char [NNL/8] > {
    vector< FILE *> fV;
    struct KIDpair {
        keyType k;
        uint32_t id;
        bool finished;
        bool friend operator <( const KIDpair &a, const KIDpair &b) {
            if (a.finished != b.finished) return (((int) a.finished) > ((int) b.finished));
            return a.k>b.k;
        }
    };
public:
    void finish() {
        for (auto f: fV) fclose(f);
    }
    void reset() {
        printf(" Do not support reset() \n");
    }
    vector<KmerReader<uint64_t> *> readers;
    vector<MultivalueFileReaderWriter<uint64_t, uint16_t> *> grpreaders; //must be 64-bit kmers, and 16-bit grpids.
    priority_queue<KIDpair> PQ;
    bool combineMode = false; //used when there are >=800 files;
    uint32_t combineCount; // split the file into combineCount groups,
    bool getFileIsSorted() {
        return true;
    }
    void groupFile(string fname, vector<string> lf, string prefix,  int32_t idshift, bool useBinaryKmerFile,uint32_t KmerLength, const char * tmpfolder) {
        vector<KmerReader<keyType> *> readers;
        priority_queue<KIDpair> PQN;
        for (string s: lf) {
            string fname = prefix + s ;
            if (useBinaryKmerFile) 
                readers.push_back(new BinaryKmerReader<keyType>(fname.c_str()));
            else {
                string tmpfname(tmpfolder); tmpfname = tmpfname + s + ".bintmp";
                readers.push_back(new SortedKmerTxtReader<keyType>(fname.c_str(),KmerLength,NULL));
            }
            keyType key; 
            readers[readers.size()-1]->getNext(&key);
            KIDpair kid = {key, idshift+readers.size()-1, false};
            PQN.push(kid);
        }
        
        MultivalueFileReaderWriter<keyType,uint16_t> * writer = new MultivalueFileReaderWriter<keyType,uint16_t> (fname.c_str(),8,2,false); 
        // Loop key for these files;
        while (true) {
            keyType key = PQN.top().k;
            uint32_t id = PQN.top().id;
            vector<uint16_t> ret;
            if (PQN.top().finished) {
                for (auto r: readers) { 
                    r->finish();
                    delete r;
                }
                writer->finish();
                delete writer;
                return;
            }
            while (PQN.top().k == key && !PQN.top().finished) {
                int tid = PQN.top().id;
                ret.push_back(tid);
                keyType nextk;
                bool finish = !readers[tid-idshift]->getNext(&nextk);
                PQN.pop();
                KIDpair kid = {nextk, tid, finish};
                PQN.push(kid);
            }
            writer->write(&key, ret); 
        }
    }
    vector< vector<uint16_t> > grpTmpValue;

    KmerGroupedReader(const char * NCBIfname, const char * fnameprefix, const char * tmpFileDirectory, uint32_t KmerLength, bool useBinaryKmerFile = true ) {
        ConstantLengthKmerHelper<keyType,uint8_t> helper(KmerLength,-1);
        FILE * fNCBI;
        string prefix ( fnameprefix);
        fNCBI = fopen(NCBIfname, "r");
        //Assuming each line of the file contains a filename.
        char buf[4096];
        readers.clear();
        vector<string> fnames;
        while (true) {
            if (fgets(buf, 4096, fNCBI) == NULL) break; // read a Species
            string fname(buf);
	    if (*fname.rbegin() == '\n') fname = fname.substr(0,fname.size()-1);
            fnames.push_back(fname);
        }

        int nn = 500;
        combineMode = (fnames.size()>nn);
        if (combineMode) {
            int curr = 0;
            int combineCount = 0;
            vector<string> * fnamesInThisgrp ;
            vector<string> grpfnames;
            while (curr < fnames.size()) {
                if (curr + nn < fnames.size())
                    fnamesInThisgrp = new vector<string> (fnames.begin()+curr, fnames.begin()+curr+nn);
                else
                    fnamesInThisgrp = new vector<string> (fnames.begin()+curr, fnames.end());
                stringstream ss;
                string tmpFolder(tmpFileDirectory);

                ss<<tmpFolder<<"TMP"<<grpfnames.size();
                string fnamegrp;
                ss>> fnamegrp;
                grpfnames.push_back(fnamegrp);
                printf("merge kmer files %d %d to grp %s\n", curr, curr+fnamesInThisgrp->size()-1, fnamegrp.c_str());
                groupFile(fnamegrp, *fnamesInThisgrp, prefix,  curr, useBinaryKmerFile,KmerLength,tmpFileDirectory);
                curr += fnamesInThisgrp->size();
                delete fnamesInThisgrp;
            }
            combineCount = grpfnames.size();
            for (string v: grpfnames) {
                grpreaders.push_back( new MultivalueFileReaderWriter<uint64_t, uint16_t>(v.c_str(), 8,2, true));
                keyType key;
                uint16_t valuebuf[1024];
                grpreaders[grpreaders.size()-1]->getNext(&key, valuebuf);
                vector<uint16_t> Vvaluebuf;
                for (int i = 0 ; grpreaders[0]->valid(valuebuf[i]); i++)
                   Vvaluebuf.push_back(valuebuf[i]);
                grpTmpValue.push_back(Vvaluebuf);
                KIDpair kid = {key, grpreaders.size()-1, false};
                PQ.push(kid);
            }
        }
        else
            for (int i = 0 ; i < fnames.size(); i++) {
                string fname = prefix + fnames[i] ;
                if (useBinaryKmerFile)
                    readers.push_back(new BinaryKmerReader<keyType>(fname.c_str()));
                else {
                    string tmpfname(tmpFileDirectory); tmpfname = tmpfname + fnames[i] + ".bintmp";
                    readers.push_back(new SortedKmerTxtReader<keyType>(fname.c_str(),KmerLength,tmpfname.c_str()));
                }
                keyType key;
                readers[readers.size()-1]->getNext(&key);
                KIDpair kid = {key, readers.size()-1, false};
                PQ.push(kid);
            }
        fclose(fNCBI);
    }
    ~KmerGroupedReader() {
        if (combineMode)  {
            for (int i = 0 ; i < grpreaders.size(); i++)
                delete grpreaders[i];
        }
        else 
            for (int i = 0 ; i < readers.size(); i++)
                delete readers[i];
    }
    bool getNext( BinaryBitSet<uint64_t, NNL> *kvpair) {
    //bool getNext(keyType *k, unsigned char *v[NNL/8]) {

        int anslevel = 0;
        keyType key = PQ.top().k;
        vector<int> ret;
        if (PQ.top().finished) {
            finish();
            return false;
        }
       // printf("Find key %llx:", key);
        while (PQ.top().k == key && !PQ.top().finished) {
            int tid;
            tid = PQ.top().id;
            keyType nextk;
            bool finish;
            if (combineMode) {
                ret.insert(ret.end(),grpTmpValue[tid].begin(),grpTmpValue[tid].end());
                int ll = grpTmpValue[tid].size();
           //     printf("   %d keys: (from %d)\t", ll, tid);
           //     for (int i: grpTmpValue[tid])
           //          printf("%x\t",i);
                uint16_t valuebuf[1024];
                finish = !grpreaders[tid]->getNext(&nextk, valuebuf);
                grpTmpValue[tid].clear();
                for (int i = 0; grpreaders[tid]->valid(valuebuf[i]); i++)
                    grpTmpValue[tid].push_back(valuebuf[i]);
            //    printf("Next Has ::%d::", grpTmpValue[tid].size());
            }
            else {
                ret.push_back(tid);
              //  printf(" %x\t",PQ.top().id);
                finish = !readers[tid]->getNext(&nextk);
            }
            PQ.pop();
            KIDpair kid = {nextk, tid, finish};
            PQ.push(kid);
        }
        kvpair -> reset();
        kvpair -> k = key;
	    for (auto a: ret) {
            kvpair->setvalue(a);
        }
        if (kvpair->k == 0x089fc04f33) {
            printf("%016llx\n",kvpair->k);
            for (int i = 0 ; i < NNL; i++)
                printf("%d",(*kvpair)[i]);
            printf("\n");
        }
        return true;
    }
};
