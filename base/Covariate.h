/**
 *  Created by lhh, 2021.3.3
 * */

#ifndef COVARIATE_H
#define COVARIATE_H

#include "../sam.h"
using namespace std;

class Covariate{
public:
    void recordValues(bam1_t *read, sam_hdr_t *readsHeader, auto keys);    //---what's the type of keys?
};

class ContextCovariate: public Covariate{
private:
    int mismatchesContextSize;
    char lowQualTail;   //---in Java, its type is byte

    /**
     * print a read alignment to stdout, for debugging
     * @param read
     */
    void printReads(bam1_t *read);
public:
    ContextCovariate();
    void recordValues(bam1_t *read, sam_hdr_t *readsHeader, auto keys);
    char* getStrandedClippedBytes(uint8_t * base, char lowQualTail);
};

#endif