/**
 * The implementation of class Covariate , ContextCovariate and CycleCovariate
 * */

#include <iostream>
#include "Covariate.h"
#include "../sam.h"

using namespace std;

ContextCovariate::ContextCovariate() {
    this->mismatchesContextSize = 2;
    this->lowQualTail = 2;
}

void ContextCovariate::recordValues(bam1_t *read, sam_hdr_t *readsHeader, auto keys) {
    //cout << "This is ContextCovariate" << endl;
    int originalReadLength = read->core.l_qseq;

    //---get the bases first
    uint8_t * base = bam_get_seq(read);
    //cout << base << endl;
    //uint8_t * strandedClippedBases = getStrandedClippedBytes(base, lowQualTail);
}

char* ContextCovariate::getStrandedClippedBytes(uint8_t * base, char lowQualTail) {
    return NULL;
}

void ContextCovariate::printReads(bam1_t *read)
{
    int ReadLength = read->core.l_qseq;
    uint8_t * base = bam_get_seq(read);

    char RealBase;
    for(int i=0; i<ReadLength; i++)
    {

        char basei = bam_seqi(base, i);
        switch (basei) {
            case 1:
                RealBase = 'A';
                break;
            case 2:
                RealBase = 'C';
                break;
            case 4:
                RealBase = 'G';
                break;
            case 8:
                RealBase = 'T';
                break;
            case 15:
                RealBase = 'N';
                break;
            default:
                RealBase = '-';
                break;
        }
        cout << RealBase;
    }
    cout << endl;
}

