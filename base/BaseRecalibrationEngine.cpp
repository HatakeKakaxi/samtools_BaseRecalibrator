/**
 * The key class to process read in BaseRecalibrator
 */
#include <iostream>
#include <assert.h>
#include "BaseRecalibrationEngine.h"
#include "../sam.h"
#include "StandardCovariateList.h"

#define NUMBEROFCOVARIATE 4
using namespace std;

BaseRecalibrationEngine::BaseRecalibrationEngine(sam_hdr_t *hdr) {
    this->readsHeader = hdr;

    int numReadGroups = sam_hdr_count_lines(hdr, "RG");     //---there may be some bugs here
    assert(numReadGroups >= 1);
    this->covariates = new StandardCovariateList(hdr);
}

void BaseRecalibrationEngine::processRead(bam1_t *read, sam_hdr_t *readsHeader) {

    //TODO: read = transform.apply(originalRead)
    int readLength = read->core.l_qseq;
    //int ** cachedKeys = new int[readLength][NUMBEROFCOVARIATE];
    auto cachedKeys = new int[readLength][NUMBEROFCOVARIATE];
            //TODO: cache the keys array
            //TODO: change the type of keys to vector
    //cout << typeid(readLength).name() << endl;
    delete cachedKeys;
}
