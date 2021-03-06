/**
 *  Created by lhh, 2021.3.2
 * */

#ifndef BASERECALIBRATION_ENGINE_H
#define BASERECALIBRATION_ENGINE_H

#include "../sam.h"
#include "StandardCovariateList.h"

class BaseRecalibrationEngine{
private:
    sam_hdr_t * readsHeader;
    StandardCovariateList * covariates;
public:
    BaseRecalibrationEngine(sam_hdr_t *hdr);

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     * @param read:  The key structure for one alignment, defined in sam.h
     * @param readsHeader:  The key structure for header of sam file
     */
    void processRead(bam1_t *read, sam_hdr_t *readsHeader);


};

#endif