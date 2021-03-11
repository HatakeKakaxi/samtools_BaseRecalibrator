/**
 *  Created by lhh, 2021.3.2
 * */

#ifndef BASERECALIBRATION_ENGINE_H
#define BASERECALIBRATION_ENGINE_H

#include <map>
#include "../sam.h"
#include "StandardCovariateList.h"
#include "htslib/faidx.h"

class BaseRecalibrationEngine{
private:
    sam_hdr_t * readsHeader;
    faidx_t * fai;
    StandardCovariateList * covariates;
    std::map<char, int> baseIndexMap;


    /**
    * Need a well-formed, consolidated Cigar string so that the left aligning code works properly.
    * For example, 1M1M1M1D2M1M --> 3M1D3M
    * If the given cigar is empty then the returned cigar will also be empty
    *
    * Note that this routine collapses cigar elements of size 0, so 2M0M => 2M
    *
    * @param c the cigar to consolidate
    * @return  a non-null cigar with consecutive matching operators merged into single operators.
    */
    bam1_t * consolidateCigar(bam1_t * originalRead);

    /**
     * Does the cigar need to be consolidated ?
     * @param originalRead
     * @return
     */
    bool needsConsolidate(bam1_t * originalRead);

    /**
     * Optional operation, depending on the arguments
     * @param originalRead
     */
    void setDefaultBaseQualities(bam1_t * originalRead);

    /**
     * Checks if a read contains adaptor sequences. If it does, hard clips them out.
     * Note: To see how a read is checked for adaptor sequence see ReadUtils.getAdaptorBoundary()
     *
     * @return a new read without adaptor sequence (Could return an empty, unmapped read)
     */
    bam1_t * hardClipAdaptorSequence(bam1_t * originalRead);

    int getAdaptorBoundary(bam1_t * read);


public:
    BaseRecalibrationEngine(sam_hdr_t *hdr, char * ref);

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     * @param read:  The key structure for one alignment, defined in sam.h
     * @param readsHeader:  The key structure for header of sam file
     */
    void processRead(bam1_t *read);

    /**
     * Transform the read
     * @param originalRead
     * @return
     */
    bam1_t * transform(bam1_t * originalRead);



    /**
     * Calculate the snpErrors
     * @param read
     * @param snpErrors : The array to contain the snp Error
     */
    void calculateSNPFractionalError(bam1_t *read, double * snpErrors);


};

#endif