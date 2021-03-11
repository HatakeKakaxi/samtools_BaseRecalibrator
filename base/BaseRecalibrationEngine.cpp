/**
 * The key class to process read in BaseRecalibrator
 */
#include <iostream>
#include <assert.h>
#include <vector>
#include <math.h>
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "BaseRecalibrationEngine.h"
#include "htslib/sam.h"
#include "StandardCovariateList.h"
#include "SamRead.h"

#define NUMBEROFCOVARIATE 4
#define CANNOT_COMPUTE_ADAPTOR_BOUNDARY -2147483648
using namespace std;

BaseRecalibrationEngine::BaseRecalibrationEngine(sam_hdr_t *hdr, char * ref) {
    this->readsHeader = hdr;
    this->fai = fai_load3_format(ref, NULL, NULL, FAI_CREATE, FAI_FASTA);
    assert(fai != NULL);

    int numReadGroups = sam_hdr_count_lines(hdr, "RG");     //---there may be some bugs here
    assert(numReadGroups >= 1);
    this->covariates = new StandardCovariateList(hdr);


}

void BaseRecalibrationEngine::processRead(bam1_t *originalRead) {

    //TODO: complete the transform operation
    bam1_t * read = transform(originalRead);
    int readLength = read->core.l_qseq;
    //int ** cachedKeys = new int[readLength][NUMBEROFCOVARIATE];
    auto cachedKeys = new int[readLength][NUMBEROFCOVARIATE];
            //TODO: cache the keys array
            //TODO: change the type of keys to vector


    double * snpErrors = new double[readLength];
    //int * snpErrors = new int[readLength];
    memset(snpErrors, 0, readLength* sizeof(double));
    calculateSNPFractionalError(read , snpErrors);

    delete cachedKeys;
    delete snpErrors;
}

bam1_t *  BaseRecalibrationEngine::transform(bam1_t * originalRead){
    bam1_t * read = consolidateCigar(originalRead);
    read = hardClipAdaptorSequence(read);
    return read;
}

bam1_t * BaseRecalibrationEngine::hardClipAdaptorSequence(bam1_t *originalRead)
{
    int adaptorBoundary = getAdaptorBoundary(originalRead);
    //cout << adaptorBoundary << endl;
    if (adaptorBoundary == CANNOT_COMPUTE_ADAPTOR_BOUNDARY || !isInsideRead(originalRead, adaptorBoundary))
    {
        return originalRead;
    }

    return originalRead;

    //TODO: complete the transform oeration 2021.3.10
    //return isReverseStrand(originalRead) ? hardClipByReferenceCoordinatesLeftTail(adaptorBoundary) : hardClipByReferenceCoordinatesRightTail(adaptorBoundary);
}

int BaseRecalibrationEngine::getAdaptorBoundary(bam1_t *read)
{
    if(!hasWellDefinedFragmentSize(readsHeader, read))
        return CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
    else if (isReverseStrand(read))
        return getMateStart(readsHeader, read) - 1;
    else {
        int insertSize = abs(getFragmentLength(read));
        return GetStart(read) + insertSize;
    }
}


void BaseRecalibrationEngine::setDefaultBaseQualities(bam1_t * originalRead){
    return;
}

bam1_t * BaseRecalibrationEngine::consolidateCigar(bam1_t *originalRead) {
    uint32_t* cigar = bam_get_cigar(originalRead);
    if(!cigar)
        throw "cigar cannot be null";

    // fast check to determine if there's anything worth doing before we create new Cigar and actually do some work
    if (!needsConsolidate(originalRead))
        return originalRead;

    return originalRead;
    //---2021.3.9 to becontinued    //---TODO: realloc a new space for a new bam1_t

    //uint32_t * NewCigar = 0;
    //vector<uint32_t> NewCigar;
    int sumLength = 0;
    int lastElement = -1;
    for(int i=0; i < originalRead->core.n_cigar; i++)
    {
        if(bam_cigar_oplen(cigar[i]) == 0)
            continue;

        if (lastElement != -1 && lastElement != bam_cigar_op(cigar[i]))
        {

        }
    }

}

bool BaseRecalibrationEngine::needsConsolidate(bam1_t * originalRead) {
    if(originalRead->core.n_cigar <= 1)
        return false;

    int lastOp = -1;
    uint32_t* cigar = bam_get_cigar(originalRead);

    for(int icig=0; icig < originalRead->core.n_cigar; icig++)
    {
        if(bam_cigar_oplen(cigar[icig]) == 0 || lastOp == bam_cigar_op(cigar[icig]))
            return true;
        lastOp = bam_cigar_op(cigar[icig]);
    }

    return false;
}

void BaseRecalibrationEngine::calculateSNPFractionalError(bam1_t *read, double *snpErrors)
{
    const char * refName = sam_hdr_tid2name(readsHeader, read->core.tid);
    hts_pos_t start = read->core.pos; //---zero based
    hts_pos_t end = bam_endpos(read)-1;

    hts_pos_t fai_ref_len;
    char * refBases;

    //---get the subsequence of fasta file
    refBases = faidx_fetch_seq64(fai, refName, start, end, &fai_ref_len);
    assert(refBases != NULL);

    int readPos = 0;
    int refPos = 0;
    int elementLength;
    for(int icig=0; icig < read->core.n_cigar; icig++)
    {
        elementLength = bam_cigar_oplen(bam_get_cigar(read)[icig]);
        int CigarOperator = bam_cigar_op(bam_get_cigar(read)[icig]);
        char c1, c2;
        switch (CigarOperator) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for (int i=0; i<elementLength; i++){
                    c1 = bam_seqi(bam_get_seq(read), readPos);
                    c2 = seq_nt16_table[(unsigned char)refBases[refPos]];
                    snpErrors[readPos] = (c1 == c2 ? 0 : 1);
                    readPos ++;
                    refPos ++;
                }
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                refPos += elementLength;
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                readPos += elementLength;
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
                break;
            default:
                throw "Unsupported cigar operator";
        }
    }


    /*
     * For debugging
     *
    for(int i=0; i<read->core.l_qseq; i++ )
    {
        if(snpErrors[i] != 0)
            printf("1");
        else
            printf("0");
        //printf("%d\t", snpErrors[i]);
    }
    printf("\n");   */

    free(refBases);
}
