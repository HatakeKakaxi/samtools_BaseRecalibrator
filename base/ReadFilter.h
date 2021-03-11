/**
 * The class used to filter reads for BaseRecalibrator, only default filters included now
 * Created by lhh, 2021.3.5
 */

#ifndef READFILTER_H
#define READFILTER_H
#include "htslib/sam.h"

#define MAPPING_QUALITY_UNAVAILABLE 255
#define NO_ALIGNMENT_REFERENCE_NAME "*"
#define NO_ALIGNMENT_START 0



class ReadFilter{
private:
    sam_hdr_t * readsHeader;
public:
    ReadFilter(sam_hdr_t * hdr);

    /**
     * Function to filter the reads
     * @param read
     * @return whether or not this read is accepted
     */
    bool AcceptOrNot(bam1_t * read );

    bool MappingQualityNotZero(bam1_t * read);
    bool MappingQualityAvailable(bam1_t * read);
    bool MappedRead(bam1_t *read);
    bool NotSecondaryAlignment(bam1_t * read);
    bool NotDuplicate(bam1_t * read);
    bool PassesVendorQualityCheck(bam1_t * read);
    bool Wellformed(bam1_t * read);

};


#endif