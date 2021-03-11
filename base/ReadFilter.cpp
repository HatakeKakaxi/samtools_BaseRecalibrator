/**
 * The implementation of ReadFilter class
 */
#include <iostream>
#include "ReadFilter.h"
#include "string.h"
#include "SamRead.h"

using namespace std;

ReadFilter::ReadFilter(sam_hdr_t * hdr)
{
    this->readsHeader = hdr;
}

bool ReadFilter::AcceptOrNot(bam1_t *read)
{
    if(!MappingQualityNotZero(read))
        return false;
    if(!MappingQualityAvailable(read))
        return false;
    if(!MappedRead(read))
        return false;
    if(!NotSecondaryAlignment(read))
        return false;
    if(!NotDuplicate(read))
        return false;
    if(!PassesVendorQualityCheck(read))
        return false;
    if(!Wellformed(read))
        return false;

    return true;
}

bool ReadFilter::MappingQualityNotZero(bam1_t *read)
{
    return (read->core.qual != 0);
}

bool ReadFilter::MappingQualityAvailable(bam1_t *read)
{
    return read->core.qual != MAPPING_QUALITY_UNAVAILABLE;
}

bool ReadFilter::MappedRead(bam1_t *read)
{
    return !IsUnmapped(readsHeader, read);
}

bool ReadFilter::NotSecondaryAlignment(bam1_t *read)
{
    return (read->core.flag & BAM_FSECONDARY) == 0;
}

bool ReadFilter::NotDuplicate(bam1_t *read)
{
    return (read->core.flag & BAM_FDUP) == 0;
}

bool ReadFilter::PassesVendorQualityCheck(bam1_t *read)
{
    return (read->core.flag & BAM_FQCFAIL) == 0;
}

bool ReadFilter::Wellformed(bam1_t *read) //TODO: To be completed
{
    return GetStart(read) > 0 ; //---how to get ReadGroup of an alignment?
}