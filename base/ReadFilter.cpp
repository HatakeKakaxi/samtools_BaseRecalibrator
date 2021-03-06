/**
 * The implementation of ReadFilter class
 */
#include <iostream>
#include "ReadFilter.h"
#include "string.h"

using namespace std;

bool ReadFilter::AcceptOrNot(sam_hdr_t *h, bam1_t *read)
{
    if(!MappingQualityNotZero(read))
        return false;
    if(!MappingQualityAvailable(read))
        return false;
    if(!MappedRead(h, read))
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

bool ReadFilter::MappedRead(sam_hdr_t *h, bam1_t *read)
{
    return !IsUnmapped(h, read);
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

bool ReadFilter::IsUnmapped(sam_hdr_t *h, bam1_t *read)
{
    return (read->core.flag & BAM_FMUNMAP) != 0 || GetReferenceName(h, read) == NULL ||
        strcmp(GetReferenceName(h, read), NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
        GetStart(read) == NO_ALIGNMENT_START;
}

const char* ReadFilter::GetReferenceName(sam_hdr_t *h, bam1_t *read)
{
    return sam_hdr_tid2name(h, read->core.tid);
}

//---this function needs to be confirmed
int ReadFilter::GetStart(bam1_t *read)
{
    return read->core.pos + 1;    //---need to +1
}