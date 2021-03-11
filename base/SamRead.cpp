/**
 * The implementation of SamRecord.h
 */

#include "SamRead.h"
#include "string.h"

#define MAPPING_QUALITY_UNAVAILABLE 255
#define NO_ALIGNMENT_REFERENCE_NAME "*"
#define NO_ALIGNMENT_START 0
#define UNSET_POSITION 0

const char* GetReferenceName(sam_hdr_t *h, bam1_t *read)
{
    return sam_hdr_tid2name(h, read->core.tid);
}

int GetStart(bam1_t *read)
{
    return read->core.pos;
}

//---how to get the end of a read?
int GetEnd(bam1_t *read)
{
    return bam_endpos(read) - 1;
}

bool IsUnmapped(sam_hdr_t *h, bam1_t *read)
{
    return (read->core.flag & BAM_FMUNMAP) != 0 || GetReferenceName(h, read) == NULL ||
           strcmp(GetReferenceName(h, read), NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
           GetStart(read) == NO_ALIGNMENT_START;
}

int getFragmentLength(bam1_t *read)
{
    return read->core.isize;
}

bool isPaired(bam1_t *read)
{
    return (read->core.flag & BAM_FPAIRED) != 0;
}

bool hasWellDefinedFragmentSize(sam_hdr_t *h, bam1_t *read)
{
    if (getFragmentLength(read) == 0)
        return false;
    if (!isPaired(read))
        return false;
    if (IsUnmapped(h, read) || mateIsUnmapped(h, read))
        return false;
    if (isReverseStrand(read) == mateIsReverseStrand(read))
        return false;
    if (isReverseStrand(read))
        return GetEnd(read) > getMateStart(h, read);
    else
        return GetStart(read) <= getMateStart(h, read) + getFragmentLength(read);
}

bool isReverseStrand(bam1_t *read)
{
    return (read->core.flag & BAM_FREVERSE) != 0;
}

bool mateIsReverseStrand(bam1_t *read)
{
    return (read->core.flag & BAM_FMREVERSE) != 0;
}

bool getMateUnmappedFlag(bam1_t *read)
{
    if (!getReadPairedFlag(read))
        throw "Inappropriate call if not paired read";
    return getProperPairFlagUnchecked(read);
}

bool mateIsUnmapped(sam_hdr_t *h, bam1_t *read)
{
    if (!isPaired(read))
        throw "Cannot get mate information for an unpaired read";
    return getMateUnmappedFlag(read) || getMateReferenceName(h, read) == NULL ||
            strcmp(getMateReferenceName(h, read), NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
            getMateAlignmentStart(read) == NO_ALIGNMENT_START;
}

const char* getMateReferenceName(sam_hdr_t *h, bam1_t *read)
{
    return sam_hdr_tid2name(h, read->core.mtid);
}

bool getReadPairedFlag(bam1_t * read)
{
    return (read->core.flag & BAM_FPAIRED) != 0;
}

bool getProperPairFlagUnchecked(bam1_t * read)
{
    return (read->core.flag & BAM_FPROPER_PAIR) != 0;
}

int getMateStart(sam_hdr_t *h, bam1_t *read)
{
    if(mateIsUnmapped(h, read))
        return UNSET_POSITION;
    return getMateAlignmentStart(read);
}

int getMateAlignmentStart(bam1_t * read)
{
    return read->core.mpos;
}

bool isInsideRead(bam1_t *read, int referenceCoordinate)
{
    return referenceCoordinate >= GetStart(read) && referenceCoordinate <= GetEnd(read);
}