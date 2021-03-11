/**
  * A wrapper for sam record
  */

#ifndef SAM_READ_H
#define SAM_READ_H

#include "htslib/sam.h"

/**
 * Only to encapsulate sam_hdr_tid2name function in sam.h
 */
const char* GetReferenceName(sam_hdr_t *h, bam1_t *read);

/**
 * zero-based start
 */
int GetStart(bam1_t *read);

/**
 * @return zero-based end
 */
int GetEnd(bam1_t *read);

int getMateStart(sam_hdr_t *h, bam1_t *read);

/**
 * Does the read have a position assigned to it for sorting purposes.
 * @return `true if this read has no assigned position or contig.
 */
bool IsUnmapped(sam_hdr_t *h, bam1_t *read);

/**
 * @return True if this read's mate is unmapped (this includes mates that have a position but are explicitly marked as unmapped,
 *         as well as mates that lack a fully-defined position but are not explicitly marked as unmapped). Otherwise false.
 * @throws IllegalStateException if the read is not paired (has no mate)
 */
bool mateIsUnmapped(sam_hdr_t *h, bam1_t *read);

/**
 * Returns the observed length of the read's fragment (equivalent to TLEN in SAM).
 *
 * Warning: the precise meaning of this field is implementation/technology dependent.
 *
 * @return The observed length of the fragment (equivalent to TLEN in SAM), or 0 if unknown.
 *         Negative if the mate maps to a lower position than the read.
 */
int getFragmentLength(bam1_t * read);

bool isPaired(bam1_t * read);

/**
 * Can the adaptor sequence of read be reliably removed from the read based on the alignment of
 * read and its mate?
 *
 * @param read the read to check
 * @return true if it can, false otherwise
 */
bool hasWellDefinedFragmentSize(sam_hdr_t *h, bam1_t * read);

bool getReadPairedFlag(bam1_t * read);

bool getProperPairFlagUnchecked(bam1_t * read);

bool getMateUnmappedFlag(bam1_t *read);

const char* getMateReferenceName(sam_hdr_t *h, bam1_t *read);

int getMateAlignmentStart(bam1_t * read);

/**
 * @return True if this read is on the reverse strand as opposed to the forward strand, otherwise false.
 */
bool isReverseStrand(bam1_t *read);

bool mateIsReverseStrand(bam1_t *read);

/**
 * Is a base inside a read?
 *
 * @param read                the read to evaluate
 * @param referenceCoordinate the reference coordinate of the base to test
 * @return true if it is inside the read, false otherwise.
 */
bool isInsideRead(bam1_t *read, int referenceCoordinate);

#endif