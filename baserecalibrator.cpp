#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include "sam_opts.h"
#include "htslib/hfile.h"
#include "htslib/faidx.h"
#include "base/BaseRecalibrationEngine.h"
#include "base/ReadFilter.h"
#include <assert.h>
#include <iostream>
using namespace std;

extern "C" int main_baserecalibrator(int argc, char *argv[]);

void error(const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fprintf(stderr, "htsfile: ");
    vfprintf(stderr, format, args);
    if (err) fprintf(stderr, ": %s\n", strerror(err));
    else fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
    //status = EXIT_FAILURE;
}

static int usage(FILE *fp, int exit_status, int is_long_help)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools baserecalibrator -I <input.bam> -O <output.table>\n");
    return exit_status;
}

static htsFile *dup_stdout(const char *mode)
{
    int fd = dup(STDOUT_FILENO);
    hFILE *hfp = (fd >= 0)? hdopen(fd, mode) : NULL;
    return hfp? hts_hopen(hfp, "-", mode) : NULL;
}

//---Only for test
void printReads(bam1_t *read)
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


int main_baserecalibrator(int argc, char *argv[])
{
    int returned=0;
    int c;
    int ret;
    char *in=0, *output = 0, *ref = 0;

    //---some bug here
    /*
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
            SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 'T', '@'),
            {"no-PG", no_argument, NULL, 1},
            { NULL, 0, NULL, 0 }
    };*/

    if (argc == 1 && isatty(STDIN_FILENO))
        return usage(stdout, EXIT_SUCCESS, 0);

    while ((c = getopt(argc, argv, "I:O:R:")) >= 0) { //TODO: make it more robust
        switch (c) {
            case 'I': in = strdup(optarg); break;
            case 'O': output = strdup(optarg); break;
            case 'R': ref = strdup(optarg); break;
            default:
                //if (parse_sam_global_opt(c, optarg, lopts, &ga) != 0)
                    return usage(stderr, EXIT_FAILURE, 0);
                //return usage();
            case '?': return usage(stderr, EXIT_FAILURE, 0); break;
        }
    }

    assert(ref != NULL);

    struct hFILE *fp = hopen(in, "r");
    assert(fp != NULL);
    htsFile *hts = hts_hopen(fp, in, "r");

    bam1_t *b = NULL;
    sam_hdr_t *hdr = NULL;


    hdr = sam_hdr_read(hts);
    if (hdr == NULL) {
        errno = 0; error("reading headers from \"%s\" failed", in);
        //goto clean;
        sam_hdr_destroy(hdr);
        bam_destroy1(b);
    }

    ReadFilter filter(hdr);
    BaseRecalibrationEngine baseRecalibrationEngine(hdr, ref);

    b = bam_init1();
    if (b == NULL)
    {
        error("can't create record");
        sam_hdr_destroy(hdr);
        bam_destroy1(b);
        return returned;
    }

    int ReadNum = 0;

    while ((ret = sam_read1(hts, hdr, b)) >= 0) {
        if(filter.AcceptOrNot(b)){
            ReadNum ++;
            //printReads(b);
            //---start of BaseRecalibrationEngine
            baseRecalibrationEngine.processRead(b);
        }

    }

    if (ret < -1) { error("reading \"%s\" failed", in); goto clean; }

    clean:
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    //if (out) hts_close(out);

    return returned;
}