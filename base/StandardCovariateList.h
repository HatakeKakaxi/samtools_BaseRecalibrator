#ifndef STANDARDCOVARIATE_LIST_H
#define STANDARDCOVARIATE_LIST_H

#include <vector>
#include "Covariate.h"
#include "../sam.h"
using namespace std;

class StandardCovariateList{
private:
    vector<Covariate *> additionalCovariates;
    sam_hdr_t * readsHeader;
public:
    StandardCovariateList(sam_hdr_t * hdr);
};

#endif
