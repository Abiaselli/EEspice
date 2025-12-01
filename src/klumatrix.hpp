#pragma once
#include <iostream>
#include <klu.h>

struct KLUmatrix
{
    klu_l_common KLUmatrixCommon ;                  /* KLU common object */
    klu_l_symbolic *KLUmatrixSymbolic ;             /* KLU symbolic object */
    klu_l_numeric *KLUmatrixNumeric ;               /* KLU numeric object */
    int *KLUmatrixAp ;                              /* KLU column pointer */
    int *KLUmatrixAi ;                              /* KLU row pointer */
    double *KLUmatrixAx ;                           /* KLU Real Elements */
    double *KLUmatrixAxComplex ;                    /* KLU Complex Elements */
    bool KLUmatrixIsComplex ;                       /* KLU Matrix Is Complex Flag */

    KLUmatrix();
    ~KLUmatrix();
};

KLUmatrix::KLUmatrix()
    : KLUmatrixSymbolic(nullptr)
    , KLUmatrixNumeric(nullptr)
    , KLUmatrixAp(nullptr)
    , KLUmatrixAi(nullptr)
    , KLUmatrixAx(nullptr)
    , KLUmatrixAxComplex(nullptr)
    , KLUmatrixIsComplex(false)
{
    klu_l_defaults(&KLUmatrixCommon);   // Initialize KLU (Use klu_l_defaults for long)
}

KLUmatrix::~KLUmatrix()
{
    if (KLUmatrixSymbolic) {
        klu_l_free_symbolic(&KLUmatrixSymbolic, &KLUmatrixCommon);
    }
    if (KLUmatrixNumeric) {
        klu_l_free_numeric(&KLUmatrixNumeric, &KLUmatrixCommon);
    }
}
