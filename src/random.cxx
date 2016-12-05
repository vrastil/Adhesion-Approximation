
#include "CBRNG_Random.h"

    //
    // Get keys for each slab in the x axis that this rank contains 
    //

    tnx1 = 0;
    tny1 = 0;
    tnz1 = 0;

    std::vector<unsigned long> slab_keys;
    slab_keys.resize(tngx);
    GetSlabKeys(&slab_keys[0], tnx1, tngx, seed);

    //
    // Generate Gaussian white noise field 
    //
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(long i=0; i<tngx; ++i) 
    {
      unsigned long ikey = slab_keys[i];
      for(long j=0; j<tngy; ++j) 
      {
        long jg = tny1 + j;
        for(long k=0; k<tngz; ++k) 
        {
          long kg = tnz1 + k;
          unsigned long index = jg*tngy + kg;
          double rn1, rn2;

          GetRandomDoublesWhiteNoise(rn1, rn2, ikey, index);
          index = (i*tngy+j)*tngz+k;
          double rn = rn1*rn1 + rn2*rn2;
          rho[index].re = rn2 * sqrt(-2.0*log(rn)/rn);
          rho[index].im = 0.0;
        }
      }
    }

