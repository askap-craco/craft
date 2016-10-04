//
//  fdmt_test.c
//  fdmt
//
//  Created by Keith Bannister on 19/07/2016.
//  Copyright (c) 2016 Keith Bannister. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "fdmt.h"
#include "cuda_fdmt.h"
#include "CudaTimer.h"

using namespace std;


int main(int argc, char* argv[])
{
    int nd = 336;
    int nf = 336;
    int nbeams = 36*6;
    
    float fmax = 1600.;
    float fmin = fmax - (float)nf;
    
    
    fdmt_dtype* din = (fdmt_dtype*) malloc(sizeof(fdmt_dtype)*nd*nf*nbeams);
    fdmt_dtype* dout = (fdmt_dtype*) malloc(sizeof(fdmt_dtype)*nd*nf*nbeams);
    printf("Starting!\n");
    
    if (argc != 3) {
        printf("Not enough arguments\n");
        exit(EXIT_FAILURE);
    }
    
    FILE* fin = fopen(argv[1], "r");
  
    //for(int i = 0; i < nd*nf; i++) {
    //din[i] = (fdmt_dtype)1;
    //}
    if (fin == NULL) {
        perror("Could not open input file");
        exit(EXIT_FAILURE);
    }
    fread(din, sizeof(fdmt_dtype), nd*nf, fin);

    // copy to all beams
    for(int b = 1; b < nbeams; b++) {
    	int idx = b*nd*nf;
    	memcpy(&din[idx], din, nd*nf*sizeof(fdmt_dtype));
    }
    fclose(fin);
    
    fdmt_t fdmt;
    fdmt_create(&fdmt, fmin, fmax, nf, nd, nbeams);
    
    CudaTimer t;
    t.start();
    for(int i = 0; i < 1; i++) {
        fdmt_execute(&fdmt, din, dout);
    }
    t.stop();
    t.sync_stop();
    cout << "FDMT took " << t << endl;


    
    FILE* fout = fopen(argv[2], "w");
    if (fout == NULL) {
        perror("Could not open output file");
        exit(EXIT_FAILURE);
    }
    fwrite(dout, sizeof(fdmt_dtype), nd*nf, fout);
    fclose(fout);
    printf("Wrote output file %s\n", argv[2]);
    
}
