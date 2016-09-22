//
//  fdmt_test.c
//  fdmt
//
//  Created by Keith Bannister on 19/07/2016.
//  Copyright (c) 2016 Keith Bannister. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "fdmt.h"


int main(int argc, char* argv[])
{
    int nd = 511;
    int nf = 511;
    int nbeams = 1;
    
    float fmax = 1600.;
    float fmin = fmax - (float)nf;
    
    
    fdmt_dtype* din = malloc(sizeof(fdmt_dtype)*nd*nf*nbeams);
    fdmt_dtype* dout = malloc(sizeof(fdmt_dtype)*nd*nf*nbeams);
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
    fclose(fin);
    
    fdmt_t fdmt;
    fdmt_create(&fdmt, fmin, fmax, nf, nd, nbeams);
    
    for(int i = 0; i < 1; i++) {
        fdmt_execute(&fdmt, din, dout);
    }
    
    FILE* fout = fopen(argv[2], "w");
    if (fout == NULL) {
        perror("Could not open output file");
        exit(EXIT_FAILURE);
    }
    fwrite(dout, sizeof(fdmt_dtype), nd*nf, fout);
    fclose(fout);
    printf("Wrote output file %s\n", argv[2]);
    
}
