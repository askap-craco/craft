
// An attempt at an FDMT but with a loop. - To illustrate the problem for David Li.
// Does not compile. Lookup table values are not correct
// But have the right general shape
// based on fdmt_d8_c8_f0.976.h

const int NC=8; // Number of channels
const float FMIN = 0.976; // Frequency of bottom channel (GHz)
const float FMAX = 0.984; // Frequency of bottom channel (GHz)
const int ND = 8; // Number of output DM trials
const int ND_IN = 2; // number of input dm trials
const float DF = 0.001; // channel interval (GHz)
const float BW = 0.008; // Total bandwidth (GHz)

// FIFO declarations

static ap_shift_reg<fdmt_t, 2> I0D0C1_fifo;
static ap_shift_reg<fdmt_t, 2> I0D1C1_fifo;
static ap_shift_reg<fdmt_t, 2> I0D0C3_fifo;
static ap_shift_reg<fdmt_t, 2> I0D1C3_fifo;
static ap_shift_reg<fdmt_t, 2> I0D0C5_fifo;
static ap_shift_reg<fdmt_t, 2> I0D1C5_fifo;
static ap_shift_reg<fdmt_t, 2> I0D0C7_fifo;
static ap_shift_reg<fdmt_t, 2> I0D1C7_fifo;
static ap_shift_reg<fdmt_t, 2> I1D0C1_fifo;
static ap_shift_reg<fdmt_t, 3> I1D1C1_fifo;
static ap_shift_reg<fdmt_t, 3> I1D2C1_fifo;
static ap_shift_reg<fdmt_t, 2> I1D0C3_fifo;
static ap_shift_reg<fdmt_t, 3> I1D1C3_fifo;
static ap_shift_reg<fdmt_t, 3> I1D2C3_fifo;
static ap_shift_reg<fdmt_t, 2> I2D0C1_fifo;
static ap_shift_reg<fdmt_t, 3> I2D1C1_fifo;
static ap_shift_reg<fdmt_t, 4> I2D2C1_fifo;
static ap_shift_reg<fdmt_t, 5> I2D3C1_fifo;

const int NITER = 3;

// LOOKUP TABLES
// Note: this is the iteration setup
// shapes are [Nchan, Ndm, Nt] - Nt is not important here.
// Size = Nchan*Ndm = size of temporary variable to pass between iteration
// Should be approximately constant.
// Iteration 1 in=[  8   2 258] size=16 out=[  4   3 259] size=12 
// Iteration 2 in=[  4   3 259] size=12 out=[  2   5 261] size=10 
// Iteration 3 in=[  2   5 261] size=10 out=[  1   8 264] size=8

const int MAX_SIZE = 16; // Maximum size of an iteration buffer - a bit wasteful

// Number of output channels for each iteration
const int NCHAN_BY_ITER = {4, 2, 1};

// Number of DMs to do for each channel, for each iteration
// Assume its the same for all output channels in the iteration
// In practice it might be different for each channel
const int NDM_BY_ITER= {3, 5, 8};

// Lookup table for previous index - I know this doesn't compile
// TODO: Flatten this into a 1D array.
const int[NITER][MAX_SIZE] PREVIDX {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // Iteration 0 is empty - 16 values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0}, // Iteration 1 - 12 valid values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0},  // Iteration 2 - 10 valid values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0} // Iteration 3 - 8 valid values
}

// Lookup table for offsets. Not complete or correct
const int[NITER][MAX_SIZE] OFFSETS {
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // Iteration 0 is empty - 16 values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0}, // Iteration 1 - 12 valid values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0},  // Iteration 2 - 10 valid values
  {0,1,1,1,2,2,2,2,3,4,5,6,0,0,0,0} // Iteration 3 - 8 valid values
}

// Lookup table for FIFOS.
// ******* THIS IS THE KEY THING THAT I DON'T KNOW HOW TO DO!
// What is the pointer type?? THe FIFOS are all diffent types
// Because the template contains the FIFO size as a template paramter
// I don't know how to address the FIFOS separately
//
// Maybe I need my own shiftreg class but make it derive from an abstract
// base class. But does HLS support this?
// https://stackoverflow.com/questions/6342464/template-class-pointer-c-declaration
static fdmt_fifo<fdmt_t, ???????> * FIFOS[NITER][MAX_SIZE] {
  {&I0D0C1_fifo, &I0D0C1_fifo, &I0D1C1_fifo, ... (etc, etc)} // Iteration 1
  {&I1D0C1_fifo, &I1D0C2_fifo, &I1D1C2_fifo, ... (etc, etc, etc)} // Iteration 2
  {&I2D0C1_fifo, &I2D0C2_fifo, &I2D1C2_fifo, ... (etc, etc, etc)} // Iteration 2

}


// Make it easier and only have 1D inputs, rather than 2D inputs
// As previously discussed
void fdmt_process_asforloop(const fdmt_t in[ND_IN*NC], fdmt_t out[ND])
{
  // Iteration buffer - could be implemented as a pingpong buffer
  // if we do the iterations sequentailly - but not if we're doing
  // a dataflow
  // Makking it MAX_SIZE is a bit wastful as the size changes with iteration
  fdmt_t buf[NITER][MAX_SIZE]; 

    // Initialise iteration buffer - could add a special case to the loop - but that's an optimisation
  for (int i = 0; i < ND_IN*NC; i++) {
    buf[0][i] = in[i];
  }
  
  for(int iter = 0; iter < NITER; iter++) {
    int nchan = NCHAN_BY_ITER[iter];
    // increment iteration pointers

    for (int ochan = 0; ochan < nchan; ochan++) {
      for (int idm = 0; idm < NDM_BY_CHAN[iter]; idm++) {
	int idx = ochan*nchan + idm;

	// Here is the caculation.
	// tyically it looks like this:
	// fdmt_t I1D2C0 = I0D1C0 + I0D1C1_fifo.read(1);
	// We need LOOKUP tables for 3 things:
	// 1. Which D? and C? to get in the first part of the sum (PREVIDX)
	// 2. Which D? and C? to get the fifo for the second part of the sum (FIFOS)
	// 3. What offset to read from the fifo (FIFOS)
	buf[iter][idx] = buf[iter-1][PREVIDX[iter][idx]] + FIFOS[iter][idx].read(OFFSETS[iter][idx])
      }
    }
  }

  // Write output from buffer. Again could be optimised in final iteration of loop
  for (int i = 0; i < ND; i++) {
    out[i] = buf[NITER-1][i];
  }
    
}
