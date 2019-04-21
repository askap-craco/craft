///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CodecCODIF.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Coding (and later decoding) of CODIF format.
//
//  Notes:      1) Only encoding of CODIF format is supported at present.
//
//              2) ASKAP data are sampled as 8 channels @32/27 MHz for each
//                 antenna.
//
//  References: [1] VLBI Data Interchange Format (VDIF) Specification
//                  Release 1.1.1, June 2014, Ratified 26-Jun-2009
//
//              [2] CSIRO Oversampled Data Interchange Format (CODIF)
//                  Draft 7 – 24 July 2017, CODIF-Spec-Draft7.pdf
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "CodecCODIF.h"
#include "VCRAFTBitMode.h"

#include <cstring>
#include <cmath>
#include <ctime>
#include <cstdlib>

#undef _VERBOSE
//#define _VERBOSE

///////////////////////////////////////////////////////////////////////////////
// CodecCODIF class implementation and associated helpers.

//////////
// Namespaces defined and/or used in this module.

using std::vector;
using std::deque;
using std::memset;
using std::memcpy;

namespace               // Anonymous namespace for internal helpers.
{
    //////////
    // Constants.

    constexpr int      iMaxPacketSize_c            = 9000;   // Maximum bytes for network transport.
    constexpr int      iBitsPerByte_c              = 8;      // Number of bits in a byte.
    constexpr int      iTimeForIntegerSamples_c    = 27;     // Period in seconds, for a whole number of samples.
    constexpr double   dTolerance_c                = 0.001;  // Tolerance for floating point calculations.
    constexpr int     DUTC                         = 37;

}                       // End anonymous namespace.

//////////
// Part of the main NCodec namespace.

namespace NCodec        // Part of the Codec namespace.
{

///////////////////////////////////////////////////////////////////////////////
// CCodecCODIF class implementation.

    CCodecCODIF::CCodecCODIF( CFileDescriptor *pFile, ECodecMode eCodecMode )
                :ICodec( pFile, eCodecMode )
    {
        m_bFirstInputBlock     = true;
        m_bSynced              = false;
        m_bInitialiseFrameTime = true;
        m_iSampleBlockSize     = 0;
        m_iDataArraySize       = 0;
        m_iDataArrayWords      = 0;
        m_iDataFrameSize       = 0;
        m_iNumberOfSyncRetries = 0;

        m_DFH.Reset();
        m_DataFrameBuffer.clear();

        m_ullBAT0 = 0;
        m_ullFrame0 = 0;
        m_iSkipSamples = 0;
        m_iSampleOffset = 0;
	m_iSamplesPerWord = 0;
	buf  = NULL;
        m_iSamplesPerWord = 0;
	mask = 0;
	convert = false;
	lookup = NULL;
    }

    //////////
    //

    CCodecCODIF::~CCodecCODIF( void )
    {
        m_ullBAT0 = 0;
        m_ullFrame0 = 0;
        m_iSkipSamples = 0;
        m_iSkipSamples = 0;
        m_DataFrameBuffer.clear();
        m_DFH.Reset();

        m_iNumberOfSyncRetries = 0;
        m_iDataFrameSize       = 0;
        m_iDataArraySize       = 0;
        m_iDataArrayWords      = 0;
        m_iSampleBlockSize     = 0;
        m_bInitialiseFrameTime = false;
        m_bSynced              = false;
        m_bFirstInputBlock     = false;
	mask                   = 0;
	if (lookup==NULL) {
	  delete [] lookup;
	  lookup = NULL;
	}
	//if (buf!=NULL) delete buf;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Decoding methods - not currently supported for CODIF format.

    bool CCodecCODIF::ReadAndValidateHeader( void )
    {
        // Not supported by the format so nothing to do.

        m_bHeaderValid = true;

        return true;
    }

    //////////
    //

    bool CCodecCODIF::DecodeHeader( void )
    {
        // Not currently supported for CODIF format hence nothing to do.
        // We therefore return true indicating a non-error condition.

        return true;
    }

    //////////
    //

    bool CCodecCODIF::DecodeChannelData( void )
    {
        // Not currently supported for CODIF format hence nothing to do.
        // We therefore return true indicating a non-error condition.

        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Encoder methods.

    bool CCodecCODIF::ConfigureEncoder( ICodec &rDecoder )
    {
        // Copy the key parameters across from the Decoder to the Encoder
        // (part of the common ICodec base-class) and reset the initial block
        // flag ahead of conversion.

        *static_cast<CCodecParameters *>(this) = *static_cast<CCodecParameters *>(&rDecoder);

        m_bFirstInputBlock = true;

        // Keep a pointer to the decoder for later use during encoding.

        m_pDecoder = &rDecoder;

        return true;
    }

    //////////
    //

    bool CCodecCODIF::EncodeHeader( void )
    {
        // Not relevant for CODIF format as CODIF format does not contain a file
        // header. Rather, Data Frame Headers (DFHs) are handled during channel
        // data conversion.

        return true;
    }

    //////////
    //

    void CCodecCODIF::DumpHeader( void )
    {
        // Meaningless for CODIF format at a file level.
    }

    //////////
    //

    bool CCodecCODIF::operator ()( void )
    {
        // Not used hence nothing to do.

        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Decoding methods.

    bool CCodecCODIF::Initialise( void )
    {
        bool bSuccess = true;

        try
        {
            // Check that the output file is open and ready for writing.

            if ( ! m_pFile->IsOpen() )
            {
                throw string{ "Output file is not open" };
            }

            // Set the time base

            // Need to determine epoch zero for correlator frames and
            // set this to epoch zero for CODIF heaaders.  When
            // processing multiples files, this logic will need to be
            // figured out for the first file then propogated to later
            // files. If processing files independently, this epoch
            // zero will need to be "stored" somewhere. Epoch zero
            // *may* be different for different antenna. Note that the
            // following essentially is a guess at the correlator
            // reset time but this approach should allow stable time
            // determination. It is essentially that data which will
            // be processed at the same time use the same "epoch zero"
            // otherwise there will be delay shifts. If the same Epoch
            // zero is used between different events (but the system
            // has not been reset then there should be no phase shift
            // between data sets. If the system has been reset this
            // approach will introduce a phase shift, even if the
            // underlying sampling clock was not reset, or was
            // syncronysly reset (unless you happen to be very luck
            // and choose the right sample to resyncronise on.

	    // The approach is to compare the BAT and frame ID of the
	    // captured data. Knowing the data sampling rate we can
	    // calculate the BAT of frame ID zero.  We then need to
	    // work out which frame(sample) starts on the next second boundary
	    // There *will* be an up to  +-1/2 sample time rounding with this approach

            unsigned long long startBAT = m_ullTriggerWriteBAT - (m_ullTriggerFrameId * (27.0/32.0));
            m_ullBAT0 = ((startBAT +5e5)/ 1e6); // Round to full second
            m_ullBAT0 *= 1e6;  // Need to do in two lines as compiler is too clever it seems (optimises it away)
	    //m_ullFrame0 = (m_ullBAT0-startBAT)*(32.0/27.0); // Number of frames(samples) from startBAT till first 1sec boundary

	    printf("DEBUG: startBAT=0x%llX\n", startBAT);
	    printf("DEBUG: BAT0=0x%llX\n", m_ullBAT0);
	      
            if ( ! ConfigureDFH() )
            {
                throw string{ "Error forming the Data Frame Header (DFH)" };
            }
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CCodecCODIF::Initialise(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }

        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::Initialise(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    bool CCodecCODIF::EncodeAndWriteChannelData( void )
    {
        bool bSuccess = false;      // Assume failure for now.

        try
        {
            // Check that the output file is open and ready for writing.

            if ( ! m_pFile->IsOpen() )
            {
                throw string{ "Output file is not open" };
            }

            // For the first block, form and output the data frames.

            if ( m_bFirstInputBlock )
            {
                // Reserve storage for a data frame.

                if ( m_DataFrameBuffer.capacity() < static_cast<size_t>( m_iDataFrameSize ) )
                {
                    m_DataFrameBuffer.reserve( m_iDataFrameSize );
                    assert( m_DataFrameBuffer.capacity() >= static_cast<size_t>( m_iDataFrameSize ) );
                }

                m_bFirstInputBlock = false;
            }

            // Process the input data into CODIF-encoded frames and write to file.

            if ( ! HandleCODIFFrameData() )
            {
                throw string{ "Error processing the initial data frame" };
            }

            bSuccess = true;
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CCodecCODIF::EncodeAndWriteChannelData(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::EncodeAndWriteChannelData(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    bool CCodecCODIF::ReadNextBlock( void )
    {
        // Decoding not currently supported for CODIF hence nothing to do.

        return true;
    }

    //////////
    //

    bool CCodecCODIF::Flush( void )
    {
        // The encoder is about to be destroyed so flush any unwritten queued
        // data to disk.

#ifdef _VERBOSE
        fprintf( stderr, "In Flush()\n" );
#endif

        return WriteDataFrames( true );
    }

    //////////
    //

    int CCodecCODIF::DataArraySize( void )
    {
        return m_iDataArraySize;
    }

    //////////
    //

    int CCodecCODIF::SkipBytes( bool *preload )
    {
      int bytestoskip = m_iSkipSamples * m_iSampleBlockSize;
      int vcraftBlock = 4*m_iNumberOfChannels*m_iNumberofPol;
      *preload = (bytestoskip % vcraftBlock) != 0;
      bytestoskip /= vcraftBlock;
      bytestoskip *= vcraftBlock;
      return bytestoskip;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Private internal helpers.

    // Return expected number of samples in a voltage dump - needed to correct time because VCRAFT file contains time
    // of LAST sample. Logic taken from vcraft.py
    unsigned long long voltage_samples(int mode) {
      int samples_per_word[] = {1,2,4,16,1,2,4,16};
      int mode_beams[] = {72,72,72,72,2,2,2,2};

      if (mode>7 || mode < 0) {
	fprintf(stderr, "Error: Mode must be 0..7, not %d\n", mode);
	exit(1);
      }
      int max_sets = 29172;

      //  Number of sample vs mode
      return (max_sets*32*samples_per_word[mode]*72/mode_beams[mode]);
    }

  
    bool CCodecCODIF::ConfigureDFH( void )
    {
        bool bSuccess = false;      // Assume failure for now.

        // Miscellaneous VCRAFT identifiers not presently used:
        // m_iBeamId     - The beam ID requested for the download *0-71.
        // m_iFPGAId     - The FPGA from which this data was downloaded.
        // m_iCardNumber - The card number from which this data was downloaded.

        try
        {
            const uint16_t ui16ThreadId     = 0;     // Not used.
            const uint16_t ui16GroupId      = 0;     // Not used.
            const int      iIsComplex       = 1;     // 1 = complex.

            uint16_t       ui16Channels     = static_cast<uint16_t>( m_iNumberOfChannels );
            char           caStationId[ 3 ] = { static_cast<char>( m_iAntennaNumber ), 0, 0 };

            // Following an initial recommendation from Chris Philips to choose a
            // Data Frame close to (but less than) 9000 bytes for convenience of decoding,
            // it was decided, in discussion with Adam Deller, to reduce this further,
            // to minimise any pre-time-syncrhonised samples being potentially dropped
            // due to a partially filled frame.

            m_iSampleBlockSize = BytesPerTimeSampleForChannels();

            // NB: Period in CODIF is the number of seconds during which there is
            // an integral number of sample periods. Given a sample rate of
            // 32/27 MHz, this means that in 27s there should be 32 x 10^6 samples.

            uint64_t uiSampleIntervalsPerPeriod = llrint( m_dSampleRate * iTimeForIntegerSamples_c );
            assert( fabs( ( uiSampleIntervalsPerPeriod / 1.0e+6 ) - 32.0 ) < dTolerance_c );

	    // Maximum number of samples/frame given the max frame size
	    int samplesPerFrame = iMaxPacketSize_c / m_iSampleBlockSize;

	    // Frame size needs to conform to a number of criteria - keep reducing size until all there are met. Specifically:
	    //    Size multiple of 8 bytes
	    //    Integer number of frames/period
	    //    Original samples packed into 32bit words, so DataArraySize in bytes must equal 4*nchan
	    assert (m_iBitsPerSample <= 16);
	    m_iSamplesPerWord = sizeof(uint32_t)*8/(m_iBitsPerSample*2);
	    
	    while (samplesPerFrame >0 && (samplesPerFrame*m_iSampleBlockSize%8
					  || uiSampleIntervalsPerPeriod%samplesPerFrame
					  || samplesPerFrame % (16/m_iBitsPerSample))) {
	      samplesPerFrame--;
	    }
	    assert (samplesPerFrame > 0);

            m_iDataArraySize   = samplesPerFrame * m_iSampleBlockSize;
            m_iDataFrameSize   = iDFHSize_c + m_iDataArraySize;
	    m_iDataArrayWords  = m_iDataArraySize/sizeof(uint32_t);
            assert( m_iDataFrameSize <= iMaxPacketSize_c );

	    // Get a pointer to the underlying DFH structure.

            CODIFDFH_t *pDFH = static_cast<CODIFDFH_t *>( m_DFH );

            // Create the header anew prior to encoding.

	    if (convert) {
	      if (m_iBitsPerSample==4) {
		struct {signed int x:4;} s;
		if (lookup!=NULL) delete[] lookup;
		lookup = new uint8_t [256];
		for (int i=0; i<256; i++) {
		  int a = s.x = (i & 0xF);      // 4 bits lowest bits - bit extend
		  int b = s.x = ((i>>4) & 0xF); // 4 bits higest bits - bit extend
		  if (a==-8) {
		    a = 7;
		  } else {
		    a = -a;
		  }
		  if (b==-8) {
		    b = 7;
		  } else {
		    b = -b;
		  }
		  lookup[i] = (a&0xF) | ((b<<4)&0xF0);
		}
	      }
	    }

            m_DFH.Reset();

#ifdef _VERBOSE

            printf("createCODIFHeader: %d %d %d %d %d %d %d %d %ld %d %s\n",
                        m_iDataArraySize, ui16ThreadId,
                        ui16GroupId, m_iBitsPerSample,
                        ui16Channels, m_iNumberofPol, m_iSampleBlockSize,
                        iTimeForIntegerSamples_c,
                        uiSampleIntervalsPerPeriod,
                        iIsComplex, caStationId);

#endif

            if ( createCODIFHeader( pDFH, m_iDataArraySize, ui16ThreadId,
                                        ui16GroupId, m_iBitsPerSample,
                                            ui16Channels*m_iNumberofPol, m_iSampleBlockSize,
                                                iTimeForIntegerSamples_c,
                                                    uiSampleIntervalsPerPeriod,
                                                        iIsComplex, caStationId )
                                  != CODIF_NOERROR )
            {
                throw string { "createCODIFHeader() failed" };
            }

            // Set the epoch based on our BAT0

            if ( setCODIFEpochMJD( pDFH, (m_ullBAT0/1e6-DUTC)/(24*60*60)) != CODIF_NOERROR )
            {
                throw string { "setCODIFEpochMJD() failed" };
            }

            // Two's complement
            setCODIFRepresentation(pDFH, 1);

            // Clear the data frame counter and prepare to count data frames in a
            // CODIF Period.

            m_DFH.SetMaxDataFrameNumber( uiSampleIntervalsPerPeriod / samplesPerFrame );

            // Figure out the "previous" Period restart

	    unsigned long long bufferSamples = voltage_samples(m_iMode);
	    printf("DEBUG: Assuming %lld samples per voltage dump\n", bufferSamples);

	    unsigned long long startFrameId = m_ullTriggerFrameId - bufferSamples + m_iSamplesPerWord;
	    // StopFrameId is the time of the last word, so need to allow for the number of samples/32bit word
	    
            unsigned long long EpochMJDSec = getCODIFEpochMJD(pDFH) * 24*60*60;
            unsigned long long BAT0MJDSec = m_ullBAT0/1e6;
            unsigned long long PeriodsSinceBAT0 = startFrameId / uiSampleIntervalsPerPeriod; // Will round down
            int frameseconds = (BAT0MJDSec -DUTC - EpochMJDSec) + PeriodsSinceBAT0 * iTimeForIntegerSamples_c;
            unsigned long framenumber = (startFrameId % uiSampleIntervalsPerPeriod) / samplesPerFrame;

	    mask = (1<<(m_iBitsPerSample*2))-1;

            // Finally set the  seconds and frame number
            setCODIFFrameEpochSecOffset(pDFH, frameseconds);
            setCODIFFrameNumber(pDFH, framenumber);

            // There will be samples that don't fit in a frame (ie first sample may start between frames)
            // These will need to be calculated and eventually discarded

            int initialSamples = startFrameId % samplesPerFrame;

            if ( initialSamples!=0 )
            {
              m_iSkipSamples = samplesPerFrame
		- initialSamples;
              printf("Warning: Will skip %d samples for frame alignment\n", m_iSkipSamples);

	      m_iSampleOffset = m_iSkipSamples%m_iSamplesPerWord;
	      m_DFH.NextFrame(); // Allow for the sample skipping which will happen
            }
            else
            {
                m_iSkipSamples = 0;
		m_iSampleOffset = 0;
            }

            bSuccess = true;
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CCodecCODIF::ConfigureDFH(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::ConfigureDFH(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    bool CCodecCODIF::HandleCODIFFrameData( void )
    {
        // Process the input data blocks from the specified file being converted.

        bool bSuccess = false;      // Assume failure for now.

        try
        {
            if ( ! WriteDataFrames() )
            {
              throw string{ "Error writing data frames" };
            }

                bSuccess = true;
        }
        catch ( string Message )
        {
            fprintf( stderr, "CCodecCODIF::HandleCODIFFrameData(), %s.\n", Message.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::HandleCODIFFrameData(), Unspecified exception caught.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }


    //////////
    //

    WordDeque_t & CCodecCODIF::SampleData( void ) const
    {
        assert ( m_pDecoder != nullptr );
        return m_pDecoder->GetSampleData().GetSamples();
    }

    //////////
    //

    // Multiple time samples per output CODIF word
    void CCodecCODIF::decodeVCRAFTBlock(WordDeque_t & rInput, vector<uint32_t>& vcraftData, vector<uint32_t>& codifData,
				       int wordstoUnpack, int samplesPerWord, int *iWordCount) {
      // Grab next set of original samples
      for (int c=0; c< (wordstoUnpack) && ( ! rInput.empty()); c++) {
	vcraftData[c] = rInput.front();
	if (convert) {
	  if (m_iBitsPerSample==1) {
	    vcraftData[c] ^= 0xCCCCCCCC; // Flip every second sample
	  } else if (m_iBitsPerSample==4) {
	    uint32_t x = vcraftData[c];
	    vcraftData[c] = (x&0x00FF00FF) | lookup[(x&0xFF00)>>8]<<8 | lookup[(x&0xFF000000)>>24]<<24;
	  } else {
	    fprintf(stderr, "Error: Cannot convert %d bit data\n", m_iBitsPerSample);
	    exit(1);
	  }
	}
	
	rInput.pop_front();
	(*iWordCount)++;
      }

      if (samplesPerWord>wordstoUnpack) {
	int samplePerOutword = samplesPerWord/wordstoUnpack;
	for (int i=0; i< samplesPerWord/samplePerOutword; i++) {
	  codifData[i] = 0;
	  for (int j=0; j<samplePerOutword; j++) {
	    for (int k=0; k<wordstoUnpack; k++) {
	      codifData[i] |= ((vcraftData[k]>>(j+i*2)*m_iBitsPerSample*2)&mask)<<(k+j*wordstoUnpack)*m_iBitsPerSample*2;
	    }
	  }
	}
      } else {                                           // 4bit
	int wordPerGroup = wordstoUnpack/samplesPerWord; // 2          4
	for (int i=0; i< samplesPerWord; i++) {          // 0..3       0..1
	  for (int j=0; j<wordPerGroup; j++) {           // 0..1       0..3
	    int c = i*wordPerGroup + j;
	    codifData[c] = 0;
	    for (int k=0; k<samplesPerWord; k++) {
	      codifData[c] |= ((vcraftData[k+j*samplesPerWord]>>(i*m_iBitsPerSample*2))&mask)<<(k*m_iBitsPerSample*2);
	    }
	  }
	}
      }
    }

  
    bool CCodecCODIF::WriteDataFrames( bool bForceFlush )
    {
        bool bSuccess = false;      // Assume failure for now.
	int wordstoUnpack = m_iNumberOfChannels*m_iNumberofPol;

	vector<uint32_t> vcraftData(wordstoUnpack); // +1 in case input words dont align with output words
	vector<uint32_t> codifData(wordstoUnpack); 

        try
        {
            WordDeque_t & rInput = SampleData();

            uint32_t *pbyDataFrame = m_DataFrameBuffer.data();
            assert( pbyDataFrame != nullptr );
	    char * frameptr;
	    

            // Process the frames in each period. Ensure we have at least enough
            // samples to form a complete frame unless we are required to perform
            // a forced flush.

            int iWordCount = 0;
            int iWordsToProcess = static_cast<int>( rInput.size() );

	    // CJP - this outer loop seems redundant
            while ( ( iWordsToProcess >= m_iDataArrayWords    ) || // Need to check edge cases for bits < 16
                    ( bForceFlush && ( iWordsToProcess > 0 ) )  )
            {
                // Start of a new CODIF period hence set the time parameters for
                // this frame and reset the frame counter. We also need to set the
                // frame-time reference now that the preamble frames have been handled.

                int iTotalFrames = iWordsToProcess / m_iDataArrayWords; // Need to check edge cases for bits < 16
		if (iTotalFrames==0)
		  return(true); // EOF

                for ( int iFrame = 0; iFrame < iTotalFrames; iFrame++ )
                {
                    // Zero the output data frame, which covers the case where there
                    // are insufficient samples to form an entire frame. CJP - probably better to discard
		  
		  memset( pbyDataFrame, 0x00, m_iDataFrameSize );

                    // Copy over the formatted Data Frame Header.

		  frameptr = (char*)pbyDataFrame;
		  
		  memcpy( frameptr, static_cast<CODIFDFH_t *>(m_DFH), iDFHSize_c );
		  frameptr += iDFHSize_c;

		    // Copy over any left over bytes from the last frame
		    if (m_iSampleOffset) {
		      void *bytesptr;
		      int ncopy = (m_iSamplesPerWord-m_iSampleOffset)*m_iSampleBlockSize;
		      if (iFrame==0) {
			if (buf==NULL) {
			  // First time through - copy in a block
			  decodeVCRAFTBlock(rInput, vcraftData, codifData, wordstoUnpack, m_iSamplesPerWord, &iWordCount);
			  bytesptr = (char*)&codifData[0] + m_iSampleOffset*m_iSampleBlockSize;

			  buf = new char[ncopy];
			} else {
			  bytesptr = buf;
			}
		      } else {
			// Next frame in this batch, codifdata was loaded in last iteration
			bytesptr = (char*)&codifData[0] + m_iSampleOffset*m_iSampleBlockSize;
		      }
		      memcpy(frameptr, bytesptr, ncopy);
		      frameptr += ncopy;
		    }

                    // Copy over the sample data where it exists. Note that this approach
                    // should leave the last partial frame zeroed out.
		    //  Swizzle samples if necessary

		    if (m_iBitsPerSample==16) { // No change
		      for ( int iWord = 0; ( iWord < m_iDataArrayWords ) && ( ! rInput.empty() ); iWord++ )
			{
			  pbyDataFrame[ iDFHWord_c + iWord ] = rInput.front();
			  rInput.pop_front();
			  iWordCount++;
			}
		    } else if (m_iBitsPerSample==1 || m_iBitsPerSample==4  || m_iBitsPerSample==8) {
		      int samplesPerWord = sizeof(uint32_t)*8/(m_iBitsPerSample*2); // 16
		      int nBlock = m_iDataArrayWords/wordstoUnpack;
		      for ( int iBlock = 0; ( iBlock < nBlock ) && ( ! rInput.empty() ); iBlock++ ) {
			// Grab next set of original samples and convert to codif ordering
			decodeVCRAFTBlock(rInput, vcraftData, codifData, wordstoUnpack, samplesPerWord, &iWordCount);

			if (m_iSampleOffset && iBlock==nBlock-1) {
			  // Just copy the start of this block;
			  memcpy(frameptr, &codifData[0], m_iSampleOffset*m_iSampleBlockSize);
			  frameptr += m_iSampleOffset*m_iSampleBlockSize;
			} else {
			  memcpy(frameptr, &codifData[0], m_iSamplesPerWord*m_iSampleBlockSize);
			  frameptr += m_iSamplesPerWord*m_iSampleBlockSize;
			}
		      }
		    } else {
		      fprintf(stderr, "Do not support %d bits\n", m_iBitsPerSample);
		      exit(1);
		    }

                    // Write the data frame to disk

		    if ( ! m_pFile->Write( pbyDataFrame, m_iDataFrameSize/sizeof(uint32_t) ) )
		      {
			throw string{ "Error writing to the encoded file" };
		      }

                    // Increment the frame counter and handle time increment.

                    m_DFH.NextFrame();

                }
                iWordsToProcess -= iWordCount;
            }

	    if (m_iSampleOffset) {  // Copy over any unused samples
	      int ncopy = (m_iSamplesPerWord-m_iSampleOffset)*m_iSampleBlockSize;
	      if (buf==NULL) {
	      }
	      void * bytesptr = (char*)&codifData[0] + m_iSampleOffset*m_iSampleBlockSize;
	      memcpy(buf, bytesptr, ncopy);
	    }
            bSuccess = true;
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CCodecCODIF::WriteDataFrames(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::WriteDataFrames(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    //////////
    //

    int CCodecCODIF::BytesPerTimeSampleForChannels( void ) const
    {
        // Assume Complex samples - I and Q per sample
        assert ( m_iNumberOfChannels * m_iNumberofPol * m_iBitsPerSample * 2 % 8 ==0);
        return ( m_iNumberOfChannels * m_iNumberofPol * m_iBitsPerSample * 2 / iBitsPerByte_c );
    }

    //////////
    //

    bool CCodecCODIF::SetPartialFrameParams( const int &riFrames )
    {
        m_DFH.SetFrameNumber( riFrames );

        if ( ! m_DFH.SetFrameTime( m_dMJDNow - 1.0 ) )
        {
            return false;
        }

        return true;
    }

    //////////
    //

}       // end namespace NCodec
