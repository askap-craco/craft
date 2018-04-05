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
  //    constexpr int      iSamplesPerFrame_c          = 256;    // Samples per frame (really should be calculated)   //CJP TODO CALCULATE
    constexpr int      iBitsPerByte_c              = 8;      // Number of bits in a byte.
    constexpr int      iTimeForIntegerSamples_c    = 27;     // Period in seconds, for a whole number of samples.
    constexpr double   dTolerance_c                = 0.001;  // Tolerance for floating point calculations.

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
    }

    //////////
    //

    CCodecCODIF::~CCodecCODIF( void )
    {
        m_ullBAT0 = 0;
        m_ullFrame0 = 0;
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
            m_ullBAT0 = ((startBAT +5e5)/ 1e6); // Round up to full second
            m_ullBAT0 *= 1e6;  // Need to do in two lines as compiler is too clever it seems (optimises it away)
	    m_ullFrame0 = (m_ullBAT0-startBAT)*(32.0/27.0); // Number of frames(samples) from startBAT till first 1sec boundary

	    printf("DEBUG: startBAT=0x%llX\n", startBAT);
	    printf("DEBUG: BAT0=0x%llX\n", m_ullBAT0);
	    printf("DEBUG: Frame0=%llu\n", m_ullFrame0);
	      
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

    int CCodecCODIF::SkipBytes( void )
    {
      printf("CCodecCODIF::SkipBytes: Skipping %d blocks of %d bytes\n", m_iSkipSamples, m_iSampleBlockSize);
      int bytestoskip = m_iSkipSamples * m_iSampleBlockSize;
      printf("Skipping %d->", bytestoskip);
      int vcraftBlock = 4*m_iNumberOfChannels*m_iNumberofPol;
      bytestoskip /= vcraftBlock;
      bytestoskip *= vcraftBlock;
      printf("%d\n", bytestoskip);
      
      return bytestoskip;
    }

  // Update the number of samples which need to be skipped, based on a passed number of bytes
  // Needed because raw data from telescopes groups data into blocks of 32bits (per coarse channel)
    void CCodecCODIF::updateSkip(int skipBytes )
    {
      printf("CCodecCODIF::updateSkip: %d\n", skipBytes);
      m_iSkipSamples -= skipBytes ;
    }

  
    ///////////////////////////////////////////////////////////////////////////
    // Private internal helpers.

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

	    // CJP Need to set based on # bits

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

	    printf("********samplesPerFrame=%d  DataArray=%d\n", samplesPerFrame, m_iDataFrameSize);

	    // Get a pointer to the underlying DFH structure.

            CODIFDFH_t *pDFH = static_cast<CODIFDFH_t *>( m_DFH );

            // Create the header anew prior to encoding.

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

            if ( setCODIFEpochMJD( pDFH, m_ullBAT0/(24*60*60*1e6)) != CODIF_NOERROR )
            {
                throw string { "setCODIFEpochMJD() failed" };
            }

            // Two's complement
            setCODIFRepresentation(pDFH, 1);

            // Clear the data frame counter and prepare to count data frames in a
            // CODIF Period.

            m_DFH.SetMaxDataFrameNumber( uiSampleIntervalsPerPeriod / samplesPerFrame );

            // Figure out the "previous" Period restart

            unsigned long long EpochMJDSec = getCODIFEpochMJD(pDFH) * 24*60*60;
            unsigned long long BAT0MJDSec = m_ullBAT0/1e6;
            unsigned long long PeriodsSinceBAT0 = (m_ullTriggerFrameId-m_ullFrame0) / uiSampleIntervalsPerPeriod; // Will round down
            int frameseconds = (BAT0MJDSec - EpochMJDSec) + PeriodsSinceBAT0 * iTimeForIntegerSamples_c;
            unsigned long framenumber = ((m_ullTriggerFrameId-m_ullFrame0) % uiSampleIntervalsPerPeriod) / samplesPerFrame;

            // Finally set the  seconds and frame number
            setCODIFFrameEpochSecOffset(pDFH, frameseconds);
            setCODIFFrameNumber(pDFH, framenumber);

            // There will be samples that don't fit in a frame (ie first sample may start between frames)
            // These will need to be calculated and eventually discarded

            int initialSamples = (m_ullTriggerFrameId-m_ullFrame0) % samplesPerFrame;

            if ( initialSamples!=0 )
            {
              m_iSkipSamples = samplesPerFrame
		- initialSamples;
              printf("Warning: Will skip %d samples for frame alignment\n", m_iSkipSamples);
	      m_DFH.NextFrame(); // Allow for the sample skipping which will happen
            }
            else
            {
                m_iSkipSamples = 0;
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
            // Attempt to sync to the first whole second. If we are unsuccessful,
            // then wait for subsequent input data blocks until we have sufficient samples
            // before flagging an error.

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

    bool CCodecCODIF::WriteDataFrames( bool bForceFlush )
    {
        bool bSuccess = false;      // Assume failure for now.

	vector<uint32_t> vcraftData(m_iNumberOfChannels*m_iNumberofPol);
	vector<uint32_t> codifData(m_iNumberOfChannels*m_iNumberofPol);

        try
        {
            WordDeque_t & rInput = SampleData();

            uint32_t *pbyDataFrame = m_DataFrameBuffer.data();
            assert( pbyDataFrame != nullptr );

            // Process the frames in each period. Ensure we have at least enough
            // samples to form a complete frame unless we are required to perform
            // a forced flush.

            int iWordCount = 0;
            int iWordsToProcess = static_cast<int>( rInput.size() );

#ifdef _VERBOSE
            fprintf( stderr, "WriteDataFrames(), %d words to process\n",
                              iWordsToProcess );
#endif

            while ( ( iWordsToProcess >= m_iDataArrayWords    ) ||
                    ( bForceFlush && ( iWordsToProcess > 0 ) )  )
            {
                // Start of a new CODIF period hence set the time parameters for
                // this frame and reset the frame counter. We also need to set the
                // frame-time reference now that the preamble frames have been handled.

                int iTotalFrames = iWordsToProcess / m_iDataArrayWords;

                if ( ( iWordsToProcess % m_iDataArrayWords ) != 0 )
                {
                    iTotalFrames++;
                }

                for ( int iFrame = 0; iFrame < iTotalFrames; iFrame++ )
                {
                    // Zero the output data frame, which covers the case where there
                    // are insufficient samples to form an entire frame.

                    memset( pbyDataFrame, 0x00, m_iDataFrameSize );

                    // Copy over the formatted Data Frame Header.

                    memcpy( pbyDataFrame, static_cast<CODIFDFH_t *>(m_DFH), iDFHSize_c );

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
		    } else if (m_iBitsPerSample==1) {
		      printf("DEBUG: Processing %d words\n", m_iDataArrayWords);
		      
		      int mask = (1<<(m_iBitsPerSample*2))-1;
		      assert(mask==0x3); // Temp check - remove

		      printf("DEBUG:  nbit=%d nchan= %d  nPol=%d\n", m_iBitsPerSample, m_iNumberOfChannels, m_iNumberofPol);
		      int samplesPerWord = sizeof(uint32_t)*8/(m_iBitsPerSample*2); // 16
		      printf("DEBUG: samplesPerWord = %d\n", samplesPerWord);
		      assert(samplesPerWord % m_iNumberOfChannels*m_iNumberofPol == 0); // Must have an integer number of channels per word
		      int samplePerOutword = samplesPerWord/  m_iNumberOfChannels*m_iNumberofPol;
		      printf("DEBUG: samplesPerOutWord = %d\n", samplePerOutword);
		      
		      for ( int iWord = 0; ( iWord < m_iDataArrayWords/(m_iNumberOfChannels*m_iNumberofPol) ) && ( ! rInput.empty() ); iWord++ ) {
			// Grab next set of original samples
			for (int c=0; c< (m_iNumberOfChannels*m_iNumberofPol) && ( ! rInput.empty()); c++) {
			  vcraftData[c] = rInput.front();
			  rInput.pop_front();
			  iWordCount++;
			}

			for (int i=0; i< samplesPerWord/samplePerOutword; i++) {
			  codifData[i] = 0;
			  for (int j=0; j<samplePerOutword; j++) {
			    for (int k=0; k<m_iNumberOfChannels*m_iNumberofPol; k++) {
			      codifData[i] |= ((vcraftData[k]>>(j+i*2)*m_iBitsPerSample*2)&mask)<<(k+j*m_iNumberOfChannels*m_iNumberofPol)*m_iBitsPerSample*2;
			    }
			  }
			}

			for (int i=0; i<m_iNumberOfChannels*m_iNumberofPol; i++) {
			  pbyDataFrame[ iDFHWord_c + iWord*(m_iNumberOfChannels*m_iNumberofPol) + i ] = codifData[i];
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
