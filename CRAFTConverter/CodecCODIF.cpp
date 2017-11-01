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

namespace               // Anonymous namespace for internal helpers.
{
    //////////
    // Constants.

    constexpr int      iMaxPacketSize_c            = 9000;   // Maximum bytes for network transport.
    constexpr int      iSamplesPerFrame_c          = 256;    // Samples per frame (really should be calculated)
    constexpr int      iMaxSyncTries_c             = 3;      // Maximum number of tries before sync deemed to have failed.
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
        m_iDataFrameSize       = 0;
        m_iNumberOfSyncRetries = 0;

        m_DFH.Reset();
        m_DataFrameBuffer.clear();
	ull_BAT0 = 0;
	m_iSkipSamples = 0;
    }

    //////////
    //

    CCodecCODIF::~CCodecCODIF( void )
    {
        m_DataFrameBuffer.clear();
        m_DFH.Reset();

        m_iNumberOfSyncRetries = 0;
        m_iDataFrameSize       = 0;
        m_iDataArraySize       = 0;
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
	    // *will* be different for different antenna Note that the
	    // following essentially is a guess at the correlator
	    // reset time but this approach Should allow stable time
	    // determination. It is essentially that data which will
	    // be processed at the same time use the same "epoch zero"
	    // otherwise there will be delay shifts

	    unsigned long long startBAT = m_ullStartWriteBAT - (m_ulStartWriteFrameId * (27.0/32.0));
	    ull_BAT0 = ((startBAT +5e5)/ 1e6); // Round to full second
	    ull_BAT0 *= 1e6;  // Need tp do in two lines and compiler is too clever it seems
	    long long BATerr = startBAT - ull_BAT0;
	    if (BATerr>5e5) BATerr -= 1e6;
	    printf("Note: BAT offset for implied correlator reset is %.2f msec\n", BATerr/1.0e3);

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

    int CCodecCODIF::DataArraySize()
    {
      return m_iDataArraySize;
    }

    int CCodecCODIF::skipBytes()
    {
      return m_iSkipSamples * m_iSampleBlockSize;
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

            m_iSampleBlockSize = BytesPerTimeSampleForChannels();
            m_iDataArraySize   = iSamplesPerFrame_c  * m_iSampleBlockSize;
            m_iDataFrameSize   = iDFHSize_c + m_iDataArraySize;
            assert( m_iDataFrameSize <= iMaxPacketSize_c );

            // Initialise the remaining DFH parameters. Remember:
            // Data Frame = Data Frame Header + Data Array.

            // NB: Period in CODIF is the number of seconds during which there is
            // an integral number of sample periods. Given a sample rate of
            // 32/27 MHz, this means that in 27s there should be 32 x 10^6 samples.

            uint64_t uiSampleIntervalsPerPeriod =
                static_cast<uint64_t>( floor( m_dSampleRate * iTimeForIntegerSamples_c ) );
            assert( fabs( ( uiSampleIntervalsPerPeriod / 1.0e+6 ) - 32.0 ) < dTolerance_c );

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
            if ( createCODIFHeader(pDFH, m_iDataArraySize, ui16ThreadId,
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
            if ( setCODIFEpochMJD(pDFH, ull_BAT0/(24*60*60*1e6)) != CODIF_NOERROR )
            {
                throw string { "setCODIFEpochMJD() failed" };
            }

            // Two's complement

            setCODIFRepresentation(pDFH, 1);

            // Clear the data frame counter and prepare to count data frames in a
            // CODIF Period.

            m_DFH.SetMaxDataFrameNumber( uiSampleIntervalsPerPeriod / iSamplesPerFrame_c );
	    
	    // Figure out the "previous" Period restart
	    unsigned long long EpochMJDSec = getCODIFEpochMJD(pDFH) * 24*60*60;
	    unsigned long long BAT0MJDSec = ull_BAT0/1e6;
	    unsigned long long PeriodsSinceBAT0 = m_ulStartWriteFrameId / uiSampleIntervalsPerPeriod; // Will round down
	    int frameseconds = (BAT0MJDSec - EpochMJDSec) + PeriodsSinceBAT0 * iTimeForIntegerSamples_c;
	    unsigned long framenumber = (m_ulStartWriteFrameId % uiSampleIntervalsPerPeriod) / iSamplesPerFrame_c;
	      
	    // Finally set the  seconds and frame number
	    setCODIFFrameEpochSecOffset(pDFH, frameseconds);
	    setCODIFFrameNumber(pDFH, framenumber);

	    // There will be samples that don't fit in a frame (ie first sample may start between frames)
	    // These will need to be calculated and eventually discarded
	    int initialSamples = m_ulStartWriteFrameId % iSamplesPerFrame_c;
	    if (initialSamples!=0) {
	      m_iSkipSamples = iSamplesPerFrame_c - initialSamples;
	      printf("Warning: Will skip %d samples for frame alignment\n", m_iSkipSamples);
	    } else {
	      m_iSkipSamples = 0;
	    }
	    
            m_DFH.NextFrame(); // Allow for the sample skipping which will happen

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

    bool CCodecCODIF::SyncAndHandlePreSyncSamples( void )
    {
        // Locate the first whole second and process all preamble samples
        // into a set of frames in a partial-period, if/as required.

        // Upon exit we should have removed leading samples so that the remaining
        // queued data starts at the boundary of the first CODIF Period.

        bool bSuccess = false;  // Assume failure for now.

	
        try
        {
            // If we are already syncrhonised then there is nothing to do.

            if ( ! m_bSynced )
            {
                // Decompose the time into whole and fractional parts.

                double  dMJDWholeSec    = 0.0;
                double  dMJDFracSec     = 0.0;

                dMJDFracSec = std::modf( m_dMJDNow, &dMJDWholeSec );

                assert( m_dSampleRate > 0.0 );
                assert( m_iBitsPerSample > 0 );

                // Compute the number of samples to the next whole seconds.

                const int iTimeSamplesToWholeSec = std::abs( ( 1.0 - dMJDFracSec ) * m_dSampleRate );

                // Handle the preamble samples, if any, then forego further processing
                // here and defer to normal sample-aligned processing.

#ifdef _VERBOSE
                fprintf( stderr, "SyncAndHandlePreSyncSamples(), %d Samples to synchronise \n", iTimeSamplesToWholeSec );
#endif

                if ( ! WritePreambleFrames( iTimeSamplesToWholeSec ) )
                {
                    throw string{ "Can't write preamble samples to data frames" };
                }

                m_bSynced = true;
            }

            bSuccess = true;
        }
        catch ( string Message )
        {
            fprintf( stderr, "CCodecCODIF::SyncAndHandlePreSyncSamples(), %s.\n", Message.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::SyncAndHandlePreSyncSamples(), Unspecified exception caught.\n" );
            bSuccess = false;
        }

        if ( ! bSuccess )
        {
            m_iNumberOfSyncRetries++;
        }

        return bSuccess;
    }

    //////////
    //

    ByteDeque_t & CCodecCODIF::SampleData( void ) const
    {
        assert ( m_pDecoder != nullptr );
        return m_pDecoder->GetSampleData().GetSamples();
    }

    //////////
    //

    bool CCodecCODIF::WriteDataFrames( bool bForceFlush )
    {
        bool bSuccess = false;      // Assume failure for now.


        try
        {
            ByteDeque_t & rInput = SampleData();

            byte_t *pbyDataFrame = m_DataFrameBuffer.data();
            assert( pbyDataFrame != nullptr );

            // Process the frames in each period. Ensure we have at least enough
            // samples to form a complete frame unless we are required to perform
            // a forced flush.

            int iByteCount = 0;
            int iBytesToProcess = static_cast<int>( rInput.size() );

#ifdef _VERBOSE
            fprintf( stderr, "WriteDataFrames(), %d bytes to process\n",
                              iBytesToProcess );
#endif

            while ( ( iBytesToProcess >= m_iDataArraySize    ) ||
                    ( bForceFlush && ( iBytesToProcess > 0 ) )  )
            {
                // Start of a new CODIF period hence set the time parameters for
                // this frame and reset the frame counter. We also need to set the
                // frame-time reference now that the preamble frames have been handled.

                int iTotalFrames = iBytesToProcess / m_iDataArraySize;

                if ( ( iBytesToProcess % m_iDataArraySize ) != 0 )
                {
                    iTotalFrames++;
                }

                for ( int iFrame = 0; iFrame < iTotalFrames; iFrame++ )
                {
                    // Zero the output data frame, which covers the case where there
                    // are insufficient samples to form an entire frame.

                    std::memset( pbyDataFrame, 0x00, m_iDataFrameSize );

                    // Copy over the formatted Data Frame Header.


		    std::memcpy( pbyDataFrame, static_cast<CODIFDFH_t *>(m_DFH), iDFHSize_c );
                    // Copy over the sample data where it exists. Note that this approach
                    // should leave the last partial frame zeroed out.

                    for ( int iByte = 0; ( iByte < m_iDataArraySize ) && ( ! rInput.empty() ); iByte++ )
                    {
                        pbyDataFrame[ iDFHSize_c + iByte ] = rInput.front();
                        rInput.pop_front();
                        iByteCount++;
                    }

                    // Write the data frame to disk.

                    if ( ! m_pFile->Write( pbyDataFrame, m_iDataFrameSize ) )
                    {
                        throw string{ "Error writing to the encoded file" };
                    }

                    // Increment the frame counter and handle time increment.

                    m_DFH.NextFrame();
                }
                iBytesToProcess -= iByteCount;
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

    bool CCodecCODIF::WritePreambleFrames( const int &riPreambleTimeSamples )
    {
        bool bSuccess = false;       // Assume failure for now.

        try
        {
            ByteDeque_t &rInput = SampleData();

            byte_t *pbyDataFrame = m_DataFrameBuffer.data();
            assert( pbyDataFrame != nullptr );

            // Process the preamble frames in this partial period.

            int iByteCount = 0;
            int iBytesToProcess = riPreambleTimeSamples * BytesPerTimeSampleForChannels();

            bool bInitialiseFrameTime = true;

#ifdef _VERBOSE
            fprintf( stderr, "WritePreambleFrames(), %d samples to synchronise (%d bytes), Queued: %ld\n",
                              riPreambleTimeSamples,
                              iBytesToProcess,
                              rInput.size() );
#endif

            while ( rInput.size() >= static_cast<size_t>( m_iDataArraySize ) )
            {
                // Work out how many frames there will be in this partial period
                // and the header parameters to suit.

                int iTotalFrames = iBytesToProcess / m_iDataArraySize;

                if ( ( iBytesToProcess % m_iDataArraySize ) != 0 )
                {
                    iTotalFrames++;
                }

                // Set the preamble-frame's time, number etc.

                if ( bInitialiseFrameTime )
                {
                    SetPartialFrameParams( iTotalFrames );
                    bInitialiseFrameTime = false;
                }

                // Form and output the frames.

                for ( int iFrame = 0; iFrame < iTotalFrames; iFrame++ )
                {
                    // Zero the output data frame, which covers the case where there
                    // are insufficient samples to form an entire frame.

                    std::memset( pbyDataFrame, 0x00, m_iDataFrameSize );

                    // Copy over the formatted Data Frame Header.

                    std::memcpy( pbyDataFrame, static_cast<byte_t *>(m_DFH), iDFHSize_c * sizeof( byte_t ) );

                    // Copy over the sample data where it exists. Note that this approach
                    // should leave the last partial frame zeroed out.

                    for ( int iByte = 0; ( iByte < m_iDataArraySize ) && ( ! rInput.empty() ); iByte++ )
                    {
                        pbyDataFrame[ iDFHSize_c + iByte ] = rInput.front();
                        rInput.pop_front();
                        iByteCount++;
                    }

                    // Write the data frame to disk.

                    if ( ! m_pFile->Write( pbyDataFrame, m_iDataFrameSize ) )
                    {
                        throw string{ "Error writing to the encoded file" };
                    }
                }

                iBytesToProcess -= iByteCount;
            }

            // Finally, should there be any residual unread samples, clear
            // them from the input queue thereby ensuring we are synced to the whole
            // second for subsequent processing.

            const int iBytesToRemove = iBytesToProcess - iByteCount;

            if ( ( iBytesToRemove > 0 ) &&
                 ( iBytesToRemove <= static_cast<int>( rInput.size() ) ) )
            {
                rInput.erase( rInput.begin(), rInput.begin() + iBytesToRemove );
            }

            bSuccess = true;
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CCodecCODIF::WritePreambleFrames(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecCODIF::WritePreambleFrames(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    int CCodecCODIF::BytesPerTimeSampleForChannels( void ) const
    {
      // Assume Complex samples - I and Q per sample
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
