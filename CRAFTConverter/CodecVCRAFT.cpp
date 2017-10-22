///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CodecVCRAFT.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "CodecVCRAFT.h"

#include <algorithm>
#include <cstring>
#include <cmath>

//////////
// Anonymous namespace for internal helpers.

namespace
{
    constexpr int iNumberOfVCRAFTChannels_c = 8;                    // VCRAFT number of channels.

    constexpr int iInputBlockSize_c         = 2 * 1042 * 1024;      // Make 2MB for now. NB: VCRAFT
                                                                    // data is sampled at 32/27 Msps so
                                                                    // we could use a multiple/divisor of this.
}

///////////////////////////////////////////////////////////////////////////////
// CCodecVCRAFT class implementation.

namespace NCodec
{
    //////////
    // Construction and destruction.

    CCodecVCRAFT::CCodecVCRAFT( CFileDescriptor *pFile, ECodecMode eCodecMode )
                 :ICodec( pFile, eCodecMode )
    {
        if ( eCodecMode == eCodecModeDecode )
        {
            SetDecoder( this );
        }

        // Clear the byte stream initially.

        m_aHeaderStream.fill( 0x00 );

        // Ensure the std::array<> is correctly sized.

        assert( m_aHeaderStream.size() == iVCRAFTFileHeaderSizeInBytes_c );

        m_bBuffersInitialised = false;
	m_iInputBlockSize = iInputBlockSize_c;
    }

    //////////
    //

    CCodecVCRAFT::~CCodecVCRAFT( void )
    {
        m_bBuffersInitialised = false;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Public codec interface methods that must be overridden.

    bool CCodecVCRAFT::ReadAndValidateHeader( void )
    {
        m_bHeaderValid = false;
        m_eErrorCode   = eCodecErrorNone;

        try
        {
            assert( m_pFile != nullptr );

            if ( m_pFile == nullptr )
            {
                throw eCodecErrorInFileInvalid;
            }
            else if ( ! m_pFile->Read( m_aHeaderStream.data(), m_aHeaderStream.size(), 0, SEEK_SET ) )
            {
                throw eCodecErrorFileReadFailed;
            }
            else if ( ! CVCRAFTParser::ValidateHeader( m_aHeaderStream.data(), m_aHeaderStream.size() ) )
            {
                throw eCodecErrorInvalidHeader;
            }
            else
            {
                m_bHeaderValid = true;
                m_eErrorCode = eCodecErrorNone;
            }
        }
        catch ( ECodecError eErrorCode )
        {
            m_bHeaderValid = false;
            m_eErrorCode = eErrorCode;
        }
        catch ( ... )
        {
            m_bHeaderValid = false;
            m_eErrorCode = eCodecErrorUnspecified;
        }

        return ( m_eErrorCode == eCodecErrorNone );
    }

    //////////
    //

    bool CCodecVCRAFT::DecodeHeader( void )
    {
        if ( m_bHeaderValid )
        {
            if ( m_HeaderDecoder.Decode( m_aHeaderStream.data() ) )
            {
                return SetHeaderParameters();
            }
        }

        return false;
    }

    //////////
    //

    bool CCodecVCRAFT::ConfigureEncoder( ICodec &rDecoder )
    {
        // Not currently supported for the VCRAFT format so return
        // true as there is nothing to do.

        UNREFERENCED_PARAMETER( rDecoder );

        return true;
    }

    //////////
    //

    bool CCodecVCRAFT::EncodeHeader( void )
    {
        // Not currently supported for the VCRAFT format so return
        // true as there is nothing to do.

        return true;
    }

    //////////
    //

    bool CCodecVCRAFT::EncodeAndWriteChannelData( void )
    {
        // Not currently supported for the VCRAFT format so return
        // true as there is nothing to do.

        return true;
    }

    //////////
    //

    bool CCodecVCRAFT::ReadNextBlock( void )
    {
        // For a decoder, we read the input voltage data in blocks.
        //
        // NB:  a) Voltage data starts immediately past the file header;
        //         the first time this method is invoked, it should be aligned.
        //
        //      b) This method returns true if there is data to process
        //         otherwise false.

        bool bDataToProcess = false;

        try
        {
            // The first time through, we setup key parameters.

            if ( ! m_bBuffersInitialised )
            {
	      m_SampleData.SetSampleParams( m_iMode, m_iBitsPerSample,
                                              m_iNumberOfChannels );

                m_bBuffersInitialised = true;
            }

            // Read the next data block from the file directly into the codec's deque.

            ByteDeque_t &rDeque = m_SampleData.GetSamples();

            bDataToProcess = ( m_pFile->Read( rDeque, m_iInputBlockSize, 0, SEEK_CUR ) > 0 );
        }
        catch ( ... )
        {
            bDataToProcess = false;
        }

        return bDataToProcess;
    }

    bool CCodecVCRAFT::SeekForward(int skipBytes )
    {
        // Skip skipBytes forward through the file
        return m_pFile->SeekForward(skipBytes);
    }
  
    //////////
    //

    bool CCodecVCRAFT::setBlockSize(int blockSize)
    {
      m_iInputBlockSize = blockSize;
      return true;
    }
  
    void CCodecVCRAFT::DumpHeader( void )
    {
        fprintf( stdout, "VCRAFT File header:\n size: %d, mode = %d, "
                         "antenna = %d, timestamp = %s, sample rate = %f\n",
                            m_iFileHeaderSize, m_iMode, m_iAntennaNumber,
                                m_sUTC.c_str(), m_dSampleRate );

        fprintf( stdout, "beam = %d, bits/sample (complex) = %d, words = %d, "
                         "FPGA id = %d, card = %d\n",
                            m_iBeamId, m_iBitsPerSample, m_iNumberOfWords,
                                m_iFPGAId, m_iCardNumber );

    }

    //////////
    //

    bool CCodecVCRAFT::operator() ( void )
    {
        // Functor for callback if required. Not currently used.

        return true;
    }


    //////////
    // Private methods.

    template <typename tValue>
    void CCodecVCRAFT::RetrieveParameter( char const *pszKey, tValue &rtValue )
    {
        tValue Value = 0;

        if ( m_HeaderDecoder.RetrieveParameter( pszKey, Value ) )
        {
            rtValue = Value;
        }
    }

    //////////
    //

    bool CCodecVCRAFT::SetHeaderParameters( void )
    {
        bool bSuccess = false;      // Assume failure for now.

        try
        {
            assert( m_bHeaderValid );

            if ( ! m_bHeaderValid )
            {
                throw string{ "Invalid header" };
            }

            // Extract key parameters from the header. Data is used later on to extract
            // the voltage data and generate the output data file.

            RetrieveParameter( "HDR_SIZE",              m_iFileHeaderSize );
            RetrieveParameter( "SAMP_RATE",             m_dSampleRate );
            RetrieveParameter( "CRAFT_MODE",            m_iMode );
            RetrieveParameter( "NBITS",                 m_iBitsPerSample );
            RetrieveParameter( "NPOL",                  m_iNumberofPol );
            RetrieveParameter( "BEAMID",                m_iBeamId );
            RetrieveParameter( "FPGA_ID",               m_iFPGAId );
            RetrieveParameter( "CARD_NO",               m_iCardNumber );
            RetrieveParameter( "ANTENNA_NO",            m_iAntennaNumber );
            RetrieveParameter( "NCHANS",                m_iNumberOfChannels );
            RetrieveParameter( "NOW_MJD",               m_dMJDNow );
            RetrieveParameter( "NOW_BAT",               m_ullNowBAT );
            RetrieveParameter( "START_WRITE_FRAMEID",   m_ulStartWriteFrameId );
            RetrieveParameter( "STOP_WRITE_FRAMEID",    m_ulFinishWriteFrameId );
            RetrieveParameter( "TRIGGER_FRAMEID",       m_ulTriggerFrameId );
            RetrieveParameter( "START_WRITE_BAT",       m_ullStartWriteBAT );
            RetrieveParameter( "STOP_WRITE_BAT",        m_ullFinishWriteBAT );
            RetrieveParameter( "TRIGGER_BAT",           m_ullTriggerWriteBAT );

	    m_iBitsPerSample /= 2;

            // UTC string (ISO formatted) verbatim from the VCRAFT header.

            m_sUTC = m_HeaderDecoder.RetrieveParameterString( "UTC_NOW" );

            // Comma separated list of frequencies applying to each of the (8) channels in MHz.
            // For now leave lists below as comma seperated strings.

            m_sFrequencies = m_HeaderDecoder.RetrieveParameterString( "FREQS" );

            // If any, set those parameters known not to be in the VCRAFT header.

            assert( m_iNumberOfChannels == iNumberOfVCRAFTChannels_c );

            m_bBuffersInitialised = false;

            bSuccess = true;
        }
        catch ( string Message )
        {
            fprintf( stderr, "CCodecVCRAFT::SetHeaderParameters(), %s.\n", Message.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CCodecVCRAFT::SetHeaderParameters, General exception caught.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

}       // end namespace NCodec

