///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CodecVDIF1.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Coding and decoding of VDIF format.
//
//  Notes:      Refer to section 9.4 Multi-channel complex-data Data Array format
//              of [1].
//
//  References: [1] http://vlbi.org/vdif/docs/VDIF_specification_Release_1.1.1.pdf
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "CodecVDIF1.h"

#include <cstring>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Namespaces or parts thereof used and/or defined.

using namespace NCodec;

//////////
// Anonymous namespace for internal helpers.

namespace
{
    //////////
    // Aliases and constants.

    using DFH_t = CCodecVDIF1::VDIF1DFH_t;
  //constexpr size_t iDFHSize_c = sizeof( DFH_t );

    // Helper class to simplify manipulation of the VDIF1 data frame header.

    class CDataFrameHeader
    {
        public:

            //////////
            // Construction, assignment and destruction.

            CDataFrameHeader( DFH_t *pDFH, bool bInitialise = false )
            {
                m_pDFH = pDFH;

                Reset();

                if ( bInitialise )
                {
                    Initialise();
                }
            }

            //////////
            //

            virtual ~CDataFrameHeader( void )
            {
                m_pDFH = nullptr;
            }

            //////////
            // The following are not implemented for this class.

            CDataFrameHeader ( void ) = delete;
            CDataFrameHeader ( const CDataFrameHeader &rRhs ) = delete;
            CDataFrameHeader & operator = ( const CDataFrameHeader &rRhs ) = delete;

        private:

            //////////
            // Private attributes.

            DFH_t    *m_pDFH;

            //////////
            // Private methods.

            void Reset( void )
            {
                assert ( m_pDFH != nullptr );
                memset( m_pDFH, 0x00, sizeof( DFH_t ) );
            }

            //////////
            //

            void Initialise( void )
            {
                assert ( m_pDFH != nullptr );

                // Word 0.
                m_pDFH->uInvalidData            = 0;    // Invalid data (taged by source): valid=0, invalid=1.
                m_pDFH->uLegacyMode             = 0;    // 0 - standard VDIF DFH; 1 - legacy header length mode.
                m_pDFH->uSecondsFromEpoch       = 0;    // Seconds from reference epoch; see ref. [1], Note 2, pp 6.

                // Word 1.
                m_pDFH->uReserved               = 0;    // Unassigned - set to 0.
                m_pDFH->uRefEpoch               = 0;    // Reference Epoch for second count.; see ref. [1], Note 2, pp 6.
                m_pDFH->uDataFrameNumber        = 0;    // Data Frame # starting at zero; integral number of Data Frames per second

                // Word 2.
                m_pDFH->uVDIFVersion            = 1;    // VDIF version number.
                m_pDFH->uLog2Channels           = 0;    // Representation (N) of number of channels s.t., # Channels = 2^N.
                m_pDFH->uDataFrameLength        = 0;    // Includes the DFH and must be a multiple of 8 bytes. Max <= 2^27 bytes.

                // Word 3.
                m_pDFH->uDataType               = 1;    // 0 - Real; 1 - Complex data.
                m_pDFH->uBitsPerSample          = 0;    // # bits in each sample. Complex => 2 x bits (I & Q).
                m_pDFH->uThreadId               = 0;    // Number of each time-series of data DFs from the same set of subbands.
                m_pDFH->uStationId              = 0;    // 2 ASCII ID chars or unsigned numeric. If first 8 bits <48 (0x30) then ASCII assumed.

                // Word 4.
                m_pDFH->uExtendedDataVersion    = 0;    // EDV = 0 if not used.
                m_pDFH->uExtendedUserData1      = 0;    // NB: Unique EDV # assigned via http://www.vlbi.org/vdif/.

                // Words 5 through 7 respectively.
                m_pDFH->uExtendedUserData2      = 0;    // Set to zero if EDV = 0.
                m_pDFH->uExtendedUserData3      = 0;    // ditto
                m_pDFH->uExtendedUserData4      = 0;    // ditto
            }
    };

}       // end anonymous namespace.

///////////////////////////////////////////////////////////////////////////////
// CCodecVDIF1 class implementation; part of the NCodec namespace.

namespace NCodec
{
    //////////
    //

    CCodecVDIF1::CCodecVDIF1( void )
        :ICodec()
    {
        CDataFrameHeader DFH( &m_sDFH );
    }

    //////////
    //

    CCodecVDIF1::CCodecVDIF1( CFileDescriptor *pFile, ECodecMode eCodecMode )
        :ICodec( pFile, eCodecMode )
    {
        CDataFrameHeader DFH( &m_sDFH );
    }

    //////////
    //

    CCodecVDIF1::~CCodecVDIF1( void )
    {
        // Nothing to do.
    }

    //////////
    // Overrides for inherited virtual interface methods.

    bool CCodecVDIF1::ReadAndValidateHeader( void )
    {
        // Not currently supported. Attempting to ustilise this
        // method should result in the header being flagged as
        // invalid until supported.

        m_bHeaderValid = false;
        return false;
    }

    //////////
    //

    bool CCodecVDIF1::DecodeHeader( void )
    {
        // Not currently supported for the VDIF format but return
        // true as there is nothing to do.

        return true;
    }

    //////////
    //

    bool CCodecVDIF1::DecodeChannelData( void )
    {
        // Not currently supported for the VDIF format but return
        // true as there is nothing to do.

        return true;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Encoder methods.

    bool CCodecVDIF1::ConfigureEncoder( ICodec &rDecoder )
    {
        // Copy the key parameters across from the Decoder to the Encoder...

        *static_cast<CCodecParameters *>( this ) = *static_cast<CCodecParameters *>( &rDecoder );

        // and keep a pointer to the decoder for later use during encoding.

        m_pDecoder = &rDecoder;

        return true;
    }

    //////////
    //

    bool CCodecVDIF1::EncodeHeader( void )
    {
        // Not used. VDIF format does not contain a file header, per se. Rather
        // it has Frame Headers and this is handled in the EncodeAndWriteChannelData()
        // method below.

        return true;
    }

    //////////
    //

    bool CCodecVDIF1::ReadNextBlock( void )
    {
        return true;
    }

    //////////
    //

    bool CCodecVDIF1::EncodeAndWriteChannelData( void )
    {
        bool bSuccess = true;      // Assume success for now.

        try
        {
            // Form the VDIF1 data frames (in the encoder using the decoder) and
            // write to the m_SampleData's buffer object for subsequent file output.

            CDataFrameHeader DFH( &m_sDFH, true );

            // For each complex sample, and for the channel set, form the DFH then
            // frame the data.

            // ***WORK OUT HOW TO DRIVE THE TIMESTAMPS***

            // Cf. ref. [1], note 2, pp 6.
            m_sDFH.uSecondsFromEpoch = 0;
            m_sDFH.uRefEpoch         = 0;

            // No. bits in each sample. Complex => 2 x bits for I then Q.
            m_sDFH.uBitsPerSample    = m_iBitsPerSample;

           // The station id is 2 ASCII id chars or an unsigned numeric. If first 8 bits <48 (0x30)
           // then ASCII assumed. Here we use it to represent the antenna number as a numeric.
            m_sDFH.uStationId = m_iAntennaNumber;

            // Representation of the number of channels = { n : #channels = 2^n }.
            assert( m_iNumberOfChannels > 0 );
            m_sDFH.uLog2Channels = static_cast<unsigned int>( floor(log( m_iNumberOfChannels ) / log( 2 )) );

            // Includes the DFH and must be a multiple of 8 bytes; max <= 2^27 bytes.
            //assert( m_iWordSizeInBytes > 0 );
            assert( m_iNumberOfChannels > 0 );
//            m_sDFH.uDataFrameLength = iDFHSize_c + ( m_iNumberOfChannels * m_iWordSizeInBytes );
//            assert( ( m_sDFH.uDataFrameLength % 8 ) == 0 );

            // Number of each time-series of DFs from the same set of channels.
            // For now, assume only a single Data Thread will be required hence
            // leave the thread id set to zero.
            m_sDFH.uThreadId = 0;

            // Data Frame Number is being used as a de-facto sample number. Starting at zero and integral number of
            // data frames per second.
            m_sDFH.uDataFrameNumber = 0;

            // Set the sample data parameters, the total storage size required and then
            // allocate the storage. Recall: encoder uses the m_SampleData object and the
            // decoder's input is accessed via the m_pEncoder->m_SampleData.
//            m_SampleData.SetSampleParams( m_iMode, m_iBitsPerComplexSample, m_iNumberOfChannels, m_iNumberOfSamples );

//            const int iSampleSizeInBytes = 2 * m_iWordSizeInBytes * m_iNumberOfChannels;
//            const int iTotalSizeInBytes  = ( m_sDFH.uDataFrameLength + iSampleSizeInBytes ) * m_iNumberOfSamples;

// TODO (warcus#1#): REWORK BUFFERING HERE

//            if ( m_SampleData.SetStorage( iTotalSizeInBytes ) )
//            {
                // For each time-sample, copy the DFH then all channels for that time-sample.
//                assert( m_SampleData.Size() == iTotalSizeInBytes );

//                ICodec *pDecoder = GetDecoder();
//                assert( pDecoder != nullptr );
//
//                byte_t *pbyInputSamples = pDecoder->GetSampleData().GetSamples();

//                for ( int iSample = 0; iSample < m_iNumberOfSamples; iSample++ )
//                {
//                    // Add the header...
//                    m_SampleData.SetSamples( reinterpret_cast<byte_t *>( &m_sDFH ), iDFHSize_c );
//
//                    // then the channels for this time-sample.
////                    m_SampleData.SetSamples( &pbyInputSamples[ iSample * iSampleSizeInBytes ], iSampleSizeInBytes );
//
//                    ++m_sDFH.uDataFrameNumber;
//                }
//            }
        }
        catch ( ... )
        {
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    bool CCodecVDIF1::WriteConvertedData( void )
    {
        // Extract the queued binary data and write the records to disc.
        // NB: Only an encoder should have final-state data to write.

        bool bSuccess = false;
// TODO (warcus#1#): REWORK BUFFERING HERE
//        try
//        {
//            if ( m_eCodecMode != eCodecModeEncode )
//            {
//                throw string{ "Unexpected codec mode for encoding" };
//            }
//            else
//            {
//                const int iBytes = m_SampleData.Size();
//                assert( iBytes > 0 );
//
//                byte_t *pbySamples = m_SampleData.GetSamples();
//                assert( pbySamples != nullptr );
//
//                if ( ( pbySamples == nullptr ) || ( iBytes <= 0 ) )
//                {
//                    throw string{ "Invalid data access" };
//                }
//                else if ( ! m_pFile->IsOpen() )
//                {
//                    throw string{ "Output file is not open" };
//                }
//                else if ( ! m_pFile->Write( pbySamples, iBytes ) )
//                {
//                    throw string{ "Error trying to write data to file" };
//                }
//
//                bSuccess = true;
//            }
//        }
//        catch ( string rErrorMessage )
//        {
//            fprintf( stderr, "CCodecVDIF1::WriteConvertedData(), %s.\n", rErrorMessage.c_str() );
//            bSuccess = false;
//        }
//        catch ( ... )
//        {
//            fprintf( stderr, "CCodecVDIF1::WriteConvertedData(), Unspecified exception caught.\n" );
//            bSuccess = false;
//        }

        return bSuccess;
    }

    //////////
    //

    void CCodecVDIF1::DumpHeader ( void )
    {
        // Not used.
    }

    //////////
    //

    bool CCodecVDIF1::operator ()( void )
    {
        // Not used.

        return true;
    }

}       // end namespace NCodec


