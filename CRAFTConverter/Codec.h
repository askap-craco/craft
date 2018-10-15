///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   Codec.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Defines the data definitions for a base-class codec.
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CCODEC_H
#define CCODEC_H

#include "SampleData.h"
#include "FileDescriptor.h"

///////////////////////////////////////////////////////////////////////////////
// ICodec base-class definition and related helpers.

namespace NCodec
{
    //////////
    // CCodecParameters helper class for generic Codec state.

    class CCodecParameters
    {
        public:

            //////////
            // Public attributes.

            int                 m_iFileHeaderSize;
            int                 m_iBitsPerSample;
            int                 m_iNumberofPol;
            int                 m_iNumberOfWords;
            int                 m_iNumberOfChannels;
            int                 m_iMode;
            int                 m_iBeamId;
            int                 m_iFPGAId;
            int                 m_iCardNumber;
            int                 m_iAntennaNumber;

            double              m_dSampleRate;
            double              m_dMJDNow;

            unsigned long long  m_ullStartWriteFrameId;
            unsigned long long  m_ullStopWriteFrameId;
            unsigned long long  m_ullTriggerFrameId;

            unsigned long long  m_ullStartWriteBAT;
            unsigned long long  m_ullStopWriteBAT;
            unsigned long long  m_ullTriggerWriteBAT;
            unsigned long long  m_ullNowBAT;

            std::string         m_sUTC;
            std::string         m_sFrequencies;

            //////////
            //

            CCodecParameters( void )
            {
                ResetContent();
            }

            //////////
            //

            virtual ~CCodecParameters( void )
            {
            }

            //////////
            //

            CCodecParameters( const CCodecParameters &rRhs ) = default;
            CCodecParameters & operator = ( const CCodecParameters &rRhs ) = default;

        private:

            void ResetContent( void )
            {
                m_iFileHeaderSize       = 0;
                m_iBitsPerSample        = 0;
                m_iNumberOfWords        = 0;
                m_iNumberOfChannels     = 0;
                m_iMode                 = 0;
                m_iBeamId               = 0;
                m_iFPGAId               = 0;
                m_iCardNumber           = 0;
                m_iAntennaNumber        = 0;

                m_dSampleRate           = 0.0;
                m_dMJDNow               = 0.0;

                m_ullStartWriteFrameId   = 0;
                m_ullStopWriteFrameId  = 0;
                m_ullTriggerFrameId      = 0;
                m_ullStartWriteBAT      = 0;
                m_ullStopWriteBAT     = 0;
                m_ullTriggerWriteBAT    = 0;
                m_ullNowBAT             = 0;

                m_sUTC.clear();
                m_sFrequencies.clear();
            }
    };

    ///////////////////////////////////////////////////////////////////////////
    // An abstract base-class to form a common interface for a generic codec
    // and supporting core functionality. Derived codecs need to
    // implement the pure virtual methods hereinafter.

    class ICodec : public virtual CCodecParameters
    {
        public:

            //////////
            // Public attributes.

            enum ECodecMode
            {
                eCodecModeDecode,
                eCodecModeEncode,
            };

            enum ECodecError
            {
                eCodecErrorNone                 = 0,
                eCodecErrorInvalidHeader,
                eCodecErrorInFileInvalid,
                eCodecErrorOutFileInvalid,
                eCodecErrorFileReadFailed,
                eCodecErrorFileWriteFailed,
                eCodecErrorDecodeHeaderFailed,
                eCodecErrorUnspecified
            };

            //////////
            // Construction, destruction and assignment.

            ICodec( CFileDescriptor *pFile      = nullptr,
                    ECodecMode       eCodecMode = eCodecModeDecode )
                :CCodecParameters()
            {
                m_bHeaderValid = false;
                m_eCodecMode   = eCodecMode;
                m_eErrorCode   = eCodecErrorNone;
                m_pDecoder     = nullptr;
                m_pFile        = pFile;
            }

            //////////
            //

            virtual ~ICodec( void )
            {
                m_pFile        = nullptr;
                m_pDecoder     = nullptr;
                m_eErrorCode   = eCodecErrorNone;
                m_eCodecMode   = eCodecModeDecode;
                m_bHeaderValid = false;
            }

            //////////
            // Copy and assignment not implemented.

            ICodec( const ICodec &rRhs ) = delete;
            ICodec & operator = ( const ICodec &rRhs ) = delete;

            //////////
            //

            void SetFileDescriptor( CFileDescriptor &rFile )
            {
                m_pFile = &rFile;
            }

            //////////
            //

            void SetDecoder( ICodec *pDecoder )
            {
                m_pDecoder = pDecoder;
            }

            //////////
            // Casts an object of this type as a pointer to the parameter base-class
            // for direct access thereof.

            operator CCodecParameters * ( void )
            {
                return static_cast<CCodecParameters *>( this );
            }

            //////////
            // Used to test whether a decoded header has been deemed to be valid.

            bool IsDecodedHeaderValid( void ) const
            {
                return m_bHeaderValid;
            }

            //////////
            //

            ECodecError GetLastError( void ) const
            {
                return m_eErrorCode;
            }

            //////////
            // Whether the codec is acting as a decoder or encoder.

            ECodecMode GetMode( void ) const
            {
                return m_eCodecMode;
            }

            //////////
            //

            CFileDescriptor * GetFileDescriptor( void ) const
            {
                return m_pFile;
            }

            //////////
            //

            CSampleData & GetSampleData( void )
            {
                return m_SampleData;
            }

            //////////
            //

            ICodec * GetDecoder( void )
            {
                return m_pDecoder;
            }

            //////////
            // Virtual methods to be overridden.

            //////////
            // Decoder methods.

            virtual bool ReadAndValidateHeader( void )        = 0;
            virtual bool DecodeHeader( void )                 = 0;
            virtual bool ReadNextBlock( void )                = 0;

            virtual bool DecodeChannelData( void )      // Obsolete.
            { return true; }

            //////////
            // Encoder methods.

            virtual bool ConfigureEncoder( ICodec &rDecoder ) = 0;
            virtual bool EncodeHeader( void )                 = 0;
            virtual bool EncodeAndWriteChannelData( void )    = 0;

            //////////
            // Virtual methods that may (optionally) be overridden.

            virtual void DumpHeader( void )            {}
            virtual bool operator () ( void )          { return true; }
            virtual bool Flush( void )                 { return true; }
            virtual int  DataArraySize( void )	       { return -1;   }
            virtual bool Initialise( void )            { return true; }
            virtual int  SkipBytes( bool *preload )    { UNREFERENCED_PARAMETER( preload ); return 0; }
            virtual void setPreload( bool preload )    { UNREFERENCED_PARAMETER( preload ); }
            virtual bool SetBlockSize( int iBlockSize )
            {
                UNREFERENCED_PARAMETER( iBlockSize );
                return true;
            }

            virtual bool SeekForward( int iSkipBytes )
            {
                UNREFERENCED_PARAMETER( iSkipBytes );
                return true;
            }

        protected:

            //////////
            // Protected attributes.

            bool                m_bHeaderValid;
            ECodecMode          m_eCodecMode;
            ECodecError         m_eErrorCode;

            CSampleData         m_SampleData;
            ICodec             *m_pDecoder;
            CFileDescriptor    *m_pFile;
    };

}       // end namespace NCodec.

#endif // CCODEC_H

