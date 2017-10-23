///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CodecVCRAFT.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Defines the VCRAFT codec class (part of the NCodec namespace).
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CODECVCRAFT_H
#define CODECVCRAFT_H

#include "Codec.h"
#include "VCRAFTParser.h"

#include <array>
#include <memory>

///////////////////////////////////////////////////////////////////////////////
// CCodecVCRAFT class definition (part of the NCodec namespace).

namespace NCodec
{
    class CCodecVCRAFT : public virtual ICodec
    {
        public:

            //////
            // Construction and destruction.

            CCodecVCRAFT( CFileDescriptor *pFile,
                          ECodecMode       eCodecMode = eCodecModeDecode );
            virtual ~CCodecVCRAFT( void );

            //////
            // Default and copy construction and assignment not supported.

            CCodecVCRAFT( void ) = delete;
            CCodecVCRAFT( const CCodecVCRAFT &rRhs ) = delete;
            CCodecVCRAFT & operator = ( CCodecVCRAFT &rRhs ) = delete;

            //////////
            // Overrides for virtual methods from the interface base-class.

            // Decoding - pure virtuals, must be overridden.

            bool ReadAndValidateHeader( void );
            bool DecodeHeader( void );

            // Encoding - pure virtuals, must be overridden.

            bool ConfigureEncoder( ICodec &rDecoder );
            bool EncodeHeader( void );
            bool EncodeAndWriteChannelData( void );
            bool ReadNextBlock( void );

            // Optional overrides (non-pure virutals).

            void DumpHeader ( void );
            bool operator ()( void );
            bool SetBlockSize( int iBlockSize );
            bool SeekForward( int iSkipBytes );

        private:

            //////////
            // Private attributes.

            using HeaderStream_t = std::array<byte_t, iVCRAFTFileHeaderSizeInBytes_c>;

            HeaderStream_t  m_aHeaderStream;        // Raw file header as a byte stream.
            CVCRAFTParser   m_HeaderDecoder;        // Helper for header decoding.
            bool            m_bBuffersInitialised;  // Buffering configured flag.
            int             m_iInputBlockSize;      // Size of read chunk

            //////////
            // Private methods.

            template <typename tValue>
            void RetrieveParameter( char const *pszKey, tValue &rtValue );

            bool SetHeaderParameters( void );

    };

}       // end namespace NCodec.


#endif // CODECVCRAFT_H

