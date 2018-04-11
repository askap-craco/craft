///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CCodecCODIF.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Defines the CODIF Coder/decoder class.
//
//  Notes:      Only encoding of CODIF formatted output is presently supported.
//
//  References: see CodecCODIF.cpp header.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CODECCODIF_H
#define CODECCODIF_H

#include "Codec.h"
#include "DFH.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////
// CCodecCODIF class definition, part of the NCodec namespace.

namespace NCodec
{
    //////////
    // Forward definitions - internal helper definitions.

    //////////
    // CCodecCODIF class definition.

    class CCodecCODIF : public virtual ICodec
    {
        public:

            //////
            // Construction and destruction.

            CCodecCODIF( CFileDescriptor *pFile, ECodecMode eCodecMode );
            virtual ~CCodecCODIF( void );

            //////
            // Copy construction and assignment not supported.

            CCodecCODIF( void ) = delete;
            CCodecCODIF( const CCodecCODIF &rRhs ) = delete;
            CCodecCODIF & operator = ( CCodecCODIF &rRhs ) = delete;

            //////////
            // Overrides of the virtual and pure virtual interface methods of ICodec.

            //////////
            // Decoding.

            bool ReadAndValidateHeader( void );
            bool DecodeHeader( void );
            bool DecodeChannelData( void );

            //////////
            // Encoding.

            bool ConfigureEncoder( ICodec &rDecoder );
            bool EncodeHeader( void );
            bool EncodeAndWriteChannelData( void );
            bool ReadNextBlock( void );

            //////////
            // Optional overrides.

            void DumpHeader ( void );
            bool operator ()( void );
            bool Flush( void );
            int  DataArraySize( void );
            bool Initialise( void );
            int  SkipBytes( bool* );

        private:

            //////////
            // Private definitions and attributes.

            using Buffer_t = std::vector<uint32_t>;

            bool     m_bFirstInputBlock;
            bool     m_bSynced;
            bool     m_bInitialiseFrameTime;

            int      m_iSampleBlockSize;
            int      m_iDataArraySize;
            int      m_iDataArrayWords;
            int      m_iDataFrameSize;
            int      m_iNumberOfSyncRetries;

            CDFH     m_DFH;
            unsigned long long m_ullBAT0;
            unsigned long long m_ullFrame0;
            int      m_iSkipSamples;
            int      m_iSampleOffset;
	    int      m_iSamplesPerWord;
            Buffer_t m_DataFrameBuffer;
	    char     *buf;
	    int      mask;  

            //////////
            // Private methods.

            bool SyncAndHandlePreSyncSamples( void );
            bool HandleCODIFFrameData( void );
            bool ConfigureDFH( void );
            bool SetBufferSize( Buffer_t &rBuffer, const int &riLength );
            WordDeque_t & SampleData( void ) const;
            bool WriteDataFrames( bool bForceFlush = false );
            bool WritePreambleFrames( const int &riPreambleSamples );
            int  BytesPerTimeSampleForChannels( void ) const;
            bool SetPartialFrameParams( const int &riFrames );
	    void decodeVCRAFTBlock(WordDeque_t & rInput, std::vector<uint32_t>& vcraftData, std::vector<uint32_t>& codifData,
				   int wordstoUnpack, int samplePerOutword, int samlpesPerWord, int *iWordCount);


    };

}       // end namespace NCodec

#endif // CODECCODIF_H

