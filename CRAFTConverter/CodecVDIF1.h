///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CodecVDIF1.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Definitions for coding and decoding of VDIF format.
//
//  Notes:      Refer to section 9.4 Multi-channel complex-data Data Array format of
//              [1] re writing multi-channel complex data.
//
//  References: [1] http://vlbi.org/vdif/docs/VDIF_specification_Release_1.1.1.pdf
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CODECVDIF1_H
#define CODECVDIF1_H

#include "Codec.h"

///////////////////////////////////////////////////////////////////////////////
// CCodecVDIF1 class definition as part of the NCodec namespace.

namespace NCodec
{
    class CCodecVDIF1 : public virtual ICodec
    {
        public:

            //////////
            // A bit field to mimic the VDIF1 Data Frame Header. NB: it is a requirement,
            // per ref. [1], that data is written in little endian format.

#pragma pack ( push, 1 )                                // Set single byte alignment.

            typedef struct tagVDIF1DFH
            {
                // Word 0:
                unsigned uInvalidData           : 1;    // Invalid data (taged by source): valid=0, invalid=1.
                unsigned uLegacyMode            : 1;    // 0 - standard VDIF DFH; 1 - legacy header length mode.
                unsigned uSecondsFromEpoch      : 30;   // Seconds from reference epoch; see ref. [1], Note 2, pp 6.

                // Word 1:
                unsigned uReserved              : 2;    // Unassigned - set to 0.
                unsigned uRefEpoch              : 6;    // Reference Epoch for second count.; see ref. [1], Note 2, pp 6.
                unsigned uDataFrameNumber       : 24;   // Data Frame # starting at zero; integral number of Data Frames per second

                // Word 2:
                unsigned uVDIFVersion           : 3;    // VDIF version number.
                unsigned uLog2Channels          : 5;    // Representation (N) of number of channels s.t., # Channels = 2^N.
                unsigned uDataFrameLength       : 24;   // Includes the DFH and must be a multiple of 8 bytes. Max <= 2^27 bytes.

                // Word 3:
                unsigned uDataType              : 1;    // 0 - Real; 1 - Complex data.
                unsigned uBitsPerSample         : 5;    // # bits in each sample. Complex => 2 x bits (I & Q).
                unsigned uThreadId              : 10;   // Number of each time-series of data DFs from the same set of subbands.
                unsigned uStationId             : 16;   // 2 ASCII ID chars or unsigned numeric. If first 8 bits <48 (0x30) then ASCII assumed.

                // Word 4:
                unsigned uExtendedDataVersion   : 8;    // Extended user data from here down - utilisation TBD.
                unsigned uExtendedUserData1     : 24;   // Unique EDV # (assigned via http://www.vlbi.org/vdif/)

                // Words 5 through 7 respectively:
                unsigned uExtendedUserData2     : 32;   // Set to zero if EDV = 0.
                unsigned uExtendedUserData3     : 32;   // Set to zero if EDV = 0.
                unsigned uExtendedUserData4     : 32;   // Set to zero if EDV = 0.
            }
            VDIF1DFH_t;

#pragma pack ( pop )                                    // Restore the previous byte alignment.

            //////////
            // Compile-time size and alignment checks.

            static_assert( alignof( VDIF1DFH_t ) == 1,  "VDIF1DFH_t does not have single byte alignment." );
            static_assert( sizeof ( VDIF1DFH_t ) == 32, "VDIF1DFH_t has an incorrect size." );

            //////
            // Construction, destruction and assignment.

            CCodecVDIF1( void );
            CCodecVDIF1( CFileDescriptor *pFile, ECodecMode eCodecMode );
            virtual ~CCodecVDIF1( void );

            CCodecVDIF1( const CCodecVDIF1 &rRhs ) = delete;
            CCodecVDIF1 & operator = ( CCodecVDIF1 &rRhs ) = delete;

            //////////
            // Public overrides for the standard codec interface.

            bool ReadAndValidateHeader( void );
            bool DecodeHeader( void );
            bool DecodeChannelData( void );

            bool ConfigureEncoder( ICodec &rDecoder );
            bool EncodeHeader( void );
            bool EncodeAndWriteChannelData( void );
            bool ReadNextBlock( void );

            void DumpHeader ( void );
            bool operator ()( void );

        private:

            //////////
            // Private attributes.

            VDIF1DFH_t  m_sDFH;

            //////////
            // Private methods.

            bool WriteConvertedData( void );
    };

}       // end namespace NCodec

#endif // CODECVDIF1_H

