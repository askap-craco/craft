///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   VCRAFTBitMode.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:
//
//  Notes:      The VCRAFT data bit depth is encoded in the Mode setting thus:
//              0 = 16b+16b; 1 = 8b+8b; 2 = 4b+4b; 3 = 1b+1b.
//              Data is complex -- I then Q.
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CVCRAFTBITMODE_H
#define CVCRAFTBITMODE_H

#include <array>

class CVCRAFTBitMode
{
    public:

        //////////
        // Public attributes.

        enum EMode
        {
            eMode0_16bit  = 0,  // 16b & 16b
            eMode1_8bit,        //  8b &  8b
            eMode2_4bit,        //  4b &  4b
            eMode3_1bit,        //  1b &  1b

            eModeMaxModes       // Count of items in this enum (keep as last item).
        };

        //////////
        // Public Methods.

        CVCRAFTBitMode( int iMode = 0x00 );
        CVCRAFTBitMode( const CVCRAFTBitMode &rRhs );
        CVCRAFTBitMode & operator = ( const CVCRAFTBitMode &rRhs );
        virtual ~CVCRAFTBitMode( void );

        void    SetMode( int iMode );

        EMode   GetMode( void ) const;
        int     ByteMask( void ) const;
        int     WordSize( void ) const;
        int     BitWidth( void ) const;

    private:

        //////////
        // Private attributes.

        using LUT_t = const std::array<int, eModeMaxModes>;

        const int   m_iModeMask;

        LUT_t       m_aWordMasks;
        LUT_t       m_aWordSizes;
        LUT_t       m_aBitWidths;
        EMode       m_eMode;
};

#endif // CVCRAFTBITMODE_H
