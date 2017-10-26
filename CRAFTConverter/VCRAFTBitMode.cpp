///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   VCRAFTBitMode.cpp
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

#include "VCRAFTBitMode.h"

//////////
//

CVCRAFTBitMode::CVCRAFTBitMode( int iMode )
               :m_iModeMask( 0x03 ),
                m_aWordMasks{{ 0xff, 0xff, 0x0f, 0x01 }},
                m_aWordSizes{{    2,    1,    1,    1 }},
                m_aBitWidths{{   16,    8,    4,    1 }}
{
    m_eMode = static_cast<EMode>( iMode & m_iModeMask );
}

//////////
//

CVCRAFTBitMode::~CVCRAFTBitMode( void )
{
    m_eMode = eMode0_16bit;
}

//////////
//

CVCRAFTBitMode::CVCRAFTBitMode( const CVCRAFTBitMode &rRhs )
               :m_iModeMask( 0x03 ),
                m_aWordMasks{{ 0xff, 0xff, 0x0f, 0x01 }},
		m_aWordSizes{{ 2,    1,    1,    1 }},
		m_aBitWidths{{ 16,   8,    4,    1 }}
{
    if ( &rRhs != this )
    {
        m_eMode = rRhs.m_eMode;
    }
}

//////////
//

CVCRAFTBitMode & CVCRAFTBitMode::operator = ( const CVCRAFTBitMode &rRhs )
{
    if ( &rRhs != this )
    {
        m_eMode = rRhs.m_eMode;
    }

    return *this;
}

//////////
//

void CVCRAFTBitMode::SetMode( int iMode )
{
    m_eMode = static_cast<EMode>( iMode & m_iModeMask );
}

//////////
//

CVCRAFTBitMode::EMode CVCRAFTBitMode::GetMode( void ) const
{
    return m_eMode;
}

//////////
//

int CVCRAFTBitMode::ByteMask( void ) const
{
    return m_aWordMasks[ m_eMode ];
}

//////////
//

int CVCRAFTBitMode::WordSize( void ) const
{
    return m_aWordSizes[ m_eMode ];
}

//////////
//

int CVCRAFTBitMode::BitWidth( void ) const
{
    return m_aBitWidths[ m_eMode ];
}

//////////
//
