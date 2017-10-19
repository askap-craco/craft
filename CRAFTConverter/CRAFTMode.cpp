///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
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
#include "CRAFTMode.h"

#include <array>

///////////////////////////////////////////////////////////////////////////////
//

namespace               // Anonymous namespace for internal helpers.
{

    //////////
    // Constants & aliases.

    //constexpr int       iNumberOfVCRAFTChannels_c = 8;
    constexpr int       iNumberOfModes_c          = 4;

    using ModeArray_t = std::array<int, iNumberOfModes_c>;

    const ModeArray_t   aSampleSizeInBytes_c{{ 2, 1, 1, 1 }};
}

namespace NCodec
{
    CCRAFTMode::CCRAFTMode( const int &riMode,
                            const int &riNumSamples,
                            const int &riNumChannels )
    {
        // Compute the derived voltage sample sizes. Mode means:
        // 0 = 16b+16b; 1 = 8b+8b; 2 = 4b+4b; 3 = 1b+1b. Data is complex, I then Q.

        const int iMode = ( riMode & 0x03 );
        assert( iMode < iNumberOfModes_c );

        // Buffer size to fit all channels for all samples. Since the samples are complex
        // we have a factor of 2.

        m_iVCRAFTDataSize = ( riNumSamples * ( riNumChannels * aSampleSizeInBytes_c[ iMode ] * 2 ) );
        m_iWordSize       = ( aSampleSizeInBytes_c[ iMode ] );
    }

    //////////
    //

   CCRAFTMode::~CCRAFTMode( void )
    {
        m_iVCRAFTDataSize = 0;
        m_iWordSize       = 0;
    }

    //////////
    //

    int CCRAFTMode::TotalVCRAFTBufferSize( void ) const
    {
        return m_iVCRAFTDataSize;
    }

    //////////
    //

    int CCRAFTMode::WordSize( void ) const
    {
        return m_iWordSize;
    }

}       // end NCodec namespace.
