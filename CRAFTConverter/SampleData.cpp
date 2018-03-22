///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   SampleData.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:
//
//  Notes:      The internal buffer m_VoltageSamples is a vector<byte_t> that
//              stores the sample data in native VCRAFT format. We allocate the
//              maximum capacity initially then use all or part of that storage
//              thereafter.
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "SampleData.h"

#include <cstring>

///////////////////////////////////////////////////////////////////////////////
//  CSampleData class implementation (part of the NCodec namespace).

namespace NCodec
{
    CSampleData::CSampleData ( void )
    {
        m_iMode               = 0;
        m_iBitsPerSample      = 0;
        m_iNumberOfChannels   = 0;

        m_VoltageSamples.clear();
    }

    //////////
    //

    CSampleData::~CSampleData( void )
    {
        m_VoltageSamples.clear();

        m_iNumberOfChannels   = 0;
        m_iBitsPerSample      = 0;
        m_iMode               = 0;
    }

    //////////
    //

    void CSampleData::SetSampleParams( const int &riMode,
                                       const int &riBitsPerComplexSample,
                                       const int &riNumberOfChannels )
    {
        m_iMode                 = riMode;
        m_iBitsPerSample        = riBitsPerComplexSample/2;
        m_iNumberOfChannels     = riNumberOfChannels;
    }

    //////////
    //

    void CSampleData::GetSampleParams( int &riMode,
                                       int &riBitsPerComplexSample,
                                       int &riNumberOfChannels )
    {
        riMode                  = m_iMode;
        riBitsPerComplexSample  = m_iBitsPerSample*2;
        riNumberOfChannels      = m_iNumberOfChannels;
    }

    //////////
    // Retrieves a reference to the deque of the byte data.

    WordDeque_t & CSampleData::GetSamples( void )
    {
        return m_VoltageSamples;
    }

    //////////
    //

    void CSampleData::SetSamples( const byte_t *pbyData, const int &riBytes )
    {
        if ( ( pbyData != nullptr ) && ( riBytes > 0 ) )
        {
            for ( int iByte = 0; iByte < riBytes; iByte++ )
            {
                m_VoltageSamples.push_back( pbyData[ iByte ] );
            }
        }
    }

    //////////
    //

    int CSampleData::Size( void ) const
    {
        return m_VoltageSamples.size();
    }

    //////////
    //

}       // end namespace NCodec.

