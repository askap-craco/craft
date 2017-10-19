///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   DFH.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Implements a light wrapper class for the CSIRO codifio module.
//
//  Notes:      None.
//
//  References: None.
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "DFH.h"
#include "dateutils.h"

#include <cstring>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// CDFH class implementation.

namespace NCodec        // Part of the Codec namespace.
{
    CDFH::CDFH( void )
    {
        std::memset( this, 0x00, iDFHSize_c );

        m_uMaxDataFramePlusOne = 0;
    }

    //////////
    //

    CDFH::~CDFH( void )
    {
        m_uMaxDataFramePlusOne = 0;
        std::memset( this, 0x00, iDFHSize_c );
    }

    //////////
    //

    CDFH::CDFH( const CDFH & rRhs )
    {
        if ( &rRhs != this )
        {
            std::memcpy( this, &rRhs, iDFHSize_c );
            m_uMaxDataFramePlusOne = rRhs.m_uMaxDataFramePlusOne;
        }
    }

    //////////
    //

    CDFH & CDFH::operator = ( const CDFH &rRhs )
    {
        if ( &rRhs != this )
        {
            std::memcpy( this, &rRhs, iDFHSize_c );
            m_uMaxDataFramePlusOne = rRhs.m_uMaxDataFramePlusOne;
        }

        return *this;
    }

    //////////
    //

    void CDFH::Reset( void )
    {
        std::memset( this, 0x00, iDFHSize_c );
        m_uMaxDataFramePlusOne = 0;
    }

    //////////
    //

    CDFH::operator CODIFDFH_t * ( void )
    {
        return static_cast<CODIFDFH_t *>( this );
    }

    //////////
    //

    CDFH::operator byte_t * ( void )
    {
        return reinterpret_cast<byte_t *>( this );
    }

    //////////
    //

    void CDFH::SetMaxDataFrameNumber( int iMaxDataFrame )
    {
        m_uMaxDataFramePlusOne = iMaxDataFrame + 1;
    }

    //////////
    //

    void CDFH::ResetDataFrameNumber( void )
    {
        // Data Frame number within the current Period. Starts at zero.

        SetFrameNumber( 0 );
    }

    //////////
    //

    void CDFH::NextFrame( void )
    {
        int iDataFrameNumber = getCODIFFrameNumber( this );

        iDataFrameNumber = ( iDataFrameNumber + 1 ) % m_uMaxDataFramePlusOne;

        setCODIFFrameNumber( this, iDataFrameNumber );
        assert( getCODIFFrameNumber( this ) == iDataFrameNumber );

        const int iNextSecond = getCODIFFrameSecond( this ) + 1;

        setCODIFFrameEpochSecOffset( this, iNextSecond );
    }

    //////////
    //

    void CDFH::SetFrameNumber( int iFrameNumber )
    {
        // Force the data frame number without checking.

        setCODIFFrameNumber( this, iFrameNumber );
    }

    //////////
    //

    bool CDFH::SetFrameTime( const double &rdMJD )
    {
        // Set the time parameters for this sequence of data frames for the current
        // Period.

        bool bSuccess = true;       // Assume success for now.

        try
        {
            double  dMJDWholeSec = 0.0;
            double  dMJDFracSec  = 0.0;

            dMJDFracSec = std::modf( rdMJD, &dMJDWholeSec );
            UNREFERENCED_PARAMETER( dMJDFracSec );

            if ( setCODIFFrameSecond( this, static_cast<int>( dMJDWholeSec ) ) )
            {
                throw string { "setCODIFFrameSecond() failed" };
            }

            // Reset the frame number ready for frame counting within a period.

            SetFrameNumber( 0 );

            bSuccess = true;
        }
        catch ( string sMessage )
        {
            fprintf( stderr, "CDFH::SetFrameTime(), %s.\n", sMessage.c_str() );
            bSuccess = false;
        }
        catch ( ... )
        {
            fprintf( stderr, "CDFH::SetFrameTime(), Unspecified exception.\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //


}               // Part of the Codec namespace.

