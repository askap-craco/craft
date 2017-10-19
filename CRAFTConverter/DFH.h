///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   DFH.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Light-weight wrapper for the codifio module.
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DFH_H
#define DFH_H

// Temporarily disable various warnings leaked from the external codifio module
// until resolved.

#if defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wtype-limits"
    #pragma GCC diagnostic ignored "-Wsign-compare"
    #include "codifio.h"
    #pragma GCC diagnostic pop
#else
    #include "codifio.h"
#endif

///////////////////////////////////////////////////////////////////////////////
// CDFH class definition and supporting helpers.

namespace NCodec        // Part of the Codec namespace.
{
    using CODIFDFH_t = codif_header;        // Alias for the external struct.

    // At compile time, check the header is the correct size.

    constexpr int iDFHSizeExpected_c = 64;
    constexpr int iDFHSize_c         = sizeof( CODIFDFH_t );

    static_assert( iDFHSize_c == iDFHSizeExpected_c, "CODIF DFH has an unexpected size." );

    //////////
    // CDFH class definition.

    class CDFH : public CODIFDFH_t
    {
        public:

            //////////
            // Construction, destruction, copy and assignment.

            CDFH( void );
            virtual ~CDFH( void );
            CDFH & operator = ( const CDFH &rRhs );
            CDFH ( const CDFH &rRhs );

            //////////
            // Public methods.

            void Reset( void );

            void SetMaxDataFrameNumber( int iMaxDataFrame );
            void ResetDataFrameNumber( void );

            void NextFrame( void );
            void SetFrameNumber( int uiFrameNumber );
            bool SetFrameTime( const double &rdMJD );

            // Cast operators for pointer access the base-class structure.

            operator CODIFDFH_t * ( void );
            operator byte_t *     ( void );

        private:

            //////////
            // Private attributes

            unsigned int m_uMaxDataFramePlusOne;
    };

}               // End Codec namespace.

#endif // DFH_H
