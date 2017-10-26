///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:
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

#ifndef CRAFTMODE_H
#define CRAFTMODE_H

///////////////////////////////////////////////////////////////////////////////
// CCRAFTMode - part of the NCodec namespace.

namespace NCodec
{
    class CCRAFTMode
    {
        private:

            ////////////
            // Private attributes.

            int     m_iVCRAFTDataSize;
            int     m_iWordSize;

        public:

            //////////
            //

            CCRAFTMode( void ) = delete;
            CCRAFTMode( const int &riMode, const int &riNumSamples, const int &riNumChannels );
            virtual ~CCRAFTMode( void );

            //////////
            //

            int TotalVCRAFTBufferSize( void ) const;
            int WordSize( void ) const;
    };

}   // end NCodec namespace.

#endif // CRAFTMODE_H
