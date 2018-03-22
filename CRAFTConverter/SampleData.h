///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   SampleData.h
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

#ifndef SAMPLEDATA_H
#define SAMPLEDATA_H

///////////////////////////////////////////////////////////////////////////////
// CSample class definition (part of the NCodec namespace).

namespace NCodec
{
    class CSampleData
    {
        public:

            //////////
            // Public attributes.


            //////////
            // Construction & destruction.

            CSampleData( void );
            virtual ~CSampleData( void );

            //////////
            // Copy and assignment not implemented.

            CSampleData( const CSampleData &rRhs ) = delete;
            CSampleData & operator = ( const CSampleData &rRhs ) = delete;

            //////////
            // Public members.

            // Set methods.

            void SetSampleParams( const int &riMode, const int &riBitsPerComplexSample,
                                  const int &riNumberOfChannels );

            void SetSamples( const byte_t *pbyData, const int &riBytes );

            // Retrieval methods.

            void GetSampleParams( int &riMode, int &riBitsPerComplexSample,
                                  int &riNumberOfChannels );

            WordDeque_t & GetSamples( void );

            int Size( void ) const;

        private:

            //////////
            // Attributes.

            int         m_iMode;
            int         m_iBitsPerSample;
            int         m_iNumberOfChannels;

            WordDeque_t m_VoltageSamples;
    };

}       // end namespace NCodec

#endif // SAMPLEDATA_H
