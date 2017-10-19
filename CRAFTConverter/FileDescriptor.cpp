///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   FileDescriptor.cpp
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
#include "FileDescriptor.h"

#include <cstdio>

///////////////////////////////////////////////////////////////////////////////
//

CFileDescriptor::CFileDescriptor( void )
{
    m_eIOType            = eIOTypeInput;
    m_eFileFormat        = eFileFormatNone;
    m_bFileTypeSupported = false;
    m_pFile              = nullptr;

    m_FileName.clear();
}

//////////
//

CFileDescriptor::CFileDescriptor( const string &rFileName, const EFileFormat &reFileFormat )
{
    m_FileName           = rFileName;
    m_eFileFormat        = reFileFormat;
    m_bFileTypeSupported = false;
    m_eIOType            = eIOTypeInput;
    m_pFile              = nullptr;
}

//////////
//

CFileDescriptor::~CFileDescriptor( void )
{
    m_eIOType               = eIOTypeInput;
    m_eFileFormat           = eFileFormatNone;
    m_bFileTypeSupported    = false;

    Close();

    m_FileName.clear();

    m_pFile = nullptr;
}

//////////
//

bool CFileDescriptor::Open( const EFileMode &rMode )
{
    bool bSuccess = false;          // Assume failure for now.

    assert( ! m_FileName.empty() );

    if ( m_pFile == nullptr )
    {
        string Mode{ "rb" };

        if ( rMode == eFileModeReadWrite )
        {
            Mode = "wb";
        }

        //fprintf( stderr, "File: %s, Mode: %s\n", m_FileName.c_str(), Mode.c_str() );

        m_pFile = std::fopen( m_FileName.c_str(), Mode.c_str() );

        bSuccess = ( m_pFile != nullptr );
    }
    else
    {
        bSuccess = true;
    }

    return bSuccess;
}

//////////
//

bool CFileDescriptor::IsOpen( void ) const
{
    return ( m_pFile != nullptr );
}

//////////
//

CFileDescriptor::EFileFormat CFileDescriptor::FileType( void ) const
{
    return m_eFileFormat;
}

//////////
//

bool CFileDescriptor::SeekPastFileHeader( const int &riChannelSeekPosition )
{
    if ( IsOpen() )
    {
        return ( std::fseek( m_pFile, riChannelSeekPosition, SEEK_SET ) == 0 );
    }

    return false;
}

//////////
//

int CFileDescriptor::Read( byte_t *pbyReadArray, const int &riNumberOfBytes,
                           const int &riOffset,  const int &riOrigin )
{
    int iBytesRead = 0;

    if ( IsOpen() && ( fseek( m_pFile, riOffset, riOrigin ) == 0 ) )
    {
        iBytesRead = std::fread( pbyReadArray, sizeof( byte_t ), riNumberOfBytes, m_pFile );
    }

    return iBytesRead;
}

//////////
//

int  CFileDescriptor::Read( ByteDeque_t& rQueue, const int &riNumberOfBytes,
                            const int &riOffset,  const int &riOrigin )
{
    int iBytesRead = 0;

    if ( IsOpen() && ( fseek( m_pFile, riOffset, riOrigin ) == 0 ) )
    {
        for ( int i = 0; i < riNumberOfBytes; i++ )
        {
            int iChar = std::fgetc( m_pFile );  // int needed to handle EOF.

            if ( iChar == EOF )
            {
                break;
            }
            else
            {
                rQueue.push_back( static_cast<byte_t>( iChar ) );
                iBytesRead++;
            }
        }
    }

    return iBytesRead;
}

//////////
//

bool CFileDescriptor::Write( const byte_t * pbyRecord, const int &riNumBytes )
{
    bool bSuccess = true;

    if ( riNumBytes <= 0 )
    {
        bSuccess = true;
    }
    else if ( pbyRecord == nullptr )
    {
        bSuccess = false;
    }
    else if ( ! IsOpen() )
    {
        bSuccess = false;
    }
    else
    {
        std::size_t BytesWritten = std::fwrite( pbyRecord, sizeof(byte_t), riNumBytes, m_pFile );

        bSuccess = ( BytesWritten == static_cast<size_t>( riNumBytes ) );
    }

    return bSuccess;
}

//////////
//

void CFileDescriptor::Close( void )
{
    if ( m_pFile != nullptr )
    {
        fclose( m_pFile );
        m_pFile = nullptr;
    }
}

