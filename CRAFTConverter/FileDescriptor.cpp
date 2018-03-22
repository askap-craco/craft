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

bool CFileDescriptor::SeekForward( const int &SeekPosition )
{
  printf("CFileDescriptor::SeekForward: Seeking forward %d bytes\n", SeekPosition);
  if ( IsOpen() )
    {
      printf("   Currently at %ld\n", std::ftell(m_pFile));
      int status = std::fseek( m_pFile, SeekPosition, SEEK_CUR );
      printf("   Now at %ld\n", std::ftell(m_pFile));
        return (status  == 0 );
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

int  CFileDescriptor::Read( WordDeque_t& rQueue, const int &riNumberOfWords,
                            const int &riOffset,  const int &riOrigin )
{
    int iWordsRead = 0;
    uint32_t iWord;
    size_t nRead;

    
    if ( IsOpen() && ( fseek( m_pFile, riOffset, riOrigin ) == 0 ) )
    {
        for ( int i = 0; i < riNumberOfWords; i++ )
        {
	  nRead = std::fread(&iWord, sizeof(uint32_t), 1, m_pFile);
	  if (nRead !=1) {
            if ( feof(m_pFile ) )
	      {
                break;
	      }
	      else
		{
		  iWordsRead = -1;
		}
	  } else {
                rQueue.push_back( iWord );
                iWordsRead++;
	  }
        }
    }

    return iWordsRead;
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


bool CFileDescriptor::Write( const uint32_t * pbyRecord, const int &riNumWords )
{
    bool bSuccess = true;

    if ( riNumWords <= 0 )
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
        std::size_t WordsWritten = std::fwrite( pbyRecord, sizeof(uint32_t), riNumWords, m_pFile );

        bSuccess = ( WordsWritten == static_cast<size_t>( riNumWords ) );
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

