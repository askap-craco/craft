///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   FileProcessor.cpp
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

#include "FileProcessor.h"
#include "CodecVCRAFT.h"
#include "CodecVDIF1.h"
#include "CodecCODIF.h"

#include <array>

///////////////////////////////////////////////////////////////////////////////
// Used, defined and/or parts of namespaces.

using namespace NCodec;

///////////////////////////////////////////////////////////////////////////////
// CFileProcessor class implementation.

CFileProcessor::CFileProcessor( void )
{
    m_InputFile.clear();
    m_OutputFile.clear();
    m_iLastError = eErrorCodeNone;
    m_bDumpHeader = false;
}

//////////
//

CFileProcessor::~CFileProcessor( void )
{
    m_bDumpHeader = false;
    m_iLastError = eErrorCodeNone;
    m_OutputFile.clear();
    m_InputFile.clear();
}

//////////
//

CFileProcessor::CFileProcessor( const CFileProcessor &rRhs )
{
    if ( &rRhs != this )
    {
        m_InputFile   = rRhs.m_InputFile;
        m_OutputFile  = rRhs.m_OutputFile;
        m_iLastError  = rRhs.m_iLastError;
        m_bDumpHeader = rRhs.m_bDumpHeader;
    }
}

//////////
//

CFileProcessor & CFileProcessor::operator = ( const CFileProcessor &rRhs )
{
    if ( &rRhs != this )
    {
        m_InputFile  = rRhs.m_InputFile;
        m_OutputFile = rRhs.m_OutputFile;
        m_iLastError = rRhs.m_iLastError;
        m_bDumpHeader = rRhs.m_bDumpHeader;
    }

    return *this;
}

//////////
//

bool CFileProcessor::ProcessFile(   const string   &rInputFileName,
                                    const string   &rOutputFileName,
                                    int             iOutputFileType,
                                    bool            bDumpHeader )
{
    bool bSuccess = false;      // Assume failure for now.

    try
    {
        if ( ! rInputFileName.empty() )
        {
            m_InputFile = rInputFileName;
        }

        if ( ! rOutputFileName.empty() )
        {
            m_OutputFile = rOutputFileName;
        }

        m_bDumpHeader = bDumpHeader;

        // For now, input from VCRAFT and output to CODIF.
        // Later, extend this more generally by propagating the command-line options.

        CFileDescriptor InFile ( m_InputFile,  CFileDescriptor::eFileFormatVCRAFT );
        CFileDescriptor OutFile( m_OutputFile, static_cast<CFileDescriptor::EFileFormat>( iOutputFileType ) );

        if ( ! InFile.Open( CFileDescriptor::eFileModeReadOnly ) )
        {
            throw eErrorCodeOpeningInputFile;
        }
        else if ( ! OutFile.Open( CFileDescriptor::eFileModeReadWrite ) )
        {
            throw eErrorCodeOpeningOutputFile;
        }
        else if ( ! HandleConversion( InFile, OutFile ) )
        {
            throw eErrorCodeFileConversion;
        }
        else
        {
            bSuccess = true;
        }
    }
    catch ( int iErrorCode )
    {
        m_iLastError = iErrorCode;
        bSuccess = false;
        fprintf( stderr, "CFileProcessor::ProcessFile(), Error encountered (code: %d).\n", iErrorCode );
    }
    catch ( ... )
    {
        m_iLastError = eErrorCodeUnhandledException;
        bSuccess = false;
        fprintf( stderr, "CFileProcessor::ProcessFile(), Unhandled exception caught.\n" );
    }

    return bSuccess;
}

//////////
//

int CFileProcessor::GetLastError( void ) const
{
    return m_iLastError;
}

//////////
//

bool CFileProcessor::HandleConversion( CFileDescriptor &rInFile, CFileDescriptor &rOutFile )
{
    EErrorCodes eErrorCode = eErrorCodeNone;

    try
    {
        switch ( rInFile.FileType() )
        {
            // Presently unhandled inpt file types.

            case CFileDescriptor::eFileFormatNone:
            case CFileDescriptor::eFileFormatVDIF:
            case CFileDescriptor::eFileFormatCODIF:
            case CFileDescriptor::eFileFormatRAWPKT:
            case CFileDescriptor::eFileFormatDIFX:
            case CFileDescriptor::eFileFormatFITS:
            case CFileDescriptor::eFileFormatMS:
            case CFileDescriptor::eFileFormatCRAFTCAL:
            default:

                throw eErrorCodeFormatNotSupported;
                break;

            case CFileDescriptor::eFileFormatVCRAFT:

                eErrorCode = ConvertVCRAFTFile( rInFile, rOutFile );

                if ( eErrorCode != eErrorCodeNone )
                {
                    throw eErrorCode;
                }

                break;
        }
    }
    catch ( EErrorCodes eCode )
    {
        fprintf( stderr, "CFileProcessor::HandleConversion(), Error encountered (code: %d).\n", eCode );
        eErrorCode = eCode;
    }
    catch ( ... )
    {
        eErrorCode = eErrorCodeUnhandledException;
    }

    m_iLastError = static_cast< int >( eErrorCode );

    return ( eErrorCode == eErrorCodeNone );
}

//////////
//

CFileProcessor::EErrorCodes
CFileProcessor::ConvertVCRAFTFile( CFileDescriptor &rInFile,
                                   CFileDescriptor &rOutFile )
{
    EErrorCodes eErrorCode = eErrorCodeNone;

    try
    {
        // First create and utilise the decoder to validate and parse the input file.

        CCodecVCRAFT Decoder( &rInFile, CCodecVCRAFT::eCodecModeDecode );

        if ( ! Decoder.ReadAndValidateHeader() )
        {
            throw eErrorCodeInvalidHeader;
        }
        else if ( ! Decoder.DecodeHeader() )
        {
            throw eErrorCodeDecodingHeader;
        }
        else
        {
            // At this point, the input file is deemed valid. Next we use the
            // encoder to read, decode and write the voltages in blocks
            // The decoder's output is used as input to the encoder.

            // Optionally dump the decoded header, if requested.

            DumpHeader( Decoder );

            switch ( rOutFile.FileType() )
            {
                case CFileDescriptor::eFileFormatVDIF:

                    {
#if 0
                        CCodecVDIF1 Encoder( &rOutFile, CCodecVDIF1::eCodecModeEncode );

                        if ( ! EncodeAndWriteOutputFile( Encoder, Decoder ) )
                        {
                            throw Encoder.GetLastError();
                        }
#endif
                    }

                    break;


                case CFileDescriptor::eFileFormatCODIF:

                    {
                        CCodecCODIF Encoder( &rOutFile, CCodecCODIF::eCodecModeEncode );

                        if ( ! EncodeAndWriteOutputFile( Encoder, Decoder ) )
                        {
                            throw Encoder.GetLastError();
                        }
                    }

                    break;

                default:

                    throw eErrorCodeFormatNotSupported;
                    break;  // Cannot be executed given the throw above but included for future completeness.
            }

            eErrorCode = eErrorCodeNone;
        }
    }
    catch ( EErrorCodes eCode )
    {
        fprintf( stderr, "CFileProcessor::ConvertVCRAFTFile(), error code: %d.\n", eCode );
        eErrorCode = eErrorCodeInternalError;
    }
    catch ( ... )
    {
        fprintf( stderr, "CFileProcessor::ConvertVCRAFTFile(), Unhandled exception caught.\n" );
        eErrorCode = eErrorCodeUnhandledException;
    }

    m_iLastError = eErrorCode;

    return eErrorCode;
}

//////////
//

bool CFileProcessor::EncodeAndWriteOutputFile( ICodec &rEncoder,
                                               ICodec &rDecoder )
{
    // Use the decoder to feed the encoder then write the output file.

    bool bSuccess = true;

    printf("Debug: CFileProcessor::EncodeAndWriteOutputFile\n");
	  
    try
    {
        // Prepare the encoder.

        if ( ! rEncoder.ConfigureEncoder( rDecoder ) )
        {
            throw rEncoder.GetLastError();
        }

        if ( ! rEncoder.EncodeHeader() )
        {
            throw rEncoder.GetLastError();
        }

        if ( ! rEncoder.Initialise() )
        {
            throw rEncoder.GetLastError();
        }

        // Want read blocksize to be multiple of output data array. May need to be tweaked if multiple files read

        int iBlockSize = ((2 * 1024 * 1024) / rEncoder.DataArraySize()) * rEncoder.DataArraySize();
        rDecoder.SetBlockSize( iBlockSize );

        // We may need to skip some blocks from the input file to allow alignment

        if ( rEncoder.SkipBytes() > 0 )
        {
	  printf("DEBUG: SkipBytes\n");
	  rDecoder.SeekForward( rEncoder.SkipBytes() );
        }

        // Next, have the Decoder read the sample data in sucessive blocks and
        // the Encoder take this blocked-data, convert format and write to the
        // output file.

        while ( rDecoder.ReadNextBlock() )
        {
            if ( ! rEncoder.EncodeAndWriteChannelData() )
            {
                throw rEncoder.GetLastError();
            }
        }

        // Give the encoder the chance to flush any queued data before the
        // output file is closed and the encoder object destroyed.

        if ( ! rEncoder.Flush() )
        {
            throw rEncoder.GetLastError();
        }
    }
    catch ( const int &riErrorCode )
    {
        bSuccess = false;
    }
    catch ( ... )
    {
        bSuccess = false;
    }

    return bSuccess;
}

//////////
//

void CFileProcessor::DumpHeader( ICodec &rCodec )
{
    if ( m_bDumpHeader )
    {
        rCodec.DumpHeader();
    }
}

//////////
//

