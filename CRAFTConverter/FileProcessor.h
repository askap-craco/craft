///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016 - 2017, Wayne Arcus, All Rights Reserved.
//
//  Filename:   FileProcessor.h
//
//  Author(s):  Wayne Arcus
//
//  Purpose:    Defines the CFileProcessor class.
//
//  Notes:
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FILEPROCESSOR_H
#define FILEPROCESSOR_H

#include "FileDescriptor.h"

///////////////////////////////////////////////////////////////////////////////
// CFileProcessor class definition.

class CFileProcessor
{
    public:

        //////////
        // Public attributes.

        enum EErrorCodes
        {
            eErrorCodeNone                  = 0,
            eErrorCodeFileConversion,
            eErrorCodeUnhandledException,
            eErrorCodeOpeningInputFile,
            eErrorCodeOpeningOutputFile,
            eErrorCodeDecodingHeader,
            eErrorCodeEncodingHeader,
            eErrorCodeInvalidHeader,
            eErrorCodeDecodingVoltages,
            eErrorCodeEncodingChannelData,
            eErrorCodeFormatNotSupported,
            eErrorCodeInternalError
        };

        //////
        // Public methods.

        CFileProcessor( void );
        virtual ~CFileProcessor( void );
        CFileProcessor & operator = ( const CFileProcessor &rRhs );
        CFileProcessor ( const CFileProcessor &rRhs );

        ////////////
        // Methods.

        bool ProcessFile(   const string &rInputFileName  = "",
                            const string &rOutputFileName = "",
                            int           iOutputFileType = CFileDescriptor::eFileFormatCODIF,
                            bool          bDumpHeader     = false );

        int GetLastError( void ) const;

    private:

        ////////////
        // Private attributes.

        string  m_InputFile;
        string  m_OutputFile;

        int     m_iLastError;
        bool    m_bDumpHeader;

        ////////////
        // Private methods.

        void DumpHeader( NCodec::ICodec &rCodec );
        bool HandleConversion( CFileDescriptor &rInFile, CFileDescriptor &rOutFile );

        EErrorCodes ConvertVCRAFTFile( CFileDescriptor &rInFile, CFileDescriptor &rOutFile );
        bool EncodeAndWriteOutputFile( NCodec::ICodec &rEncoder, NCodec::ICodec &rDecoder );

};

#endif // FILEPROCESSOR_H
