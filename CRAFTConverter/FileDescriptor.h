///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   FileDescriptor.h
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

#ifndef FILEDESCRIPTOR_H
#define FILEDESCRIPTOR_H

#include <stdio.h>
#include <cstdio>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// CFileDescriptor class definition.

class CFileDescriptor
{
    public:

        //////////
        // Private attributes.

        enum EIOType
        {
            eIOTypeInput,           // Input for conversion.
            eIOTypeOutput           // Output following conversion.
        };

        enum EFileMode
        {
            eFileModeReadOnly,
            eFileModeReadWrite
        };

        enum EFileFormat            // Different file formats recognised.
        {
            eFileFormatNone,        // No format selected.
            eFileFormatVCRAFT,      // Raw volatages from the ASKAP/CRAFT telescope.
            eFileFormatVDIF,        // DIF format.
            eFileFormatCODIF,       // CSIRO Oversampled Data Interchange Format.
            eFileFormatRAWPKT,      // TBA.
            eFileFormatDIFX,        // DIFX - not currently supported.
            eFileFormatFITS,        // FITS - not currently supported.
            eFileFormatMS,          // MS   - not currently supported.
            eFileFormatCRAFTCAL     // CRAFT Calibration data.
        };

        //////////
        // Construction, destruction and assignment.

        CFileDescriptor( void );
        CFileDescriptor( const std::string &rFileName, const EFileFormat &reFileFormat );
        virtual ~CFileDescriptor( void );

        CFileDescriptor( const CFileDescriptor &rRhs ) = delete;
        CFileDescriptor & operator = ( const CFileDescriptor &rRhs ) = delete;

        //////////
        // General public methods.

        bool Open( const EFileMode &rMode = EFileMode::eFileModeReadOnly );
        bool IsOpen( void ) const;

        bool SeekPastFileHeader( const int &riChannelSeekPosition );

        int  Read( byte_t *pbyReadArray, const int &riNumberOfBytes,
                   const int &riOffset,  const int &riOrigin = SEEK_SET );

        int  Read( ByteDeque_t& rQueue, const int &riNumberOfBytes,
                   const int &riOffset,  const int &riOrigin = SEEK_SET );

        bool Write( const byte_t * pbyRecord, const int &riNumBytes );

        EFileFormat FileType( void ) const;

    private:

        //////////
        // Priviate attributes.

        EIOType         m_eIOType;
        EFileFormat     m_eFileFormat;
        bool            m_bFileTypeSupported;
        FILE           *m_pFile;
        std::string     m_FileName;

        //////////
        // Private methods.

        void Close( void );
};

#endif // FILEDESCRIPTOR_H
