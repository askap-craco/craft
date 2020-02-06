///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
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

#include "StdApp.h"

#include "FileProcessor.h"
#include "CommandLineParser.h"
#include "FileDescriptor.h"

//#include <sys/stat.h>
#include <vector>
//#include <fstream>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////
// Namespaces.

namespace       // Anonymous namespace for internal helpers.
{
    //////////
    // Constants.

    enum EReturnCode
    {
        eReturnCodeFailure = -1,
        eReturnCodeSuccess =  0
    };

    const string sDefaultOutputName{ "Default.codif" };

    //////////
    // Functions.

    bool ProcessFiles( string &rsInputName, string &rsOutputName, int iOutputFileType, bool bDumpHeader = false )
    {
        bool bSuccess = false;          // Assume failure for now.

        if ( rsInputName.empty() )
        {
            bSuccess = false;
            fprintf( stderr, "No input file specified.\n" );
        }
        else
        {
            if ( rsOutputName.empty() )
            {
                rsOutputName = sDefaultOutputName;
                fprintf( stderr, "No output file name specified, assuming %s\n", rsOutputName.c_str() );
            }

            CFileProcessor Processor;

            bSuccess = Processor.ProcessFile( rsInputName, rsOutputName, iOutputFileType, bDumpHeader );

            if ( ! bSuccess )
            {
                fprintf( stderr, "ProcessFiles(), File Processor failed. Errorcode: %d\n.", Processor.GetLastError() );
            }
        }

        return bSuccess;
    }
}

///////////////////////////////////////////////////////////////////////////////
// main()

int main( int argc, char **argv )
{
    bool                    bTerminate  = false;
    bool                    bDumpHeader = false;
    EReturnCode             eReturnCode = eReturnCodeFailure;

    static struct option    sLongOptions[] = {{"bat0", 1, 0, 'b'}, { NULL, 0, NULL, 0 } };
    CCommandLineParser      CmdLine;
    string                  sOptionsTemplate{ "d1chb:" };

    string                  sInputName;
    string                  sOutputName;
    int                     iOutputType = CFileDescriptor::eFileFormatCODIF;

    if ( ! CmdLine.Parse( argc, argv, sOptionsTemplate, sLongOptions ) )
    {
        // No command line options specified.

        bTerminate = true;
        fprintf( stderr, "No command line options specified\n" );
    }
    else if ( CmdLine.DoesOptionExist( '?' ) )
    {
        bTerminate  = true;
        eReturnCode = eReturnCodeSuccess;
        fprintf( stderr, "Usage: CRAFTDataConverter [-d -1 -c -? -b <bat0>] InputFile OutputFile\n" );
    }
    else
    {
        if ( CmdLine.DoesOptionExist( 'd' ) )
        {
            bDumpHeader = true;
        }

        // Give preference to VDIF1 format if both '1' and 'c' specified.

        if ( CmdLine.DoesOptionExist( '1' ) )
        {
            iOutputType = CFileDescriptor::eFileFormatVDIF;
        }
        else if ( CmdLine.DoesOptionExist( 'c' ) )
        {
            iOutputType = CFileDescriptor::eFileFormatCODIF;
        }

        if ( CmdLine.DoesOptionExist( 'b' ) )  // Not sure this is useful - leave it in for now, but this is not used elsewhere
        {
	  unsigned long long bat0;
	  std::stringstream ss;
	  ss << std::hex << CmdLine.GetOption('b');
	  ss >> bat0;
	  // printf("BAT0=0x%llX\n", bat0);
        }
	
        // Handle mandatory parameters; only the first two are used.

        if ( CmdLine.AreThereParameters() )
        {
            std::vector<string> Parameters = CmdLine.GetParams();

            if ( Parameters.size() >= 2 )
            {
                sInputName  = Parameters[ 0 ];
                sOutputName = Parameters[ 1 ];
            }
        }
    }

    // If we are successfull thus far, we can process the files else return with
    // an error code.

    if ( ! bTerminate )
    {
        eReturnCode = eReturnCodeSuccess;

        if ( ! ProcessFiles( sInputName, sOutputName, iOutputType, bDumpHeader ) )
        {
            eReturnCode = eReturnCodeFailure;
            fprintf( stderr, "ProcessFiles() failed.\n" );
        }
    }

    return eReturnCode;
}


