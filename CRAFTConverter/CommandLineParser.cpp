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
#include "CommandLineParser.h"

#include <array>
#include <cstdio>
#include <getopt.h>

///////////////////////////////////////////////////////////////////////////////
//

CCommandLineParser::CCommandLineParser( void )
{
    m_Options.empty();
    m_Params.empty();
}

//////////
//

CCommandLineParser::~CCommandLineParser( void )
{
    m_Options.empty();
    m_Params.empty();
}

//////////
//

CCommandLineParser::CCommandLineParser( const CCommandLineParser &rRhs )
{
    if ( &rRhs != this )
    {
        m_Options = rRhs.m_Options;
        m_Params  = rRhs.m_Params;
    }
}

//////////
//

CCommandLineParser& CCommandLineParser::operator = ( const CCommandLineParser &rRhs )
{
    if ( &rRhs != this )
    {
        m_Options = rRhs.m_Options;
        m_Params  = rRhs.m_Params;
    }

    return *this;
}

//////////
//

bool CCommandLineParser::Parse( int argc, char **argv, string &rTemplate, struct option *pLongOpts )
{
    bool bSuccess = false;

    try
    {
        int c;
        int iOptionIndex = 0;

        // Loop through the options as set the m_Options map accordingly.

        opterr  = 0;            // Suppress stderr output from getopt_long.

        string sValue;

        while ( ( c = getopt_long( argc, argv, rTemplate.c_str(), pLongOpts, &iOptionIndex ) ) != -1 )
        {
            if ( c != '?' )
            {
                sValue.clear();

                if ( optarg != nullptr )
                {
                    sValue = optarg;
                }

                m_Options[ static_cast<char>( c ) ] = sValue;   // Add the key/option string to the map.
            }

            opterr = 0;         // Reset for next iteration.
        }

        // Now push the non-optional parameters onto out vector for later retrival.

        for ( int iIndex = optind; iIndex < argc; iIndex++ )
        {
            string sParam{ argv[ iIndex ] };
            m_Params.push_back( sParam );
        }

        bSuccess = true;
    }
    catch ( ... )
    {
        bSuccess = false;
    }

    return bSuccess;
}

//////////
//

bool CCommandLineParser::DoesOptionExist( const char &rcValue ) const
{
    auto Search = m_Options.find( rcValue );
    return ( Search != m_Options.end() );
}

//////////
//

string CCommandLineParser::GetOption( const char &rcValue ) const
{
    string sReturnString = "";

    auto Search = m_Options.find( rcValue );

    if ( Search != m_Options.end() )
    {
        sReturnString = Search->second;
    }

    return sReturnString;
}

//////////
//

bool CCommandLineParser::AreThereParameters( void ) const
{
    return ( ! m_Params.empty() );
}

//////////
//

std::vector<std::string> CCommandLineParser::GetParams( void ) const
{
    return m_Params;
}

//////////
//

