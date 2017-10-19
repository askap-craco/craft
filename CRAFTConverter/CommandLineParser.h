///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   CommandLineParser.h
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

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include <vector>
#include <map>
#include <unistd.h>
#include <getopt.h>

///////////////////////////////////////////////////////////////////////////////
// CCommandLineParser class definition.

class CCommandLineParser
{
    public:

        //////
        // Construction, destruction and assignment.

        CCommandLineParser( void );
        virtual ~CCommandLineParser( void );
        CCommandLineParser & operator = ( const CCommandLineParser &rRhs );
        CCommandLineParser ( const CCommandLineParser &rRhs );

        ////////////
        // Public methods.

        bool Parse( int argc, char **argv, std::string &rTemplate, struct option *pLongOpts );

        std::string GetOption( const char &rcValue ) const;
        bool DoesOptionExist( const char &rcValue ) const;

        std::vector<std::string> GetParams( void ) const;
        bool AreThereParameters( void ) const;

    private:

        ////////////
        // Privte attributes.

        using Options_t = std::map<char, std::string>;
        using Params_t  = std::vector<std::string>;

        Options_t   m_Options;
        Params_t    m_Params;
};

#endif // COMMANDLINEPARSER_H
