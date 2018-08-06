///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   VCRAFTParams.cpp
//
//  Author(s):  Wayne Arcus
//
//  Purpose:
//
//  Notes:      When time permits, rewrite the parsing routines herein using C++ v11's
//              Regex Library. A possible template may be this:
//              "'(\S+) (\S+) # (\S*)\Ŕn'".
//              cf. http://en.cppreference.com/w/cpp/regex.
//
//  References:
//
///////////////////////////////////////////////////////////////////////////////

#include "StdApp.h"
#include "VCRAFTParser.h"

#include <cstring>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
//CVCRAFTParser class implementation (as part of the NCodec namespace).

namespace NCodec
{
    //////////
    // Construction and destruction.

    CVCRAFTParser::CVCRAFTParser( void )
    {
        m_ParamMap.clear();
        LoadParameterList();
    }

    //////////
    //

    CVCRAFTParser::~CVCRAFTParser( void )
    {
        m_ParamMap.clear();
    }

    //////////
    //

    string CVCRAFTParser::RetrieveParameterString( char const *pszKey )
    {
        string sKey{pszKey};
        string sParam;

        ParamMapIterator_t pItem = m_ParamMap.find( sKey );

        if ( pItem != m_ParamMap.end() )
        {
            sParam = pItem->second.m_sParameter;
        }

        return sParam;
    }

    //////////
    //

    bool CVCRAFTParser::Decode( byte_t *pbyHeader )
    {
        // Here we decode the header based on the list of keys and
        // parameter types in the map. Specifically, cycle through
        // the keys and extract the fields placing them into the map.

        bool bSuccess = true;

        try
        {
            assert( pbyHeader != nullptr );

            for ( ParamMapIterator_t pItem = m_ParamMap.begin();
                    pItem != m_ParamMap.end();
                        pItem++ )
            {
                if ( ! ParamsFromHeader( pItem, pbyHeader ) )
                {

                    // It is probable that additional parameters are added which
                    // is not an error, per se, so simply flag this fact.
                    fprintf( stderr, "Unrecognised VCRAFT parameter found.\n" );
                }
            }
        }
        catch ( ... )
        {
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    string CVCRAFTParser::GetSingleParamOnly( const string &rsIdentifier, const byte_t *pbyHeader )
    {
        // This routine extracts a single parameter (token) as a string immediately following
        // the given identifier.

        // Scan the byte stream for the pattern:
        // identifier_parameter_#_comment<lf> where '_' is a space and <lf> = 0x0A, a line feed.

        string sParameter{ "" };

        try
        {
            const char *pResult  = reinterpret_cast<const char *>( pbyHeader );
            const char *pTarget1 = rsIdentifier.c_str();

            // Find the target substring in the header. If found, decode the token,
            // delimited by a # then extract it as a string. Ignore the comment for now.

            if ( ( pResult = strstr( pResult, pTarget1 ) ) != NULL )
            {
                // Now search past the Identifier and its following space up to
                // the terminating delimiter, viz., _#_.

                const char *pStartOfSearch     = pResult + rsIdentifier.size() + 1;
                const int   iPointerOffset     = pResult - reinterpret_cast<const char *>( pbyHeader );
                const int   iRemainderOfHeader = iVCRAFTFileHeaderSizeInBytes_c - iPointerOffset;

                // Convert the header from the search position to the end to a string
                // and then look for the next terminator (i.e., _#_ ).

                string sSearchString( pStartOfSearch, iRemainderOfHeader );
                size_t iPos = sSearchString.find( " # " );

                // If found, the substring we are looking for should be the first item
                // through to one character before the terminator.

                if ( iPos != string::npos )
                {
                    sParameter = sSearchString.substr( 0, iPos );
                }
            }
        }
        catch ( ... )
        {
            sParameter.clear();
            fprintf( stderr, "CVCRAFTParser::GetSingleParamOnly(), Exception caught!\n" );
        }

        return sParameter;
    }

    //////////
    //

    void CVCRAFTParser::LoadParameterList( void )
    {
            // Parameter keys and types expected.

            static const std::vector<ParamItem_t> ParamTemplate =
            {
                { "HDR_SIZE",             eParamTypeInt,     "", ""  },
                { "DATA_TYPE",            eParamTypeString,  "", ""  },
                { "CMD_ARGS",             eParamTypeString,  "", ""  },
                { "CMD_NAMESPACE",        eParamTypeString,  "", ""  },
                { "SAMP_RATE",            eParamTypeDouble,  "", ""  },
                { "CRAFT_MODE",           eParamTypeInt,     "", ""  },
                { "ANT_RA",               eParamTypeDouble,  "", ""  },
                { "ANT_DEC",              eParamTypeDouble,  "", ""  },
                { "ANT_AZ",               eParamTypeDouble,  "", ""  },
                { "ANT_EL",               eParamTypeDouble,  "", ""  },
                { "NBITS",                eParamTypeInt,     "", ""  },
                { "NPOL",                 eParamTypeInt,     "", ""  },
                { "BEAM",                 eParamTypeInt,     "", ""  },
                { "FPGA_ID",              eParamTypeInt,     "", ""  },
                { "CARD_NO",              eParamTypeInt,     "", ""  },
                { "ANTENNA_NO",           eParamTypeInt,     "", ""  },
                { "NCHANS",               eParamTypeInt,     "", ""  },
                { "NOW_MJD",              eParamTypeDouble,  "", ""  },
                { "NOW_BAT",              eParamTypeHex,     "", ""  },
                { "START_WRITE_FRAMEID",  eParamTypeInt,     "", ""  },
                { "STOP_WRITE_FRAMEID",   eParamTypeInt,     "", ""  },
                { "TRIGGER_FRAMEID",      eParamTypeInt,     "", ""  },
                { "START_WRITE_BAT",      eParamTypeHex,     "", ""  },
                { "STOP_WRITE_BAT",       eParamTypeHex,     "", ""  },
                { "TRIGGER_BAT",          eParamTypeHex,     "", ""  },
                { "FREQS",                eParamTypeList,    "", ""  },
                { "UTC_NOW",              eParamTypeISODate, "", ""  }
            };

            // Now create the key and parameter types anew in the map ahead of decoding the parameters.

            m_ParamMap.clear();

            for ( const ParamItem_t &rItem : ParamTemplate )
            {
                m_ParamMap[ rItem.m_sKey ] = rItem;
            }
    }

    //////////
    //

    bool CVCRAFTParser::GetSingleParam( ParamMapIterator_t pItem, const byte_t *pbyHeader )
    {
        // This method extracts a single parameter (token) as a string
        // immediately following the given identifier then the comment field,
        // if set.

        bool bSuccess = false;      // Assume failure for now.

        try
        {
            // Scan the byte stream for the pattern:
            // identifier_parameter_#_comment<lf> where '_' is a space. For example, the UTC_NOW tag thus:
            // "2016-08-04T07:09:24.161595+00:00 # UTC date stamp for when the file was written (ISO format)"

            string    sIdentifier( pItem->second.m_sKey );
	    sIdentifier += ' ';
            const char     *pResult  = reinterpret_cast<const char *>( pbyHeader );
            const char     *pTarget1 = sIdentifier.c_str();

            // Find the target substring in the header. If found, decode the token,
            // delimited by a # then extract the parameter and the comment into pItem.

            if ( ( pResult = strstr( pResult, pTarget1 ) ) != NULL )
            {
                // Now search past the Identifier and its following space up to
                // the terminating delimiter, viz., '_#_'.

                const char *pStartOfSearch     = pResult + pItem->second.m_sKey.size() + 1;
                const int   iPointerOffset     = pResult - reinterpret_cast<const char *>( pbyHeader );
                const int   iRemainderOfHeader = iVCRAFTFileHeaderSizeInBytes_c - iPointerOffset;

                // Convert the header from the search position to the end to a string
                // and then look for the next terminator - i.e., '_#_'.

                string sSearchString( pStartOfSearch, iRemainderOfHeader );
                size_t iPos1 = sSearchString.find( " # " );

                // If found, the substring we are looking for should be the first
                // character before the terminator.

                if ( iPos1 != string::npos )
                {
                    pItem->second.m_sParameter = sSearchString.substr( 0, iPos1 );

                    // The comment will follow the Parameter and be between the "_#_" and
                    // line feed (i.e., <lf> = 0x0A) delimiters. NB: the offest of +3 is to point
                    // one character past the "_#_" delimiter.

                    const int iCharOffset = 3;

                    size_t iPos2 = sSearchString.find( static_cast<char>( '\n' ), iPos1 + iCharOffset );

                    // Treat the comment field as optional and if not found then ensure the corresponding
                    // item in the map is cleared.

                    if ( iPos2 == string::npos )
                    {
                        pItem->second.m_sComment.clear();
                    }
                    else
                    {
                        pItem->second.m_sComment = sSearchString.substr( iPos1 + iCharOffset,
                                                                         iPos2 - iPos1 - iCharOffset );
                    }

                    bSuccess = true;
                }
            }
        }
        catch ( ... )
        {
            fprintf( stderr, "CVCRAFTParser::GetSingleParam(), Exception caught!\n" );
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    //

    CVCRAFTParser::ParamItem_t CVCRAFTParser::RetrieveParameter( string &rsKey )
    {
        ParamItem_t sParamItem;

        if ( ! rsKey.empty() )
        {
            ParamMapIterator_t pItem = m_ParamMap.find( rsKey );

            if ( pItem != m_ParamMap.end() )
            {
                sParamItem = pItem->second;
            }
        }

        return sParamItem;
    }

    //////////
    //

    bool CVCRAFTParser::ParamsFromHeader( ParamMapIterator_t pItem, byte_t *pbyHeader )
    {
        // Perform a full decode of the parameter(s) and comment field for the map item.

        bool bSuccess = true;               // Assume success for now.

        try
        {
  	     switch ( pItem->second.m_eParamType )
            {
                case eParamTypeInt:         // Next single token as an int.
                case eParamTypeDouble:      // Next single token as a double.
                case eParamTypeHex:         // Next single token as a hexidecimal value, prefix 0x.
                case eParamTypeString:      // String between.
                case eParamTypeISODate:     // ISO Date/Time as a string.
                case eParamTypeList:        // Comma delimited string of items.

                    bSuccess = GetSingleParam( pItem, pbyHeader );
                    break;

                default:

                    bSuccess = false;       // Unexpected type.
                    break;
            }
        }
        catch ( ... )
        {
            bSuccess = false;
        }

        return bSuccess;
    }

    //////////
    // Static method.

    bool CVCRAFTParser::ValidateHeader( byte_t *pbyHeader, int iHeaderSize )
    {
        bool bSuccess = false;      // assume failure for now.

        assert( pbyHeader != nullptr );
        assert( iHeaderSize > 0 );

        // Uses a simple validation scheme, viz., that the encoded header size should
        // match the expected size and the self identifier therein (DATA_TYPE) is
        // "CRAFT_VOLTAGES".

        string  sIdentifier = "HDR_SIZE";
        string  sParameter  = GetSingleParamOnly( sIdentifier, pbyHeader );

        if ( ! sParameter.empty() )
        {
            if ( iHeaderSize == atoi( sParameter.c_str() ) )
            {
                sIdentifier = "DATA_TYPE";
                sParameter  = GetSingleParamOnly( sIdentifier, pbyHeader );

                if ( sParameter == "CRAFT_VOLTAGES" )
                {
                    bSuccess = true;
                }
            }
        }

        return bSuccess;
    }

    //////////
    //

}       // end namespace NCodec.
