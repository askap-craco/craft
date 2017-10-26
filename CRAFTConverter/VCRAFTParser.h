///////////////////////////////////////////////////////////////////////////////
//
//  Copyright © 2016, Wayne Arcus, All Rights Reserved.
//
//  Filename:   VCRAFTParser.h
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

#ifndef VCRAFTPARSER_H_INCLUDED
#define VCRAFTPARSER_H_INCLUDED

#include <map>

///////////////////////////////////////////////////////////////////////////////

namespace NCodec
{
    //////////
    //

    constexpr int iVCRAFTFileHeaderSizeInBytes_c = 4096;

    //////////
    // CVCRAFTParser class definition.

    class CVCRAFTParser
    {
        public:

            //////////
            // Construction and destruction.

            CVCRAFTParser( void );
            virtual ~CVCRAFTParser( void );

            //////////
            // Copy and assignment not supported.

            CVCRAFTParser & operator = ( const CVCRAFTParser &rRhs ) = delete;
            CVCRAFTParser ( const CVCRAFTParser &rRhs ) = delete;

            //////
            // General public methods.

            std::string RetrieveParameterString( char const *pszKey );
            bool        Decode( byte_t *pbyHeader );

            static bool ValidateHeader( byte_t *pbyHeader, int iHeaderSize );

            template <typename tParam> bool RetrieveParameter( char const *pszKey, tParam &rParam );

        private:

            //////////
            // Private attributes.

            enum EParamType
            {
                eParamTypeInt       = 0,
                eParamTypeDouble,
                eParamTypeString,
                eParamTypeList,
                eParamTypeHex,
                eParamTypeISODate
            };

            //////////
            //

            typedef struct tagParamItem
            {
                string      m_sKey;         // A local copy of the key.
                EParamType  m_eParamType;   // Hint to decoding the parameter string.
                string      m_sParameter;   // Parameter(s) as a string.
                string      m_sComment;     // Comment field.
            }
            ParamItem_t;

            //////////
            // Aliases.

            using ParamMap_t         = std::map<string, ParamItem_t>;
            using ParamMapIterator_t = ParamMap_t::iterator;

            //////////
            // Attributes.

            ParamMap_t  m_ParamMap;

            //////////
            //  Methods.

            void LoadParameterList( void );

            bool GetSingleParam( ParamMapIterator_t pItem, const byte_t *pbyHeader );

            ParamItem_t RetrieveParameter( string &rsKey );

            bool ParamsFromHeader( ParamMapIterator_t pItem, byte_t *pbyHeader );

            static string GetSingleParamOnly( const string &rsIdentifier, const byte_t *pbyHeader );
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Template public method needing to be externally availabe to other modules,
    // hence implemented in this header file.

    template <typename tParam>
    bool CVCRAFTParser::RetrieveParameter( char const *pszKey, tParam &rParam )
    {
        bool bConverted = false;

        string sKey{ pszKey };

        ParamItem_t Item = RetrieveParameter( sKey );

        switch ( Item.m_eParamType )
        {
            case eParamTypeInt:

                rParam = std::stol( Item.m_sParameter );
                bConverted = true;
                break;

            case eParamTypeDouble:

                rParam = std::stod( Item.m_sParameter );
                bConverted = true;
                break;

            case eParamTypeHex:

                rParam =  std::stoll( Item.m_sParameter, 0, 16 );
                bConverted = true;
                break;

            case eParamTypeList:
            case eParamTypeISODate:
            case eParamTypeString:
            default:

                // Nothing to do as these types are natively strings. Even the comma
                // seperated list should be treated as a string for now.
                bConverted = true;
                break;
        }

        return bConverted;
    }

}       // end namespace NCodec.

#endif // end VCRAFTPARSER_H_INCLUDED

