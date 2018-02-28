curl -XDELETE 'localhost:9200/craftscans?pretty'
curl -XPUT 'localhost:9200/craftscans?pretty' -H 'Content-Type:application/json' -d'
{
 "mappings": {
   "antscan": {
   "dynamic_templates": [
            { "notanalyzed": {
                  "match":              "*", 
                  "match_mapping_type": "string",
                  "mapping": {
                      "type":        "string",
                      "index":       "not_analyzed"
                  }
               }
            }
          ],

     "properties": {
       "ant_direction": {
         "type": "geo_point"
       },
       "beam_directions": {
         "type": "geo_point"
       },
       "scan_start": {
         "type": "date"
       },
       "hdr_utc": {
         "type":"date"
       },
       "target": {
          "type":"string",
          "index":"not_analyzed"
       },
       "field_name": {
           "type":"string",
	   "index":"not_analyzed"
       },
       "scan_intent": {
           "type":"string",
	   "index":"not_analyzed"
       },
       "antname": {
              "type":"string",
       	      "index":"not_analyzed"
        },
       "sb_owner": {
              "type":"string",
       	      "index":"not_analyzed"
        },
       "sb_alias": {
              "type":"string",
       	      "index":"not_analyzed"
        },
       "sb_template": {
              "type":"string",
       	      "index":"not_analyzed"
        },
       "filterbanks": {
         "type":"nested",
         "properties": {
           "atime": {
             "type":"date"
           },
           "mtime":{
             "type":"date"
           },
           "ctime:":{
             "type":"date"
           }
         }

       }
     }
   }
 }
}'