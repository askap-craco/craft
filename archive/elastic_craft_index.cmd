curl -XPUT 'akingest01.atnf.csiro.au:9200/craftscans?pretty' -d 'Content-Type:application/json' -d'
{
 "mappings": {
   "antscan": {
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
