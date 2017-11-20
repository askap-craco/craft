PUT craftscans
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
}
