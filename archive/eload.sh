#!/bin/bash

for f in $@ ; do
    echo $f
    #curl -s -H "Content-Type: application/x-ndjson" -XPOST localhost:9200/_bulk --data-binary "@${f}"
    wget --header "Content-Type: application/x-ndjson"  http://localhost:9200/_bulk --post-file $f -O $f.eloadlog
done
