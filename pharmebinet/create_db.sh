#!/bin/bash

mkdir -p ~/git/pharmebinet/data

docker run --rm \
  --volume ~/git/pharmebinet/data:/data \
  --volume ~/git/pharmebinet/pharmebinet:/dumps \
  neo4j:5.26.4 \
  neo4j-admin database load --from-path=/dumps pharmebinet --overwrite-destination=true
