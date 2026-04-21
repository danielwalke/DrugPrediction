docker run --detach \
  --name neo4j-pharmebinet \
  --publish 7470:7474 \
  --publish 7683:7687 \
  --volume ~/git/pharmebinet/data:/data \
  --env NEO4J_AUTH=neo4j/password \
  --env NEO4J_initial_dbms_default__database=pharmebinet \
  --env NEO4J_PLUGINS='["apoc", "graph-data-science"]' \
  neo4j:5.26.4
