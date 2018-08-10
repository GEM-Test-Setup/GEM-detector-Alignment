//Daniel DeLayo, Marisa Petrusky, Tao

1) Convert Gem data into histogram format (1 hist in X, 1 hist in Y per GEM)
2) Run Data through convertRaw.c with no offsets
3) Run geometrically incorrect tracks through makeTruthTrack.c (may need format conversion)
4) Run Data and Offsets through convertRaw.c

To generate test data,

1) Run makeTracks.c
