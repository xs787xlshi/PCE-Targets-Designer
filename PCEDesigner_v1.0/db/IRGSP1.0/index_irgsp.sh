# index the rice IRGSP1.0 using Batmis
/usr/local/bin/build_index IRGSP1.0.fa   
# index the rice IRGSP1.0 using ncbi-blast
/hellogene/soft/ncbi-blast-2.3.0+/bin/makeblastdb -in IRGSP1.0.fa -dbtype nucl -parse_seqids -out IRGSP1.0.fa
