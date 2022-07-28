#! /usr/bin/env bash
# Get list of all kinase types for H. Sapiens
cat ../datasets/phosphoELM_all_2015-04.dump |\
grep -i "homo sapiens" |\
cut -f1,6,8 | cut -f2 | sed 's/^$/UNKNOWN/' | sort | uniq

# Get all phosphosites with known kinase in humans
cat ../datasets/phosphoELM_all_2015-04.dump |\
grep -i "homo sapiens" |\
cut -f1,6,8 | grep -i "\S\sHomo sapiens" |\
#cut -f2 | sort | uniq | wc -l 


# WITH USEFUL FIELDS
cat ../datasets/phosphoELM_all_2015-04.dump |\
grep -i "homo sapiens" |\
cut -f1,3,4,6,8 | grep -i "\S\sHomo sapiens" > human_known_kin
