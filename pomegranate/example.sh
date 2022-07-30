python \
load.py ../datasets/human_known_kinases.txt ../structures/human saved_graphs/g-human-allknown-20A-0.0-3 \
 --debug --radius 20 --rsa 0.0 --format=sg \
 --nf=pdist,bfac,rsa,m1,m2,m3,m4,m5,m6,m7 \
 -N 5 --skip-download