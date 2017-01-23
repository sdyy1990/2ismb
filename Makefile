all: Builder Query


GITVERSION= "no-version"

CXXFLAGS+= -std=gnu++11 -DGITVERSION=\"$(GITVERSION)\"
CXXFLAGS+= -Wno-format -Wno-pointer-arith  -I ./

CXXFLAGS+=   -march=native   -O3


HEADERS= ../*.h

cluster: cluster.cc

reader: reader.cc

