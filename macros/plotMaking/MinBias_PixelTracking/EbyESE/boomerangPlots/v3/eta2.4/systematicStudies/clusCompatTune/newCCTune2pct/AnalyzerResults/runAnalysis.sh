N="$1"
E="$2"
PMN="$3"
PMX="$4"
T="$5"

root -l -b <<EOF
.L makeVNDet.C++ 
makeVNDet($N, $E, $PMN, $PMX, $T)
EOF
root -l -b <<EOF
.L ReadTree_normDet.C++
ReadTree_normDet($N, $E, $PMN, $PMX, $T)
EOF
cd DDResp
root -l -b <<EOF
.L makeDDResp.C++
makeDDResp($N)
EOF