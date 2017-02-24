echo "Job started on $(date)"
#rm *.root
#root -l -b <<EOF
#.x makeVNDet.C++
#EOF
root -l -b <<EOF
.x ReadTree_normDet.C++
EOF
cd DDResp
root -l -b <<EOF
.x makeDDResp.C++
EOF
echo "Analysis complete!"
echo "Job ended on $(date)"