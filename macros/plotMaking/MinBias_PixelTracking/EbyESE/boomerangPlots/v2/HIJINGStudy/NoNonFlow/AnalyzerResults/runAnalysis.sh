echo "Job started on $(date)"
root -l -b <<EOF
.x makeVNDet.C++
EOF
root -l -b <<EOF
.x ReadTree_normDet.C++
EOF
cd DDResp
root -l -b <<EOF
.x makeDDResp.C++
EOF
cd ../
bash unfold.sh
echo "Analysis complete!"
echo "Job ended on $(date)"