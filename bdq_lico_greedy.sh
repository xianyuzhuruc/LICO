#!/bin/bash
result_dir="/mnt/supplement/result"
# Modify it to the dir where your dataset is located, the structure is as followed:
# data_dir/CW12B/
#                               |_______CW12B.docs
#                               |_______CW12B.freqs  // the file containing frequency of each term is not necessary for LICO project, but essential for d2si project
#                               |_______CW12B-2.queries
#                               |_______CW12B-3.queries
#                               |_______CW12B-4.queries
#                               |_______CW12B-5.queries  // this file contains queries that the number of terms is equal or larger than 5
data_dir="/mnt/home/datasets"
code_dir="/mnt/home/codes/LICO-last/build"


source_dir="./"
index_type="lico"
decode_type="all"

#
## build index
rm -r $code_dir
cmake -B $code_dir -S $source_dir
sleep 1
cd $code_dir
make -j


compress_type="none"
epsilon=0
# add your datasets like this
for dataset in "CW12B" "CCNews" "WITD" 
do
  for lico_m in 2048
  do
    read_only="f"

    index_save_dir="$result_dir/index/$dataset/$index_type-$compress_type-$lico_m/$dataset-$index_type-$epsilon/"
    log_save_dir="$result_dir/log/$dataset/$index_type-$compress_type-$lico_m"
    mkdir -p $index_save_dir
    mkdir -p $log_save_dir

  echo "=================Building index for dataset: $dataset with epsilon: $epsilon================="
  $code_dir/lico_build_none $index_type $data_dir/$dataset/$dataset $index_save_dir $epsilon $read_only $decode_type $log_save_dir/$dataset-$index_type-$epsilon $lico_m $compress_type

   echo "=================Decode index for dataset: $dataset with epsilon: $epsilon================="
   for decode_type in "normal" "simd"
   do
     for repeat in 1 2 3 4 5
     do
       echo "————————————repeat : $repeat————————————"
       $code_dir/lico_decode_none $index_type $data_dir/$dataset/$dataset $index_save_dir $epsilon $decode_type $log_save_dir $compress_type
     done
     echo " "
   done

   echo "=================Query index for dataset: $dataset with epsilon: $epsilon================="
   for query_type in "AND" "OR"
   do
     for query_num in 5 4 3 2
     do
       for decode_type in "simd"
       do
         echo "=========================dataset : $dataset epsilon: $epsilon query_type: $query_type query_num: $query_num decode_type: $decode_type========================="
         $code_dir/lico_query_none $index_type $data_dir/$dataset/$dataset $index_save_dir $epsilon $decode_type $data_dir/$dataset/$dataset-$query_num.queries $query_type $log_save_dir/$dataset-$compress_type-$epsilon.query-$query_num-log.txt $compress_type
         echo " "
       done
     done
   done
  done
done