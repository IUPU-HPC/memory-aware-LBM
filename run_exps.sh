#!/bin/bash
WORK_DIR=`pwd`
RESULT_DIR=${WORK_DIR}/results/`date +"%y%m%d_%H%M"`
mkdir -pv ${RESULT_DIR}

BUILD_DIR=${WORK_DIR}/build/bin
cd ${RESULT_DIR}


#dimension
# smallest basic test
# lx=8
# ly=8
# export STOP=10 #time stop
# thread_block=2 #each thread will compute this number of x lines or y lines
# blk_size=2

# lx=16
# ly=16
# export STOP=10 #time stop
# thread_block=4 #each thread will compute this number of x lines or y lines
# export OMP_NUM_THREADS=4
# blk_size=2

# lx=64
# ly=64
# export STOP=10 #time stop
# thread_block=16 #each thread will compute this number of x lines or y lines
# export OMP_NUM_THREADS=4
# blk_size=4

# verifiction test
lx=256
ly=256
export STOP=100 #time stop
thread_block=64 #each thread will compute this number of x lines or y lines
blk_size=32
export OMP_NUM_THREADS=4
export KMP_AFFINITY=verbose, granularity=fine, compact
#blk_size=(1 2 4 8 16 32 64 128 256 512)
##for CASE_NAME in tight_block_papi
#for CASE_NAME in tight_block
#do
    #for ((k=0; k<${#blk_size[@]}; k++))
    #do
       #for i in `seq 3`
        #do
            #${BUILD_DIR}/unsteady_${CASE_NAME} $lx $ly ${blk_size[k]}

            #echo "No.$i exp of case ${CASE_NAME} done "
        #done
    #done
#done


# export STOP=5000 #time stop
# thread_block=16 #each thread will compute this number of x lines or y lines
# blk_size=4
# export OMP_NUM_THREADS=16

# for CASE_NAME in  origin combine quicktest tight tight_block segment
# for CASE_NAME in origin combine tight

# for CASE_NAME in origin tight origin_openmp combine_openmp tight_openmp tight_block_openmp panel panel_openmp
# for CASE_NAME in origin_openmp_papi tight_openmp_papi tight_block_openmp_papi
# for CASE_NAME in origin combine tight origin_openmp combine_openmp tight_openmp origin_papi origin_openmp_papi combine_papi combine_openmp_papi tight_papi tight_openmp_papi tight_block tight_block_openmp tight_block_papi tight_block_openmp_papi panel panel_openmp panel_papi panel_openmp_papi
for CASE_NAME in origin_papi tight_block_openmp_papi
do
    echo "result will be saved in ${RESULT_DIR}"
    for i in `seq 1`
    do
        # ${BUILD_DIR}/unsteady $lx $ly $blk_size &> ${RESULT_DIR}/${CASE_NAME}_${i}.log
        ${BUILD_DIR}/unsteady_${CASE_NAME} $lx $ly $blk_size $thread_block &>> ${RESULT_DIR}/log
        echo "No.$i exp of case ${CASE_NAME} done " &>> ${RESULT_DIR}/log
    done
done

# blk_size=32
# chunk_size=32
# for CASE_NAME in segment segment_openmp
# do
#     echo "result will be saved in ${RESULT_DIR}"
#     for i in `seq 1`
#     do
#         # ${BUILD_DIR}/unsteady $lx $ly $blk_size &> ${RESULT_DIR}/${CASE_NAME}_${i}.log
#         ${BUILD_DIR}/unsteady_${CASE_NAME} $lx $ly $blk_size $chunk_size &>> ${RESULT_DIR}/log
#         echo "No.$i exp of case ${CASE_NAME} done " &>> ${RESULT_DIR}/log
#     done
# done

#run tight_block with different block size
# blk_size=(32)
# for CASE_NAME in tight_block tight_block_openmp
# do
#     echo "result will be saved in ${RESULT_DIR}"
#     for ((k=0; k<${#blk_size[@]}; k++)); do
#         for i in `seq 1`
#         do
#             ${BUILD_DIR}/unsteady_${CASE_NAME} $lx $ly ${blk_size[k]} &>> ${RESULT_DIR}/log
#             echo "No.$i exp of case ${CASE_NAME} ${blk_size[k]} done " &>> ${RESULT_DIR}/log
#         done
#     done
# done

# #run tight_block with different segment size
# blk_size=(32)
# for CASE_NAME in segment segment_openmp
# do
#     echo "result will be saved in ${RESULT_DIR}"
#     for ((k=0; k<${#blk_size[@]}; k++)); do
#         for i in `seq 1`
#         do
#             ${BUILD_DIR}/unsteady_${CASE_NAME} $lx $ly ${blk_size[k]} &>> ${RESULT_DIR}/log
#             echo "No.$i exp of case ${CASE_NAME} ${blk_size[k]} done " &>> ${RESULT_DIR}/log
#         done
#     done
# done
