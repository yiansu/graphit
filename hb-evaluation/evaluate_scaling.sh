#!/bin/bash

# environment
ROOT_DIR=`git rev-parse --show-toplevel`
source /nfs-scratch/yso0488/heartbeatcompiler0/enable ;
source /project/extra/burnCPU/enable ;
source common.sh ;

# experiment
experiment=scaling
keyword=exectime
mkdir -p ${ROOT_DIR}/hb-evaluation/results/${experiment};

########################################################
# experiment sections
baseline=false
hbc_acc=false                       # software polling + acc
hbc_static=false                    # software polling + static chunksize
hbc_rf=false                        # rollforward + interrupt ping thread
hbc_rf_kmod=false                   # rollforward + kernel module
openmp=false                        # default openmp with different schedulers
openmp_nested_parallelism=false     # run openmp with nested parallelism enabled
openmp_chunksize_sensitivity=false  # run openmp with various chunksize

# benchmark targetted
benchmarks=(BFS CC CF PR PRDelta SSSP)

# openmp settings
omp_schedules=(STATIC DYNAMIC GUIDED)
omp_chunksizes=(0 1 2 4 8 16)  # 0 means default chunksize
########################################################

function run_and_collect {
  local technique=${1}
  local results_path=${2}
  mkdir -p ${results_path} ;
  local output=${results_path}/output.txt

  if [ ${technique} == "baseline" ] ; then
    # burnP6
    killall burnP6 &> /dev/null ;

    for coreID in `seq -s' ' 4 2 10` ; do
      taskset -c ${coreID} burnP6 &
    done

    for i in `seq 1 3` ; do
      taskset -c 2 make run_baseline >> ${output} ;
    done

    # kill burnP6
    killall burnP6 &> /dev/null ;

  elif [ ${technique} == "openmp" ] ; then
    for i in `seq 1 5` ; do
      WORKERS=28 \
      numactl --physcpubind=28-55 --interleave=all \
      make run_openmp >> ${output} ;
    done

  else
    for i in `seq 1 5` ; do
      WORKERS=28 \
      CPU_FREQUENCY_KHZ=${cpu_frequency_khz} \
      KAPPA_USECS=${heartbeat_interval} \
      numactl --physcpubind=28-55 --interleave=all \
      make run_hbc >> ${output} ;
    done

  fi

  collect ${results_path} ${output} ;
}

# main script
pushd . ;

# preparation
cd ${ROOT_DIR}/hb-benchmarks ;
make link &> /dev/null ;

# run experiment per benchmark
for benchmark in ${benchmarks[@]} ; do

  results=${ROOT_DIR}/hb-evaluation/results/${experiment}/${benchmark}
  mkdir -p ${results} ;

  cd ${benchmark} ;
  clean ;

  # baseline
  if [ ${baseline} = true ] ; then
    clean ; make baseline INPUT_SIZE=${input_size} &> /dev/null ;
    run_and_collect baseline ${results}/baseline ;
  fi

  # hbc_acc
  if [ ${hbc_acc} = true ] ; then
    clean ; make hbc INPUT_SIZE=${input_size} ACC=true CHUNKSIZE=1 &> /dev/null ;
    run_and_collect hbc_acc ${results}/hbc_acc ;
  fi

  # hbc_rf
  if [ ${hbc_rf} = true ] ; then
    clean ; make hbc INPUT_SIZE=${input_size} ENABLE_ROLLFORWARD=true CHUNK_LOOP_ITERATIONS=false &> /dev/null ;
    run_and_collect hbc_rf ${results}/hbc_rf ;
  fi

  # openmp
  if [ ${openmp} = true ] ; then
    for omp_schedule in ${omp_schedules[@]} ; do
      clean ; make openmp INPUT_SIZE=${input_size} OMP_SCHEDULE=${omp_schedule} &> /dev/null ;
      run_and_collect openmp ${results}/openmp_`echo -e ${omp_schedule} | tr '[:upper:]' '[:lower:]'` ;
    done
  fi

  # openmp_nested_parallelism
  if [ ${openmp_nested_parallelism} = true ] ; then
    for omp_schedule in ${omp_schedules[@]} ; do
      clean ; make openmp INPUT_SIZE=${input_size} OMP_SCHEDULE=${omp_schedule} OMP_NESTED_PARALLELISM=true &> /dev/null ;
      run_and_collect openmp ${results}/openmp_`echo -e ${omp_schedule} | tr '[:upper:]' '[:lower:]'`_nested_parallelism ;
    done
  fi

  # openmp_chunksize_sensitivity
  if [ ${openmp_chunksize_sensitivity} = true ] ; then
    for omp_schedule in ${omp_schedules[@]} ; do
      for omp_chunksize in ${omp_chunksizes[@]} ; do
        clean ; make openmp INPUT_SIZE=${input_size} OMP_SCHEDULE=${omp_schedule} OMP_CHUNKSIZE=${omp_chunksize} &> /dev/null ;
        run_and_collect openmp ${results}/openmp_`echo -e ${omp_schedule} | tr '[:upper:]' '[:lower:]'`_${omp_chunksize} ;
      done
    done
  fi

  clean ;
  cd ../ ;

done

popd ;
