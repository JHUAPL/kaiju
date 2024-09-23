#!/usr/bin/env bash

# usage: pinCpuCores.sh <program executable> <program arguments>

# THIS IS ONLY MEANT TO WORK ON DERECHO
# USE ON OTHER SYSTEMS AT YOUR OWN RISK

if [[ ! -v PMI_LOCAL_SIZE ]] || [[ ! -v PMI_LOCAL_RANK ]]; then
    echo "pinCpuCores.sh is only meant to be used on systems that set the"
    echo "   PMI_LOCAL_SIZE and PMI_LOCAL_RANK environment variables."
    echo "   One or both of these was not set. Please use the"
    echo "   correctOMPEnvironment.sh script in conjunction with omplace instead"
    exit 1
fi

let numCpus=`lscpu | sed --quiet "s/^CPU(s): \\+\\([0-9]\\+\\)$/\\1/p"`
let threadsPerCore=`lscpu | sed --quiet "s/^Thread(s) per core: \\+\\([0-9]\\+\\)$/\\1/p"`
let numCores=$numCpus/$threadsPerCore
let newNumThreads=$numCores/$PMI_LOCAL_SIZE
if [[ $newNumThreads -lt 1 ]]; then
   newNumThreads=1
fi
let minThread=$newNumThreads*$PMI_LOCAL_RANK
let maxThread=$minThread+$newNumThreads-1
export OMP_NUM_THREADS=$newNumThreads

#echo "I am rank $PMI_LOCAL_RANK of $PMI_LOCAL_SIZE on this node using cpus $minThread to $maxThread"

# call the arguments to this script
taskset --cpu-list $minThread-$maxThread $1 "${@:2}"

