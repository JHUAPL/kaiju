#!/bin/sh
# usage: correctOMPenvironment.sh NODEFILE prog <args>

# get my host and figure out how many other jobs are on this node
myHost=`hostname -f`  # Used on cheyenne
numJobsOnMyNode=`grep -c "^$myHost\$" $1`
if [ $numJobsOnMyNode -eq 0 ]
then
   myHost=`hostname`  # Used on pleiades.
   numJobsOnMyNode=`grep -c "^$myHost\$" $1`
fi
if [ $numJobsOnMyNode -eq 0 ]
then
   echo "Warning: this host not found in nodefile $1!"
else
   # split the CPUs on this node between me and the other jobs
   #echo "myHost = $myHost"
   echo "numJobsOnMyNode = $numJobsOnMyNode"
   let newNumThreads=$OMP_NUM_THREADS/$numJobsOnMyNode
   if [ $newNumThreads -lt 1 ]
   then
      newNumThreads=1
   fi
   export OMP_NUM_THREADS=$newNumThreads
   #echo "Setting OMP_NUM_THREADS = $OMP_NUM_THREADS"
fi

# call the arguments to this script
$2 "${@:3}"

