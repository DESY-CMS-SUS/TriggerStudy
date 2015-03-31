#!/bin/zsh
## make sure the right shell will be used
#$ -S /bin/zsh
## Job name
#$ -N HLTjob
## the cpu time for this job
#$ -l h_rt=02:59:00
## the maximum memory usage of this job
#$ -l h_vmem=3900M
## operating system
#$ -l distro=sld6
## architecture
#$ -l arch=amd64
## environment and cwd
#$ -V
#$ -cwd
## stderr and stdout are merged together to stdout
#$ -j y
##(send mail on job's end and abort)
##$ -m a
#$ -l site=hh
## define outputdir,executable,config file and LD_LIBRARY_PATH
#$ -v EXECUTABLE=hlt_template_file.py

export X509_USER_PROXY=~/.globus/my.proxy

echo job start at `date`
echo "Running job on machine " `uname -a`

# expect to be in basedir already
BASEDIR=`pwd -P`
echo "Locating in "$BASEDIR

eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
CMSRUN=`which cmsRun`
echo "CMSSW environnment: $CMSSW_VERSION"

TaskID=$((SGE_TASK_ID))
echo "SGE_TASK_ID: " $TaskID

# list containing all file names
FileList=$1

InFile=$(cat $FileList | sed ''$TaskID'q;d')
echo "Going to process file:"
echo $InFile

OutDir=./$2
JobDir=$OutDir"/HLTtuple_chunk"$TaskID
OutFile="HLTtuple_chunk"$TaskID".root"

echo "Changing to workdir" $JobDir
if [ ! -d $JobDir ]; then
    mkdir -p $JobDir
fi

cd $JobDir

if [ -f processing ] && [ ! -f failed ] ; then
    echo "Already processing!"
    echo "Aborting."
    exit 1
fi

if [ -f processed ]; then
    echo "Already processed!"
    echo "Aborting."
    exit 1
fi

HLTconfig="cfg_HLT_chunk"$TaskID".py"

cp $BASEDIR/$EXECUTABLE $HLTconfig

sed -i "s|hlt_tuple.root|$OutFile|" $HLTconfig
sed -i "s|INFILE|$InFile|" $HLTconfig

echo "Starting processing"
touch processing

memtime=/usr/bin/time
$memtime -v time $CMSRUN $HLTconfig >>cmsRun.log 2>&1
rm processing

if [ -f $OutFile ]; then
    echo "Sucessfully processed!"
    touch processed
else
    echo "Failed processing!"
    touch failed
fi

echo "Complete at" `date`
