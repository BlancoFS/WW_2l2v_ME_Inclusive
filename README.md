# Analysis of WW Inclusive production 

First of all, let's start building the programming enviroment. The next part explain the instructions to execute the code. 

First, log in gridui:

```ssh -Y ***@gridui.ifca.es -o ServerAliveInterval=240```


Set CMS enviroment for the first time:

```
bash -l

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc820

cmsrel CMSSW_10_6_10
```


Enter CMS eviroment on gridui:

```
bash -l

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd CMSSW_10_6_10/src

cmsenv
```

Compile code:

```scram b -j 8```


Get some code from gitHub:

```git clone https://github.com/BlancoFS/...```


Change condor scehduler if it's not running properly: 

```
export _condor_SCHEDD_HOST="bigbird02.cern.ch"
```

## Install MoMEMta framework to compute the Matrix Element

Complete information in the next tutorial page:

```
https://github.com/BlancoFS/Tutorials/tree/main/MoMEMta
```



## Run code:

It should be done a Run.sh script to submit jobs in Slurm, the old way is for condor.


```
sbatch -o logfile.log -e errofile.err --qos=gridui_sort --partition=cloudcms Run.sh
```

Run and make plots with LatinoAnalysis tools:

```
mkShapesMulti.py --pycfg=configuration.py --doBatch=1 --batchSplit=Samples,Files --batchQueue=testmatch

mkShapesMulti.py --pycfg=configuration.py --doHadd=1 --batchSplit=Samples,Files --doNotCleanup --nThreads=8

mkPlot.py --pycfg=configuration.py --inputFile=rootFile/plots_WW_2016.root --minLogC=0.01 --minLogCratio=0.01 --maxLogC=1000 --maxLogCratio=1000 --showIntegralLegend=1

Plot_Sig_Bkg.py --pycfg=configuration.py --inputFile=rootFile/plots_WW_2016.root

```
