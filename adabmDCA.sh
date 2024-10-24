#!/bin/bash

# Update soft link
{
  rm adabmDCA &&
  ln -s ./src/adabmDCA adabmDCA 
} || {
  ln -s ./src/adabmDCA adabmDCA 
}

# Check if the first positional argument is provided
if [ -z "$1" ]; then
  echo "Error: No command provided. Use 'train', 'sample', 'contacts', 'energies', or 'DMS'."
  exit 1
fi

# Assign the first positional argument to a variable
COMMAND=$1
shift # Remove the first positional argument, so "$@" now contains only the optional arguments

# Define short and long options
SHORT=d:,o:,m:,w:,p:,c:,l:,h
LONG=data:,output:,model:,weights:,path_params:,path_chains:,label:,alphabet:,lr:,nsweeps:,sampler:,nchains:,target:,nepochs:,pseduocount:,seed:,gsteps:,factivate:,density:,drate:,ngen:,nmeasure:,nmix:,beta:,help

# Parse command line arguments
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")

# Default values
outfolder="DCA_model"
model="bmDCA"
abc="protein"
sampler="gibbs"
lrate="0.01"
nsweeps="10"
nchains="10000"
pearson="0.95"
niter="50000"
seed="0"
gsteps="10"
factivate="0.001"
sparsity="0.98"
drate="0.01"
beta="1.0"
nmix="2"
abc="protein"
# Optional 
weights=""
params=""
chains=""
label=""
lrate=""
pcount=""
ngen=""


# Evaluate parsed options
if [ "$OPTS" ]; then
  # Extract input
  while [ $# -gt 0 ]; do
    case "$1" in
      -d | --data ) msa="$2"; shift 2 ;;
      -o | --output ) outfolder="$2"; shift 2 ;;
      -m | --model ) model="$2"; shift 2 ;;
      -w | --weights ) weights="-w $2"; shift 2 ;;
      -p | --path_params ) params="-p $2"; shift 2 ;;
      -c | --path_chains ) chains="-B $2"; shift 2 ;;
      -l | --label ) label="-k $2"; shift 2 ;;
      --alphabet ) abc="$2"; shift 2 ;;
      --lr ) lrate="-u $2 -v $2"; shift 2 ;;
      --nsweeps ) nsweeps="$2"; shift 2 ;;
      --sampler ) sampler="$2"; shift 2 ;;
      --nchains ) nchains="$2"; shift 2 ;;
      --target ) pearson="$2"; shift 2 ;;
      --nepochs ) niter="$2"; shift 2 ;;
      --pseudocount ) pcount="-d $2"; shift 2;;
      --seed ) seed="$2"; shift 2;;
      --gsteps ) gsteps="$2"; shift 2;;
      --factivate ) factivate="$2"; shift 2;;
      --density ) sparsity="$(echo 1.0 - $2 | bc -l)"; shift 2;;
      --drate ) drate="$2"; shift 2;;
      --beta ) beta="$2"; shift 2;;
      --ngen ) ngen="-s $2"; shift 2;;
      --nmix ) nmix="$2"; shift 2;;
      -h | --help ) help=true ;;
      -- ) shift; break ;;
      * ) echo "Unexpected option: $1"; exit 1 ;;
    esac
  done
fi

# Deal with alphabet and sampler
case "$sampler" in
  gibbs ) mc="-G";;
  metropolis ) mc="-M";;
  * ) echo "Unexpected sampler: $sampler"; exit 1 ;;
esac

case "$abc" in
  protein ) alpha="-b a";;
  rna ) alpha="-b n";;
  dna ) alpha="-b -ACGT";;
  * ) alpha="-b $abc";;
esac

# Perform actions based on parsed options
if [ "$help" ]; then
  echo "Usage:"
  echo "./adabmDCA.sh train -m [ model ] [ -d | --data ] [ -o | --output ] [ -h | --help ]"
  echo "./adabmDCA.sh sample [ -p | --path_params ] [ -d | --data ] [ -o | --output ] [ --ngen ] [ -h | --help ]"
  echo "./adabmDCA.sh energies [ -p | --path_params ] [ -d | --data ] [ -o | --output ] [ -h | --help ]"
  echo "./adabmDCA.sh DMS [ -p | --path_params ] [ -d | --data ] [ -o | --output ] [ -h | --help ]"
  echo "./adabmDCA.sh contacts [ -p | --path_params ] [ -o | --output ] [ -h | --help ]"
  echo "See https://github.com/spqb/adabmDCA or type ./adabmDCA -h for this C/C++ implementation"
  exit 0
fi

# Run different command lines according to model 
train=false
case "$COMMAND" in
  train ) train=true;;
  sample ) SCRIPT="./adabmDCA $params $label -f $msa -i 0 -S -L -a $outfolder $ngen -W $nmix $mc -J $beta -H $beta";;
  energies ) SCRIPT="./adabmDCA --energies $label -a $outfolder $alpha $params -f $msa";;
  DMS ) SCRIPT="./adabmDCA --dms -a $outfolder $alpha $params -f $msa";;
  contacts ) echo "Check the file $outfolder/"LABEL_frobenius.dat"";;
esac

if "$train"
then
  case "$model" in
    bmDCA ) SCRIPT="./adabmDCA -f $msa -z 1 -m 50 -L -a $outfolder $label $mc $weights $params $chains $alpha $lrate -e $nsweeps -s $nchains -c $pearson -i $niter $pcount -y $seed";;
    edDCA ) SCRIPT="./adabmDCA -f $msa -z 1 -m 50 -L -a $outfolder $label $mc $weights $params $chains $alpha $lrate -e $nsweeps -s $nchains -c $pearson -i $niter $pcount -y $seed -x $sparsity -V $drate";;
    eaDCA ) SCRIPT="./adabmDCA -f $msa -z 1 -m 50 -L -a $outfolder $label $mc $weights $chains $alpha $lrate -e $nsweeps -s $nchains -c $pearson -i $niter $pcount -y $seed -I 0. -Z -X $gsteps -U $factivate";;
  esac
fi

$SCRIPT


