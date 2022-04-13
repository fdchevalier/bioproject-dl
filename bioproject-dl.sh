#!/bin/bash
# Title: bioproject-dl.sh
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2022-04-09
# Modified in: 2022-04-11
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Download fastq files from a given BioProject."



#==========#
# Versions #
#==========#

# v0.1 - 2022-04-11: add parallel downloads
# v0.0 - 2022-04-09: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -b|--bp value -d|--dir path -p|--pl number -s|--skip -h|--help

Aim: $aim

Version: $version

Options:
    -b, --bp        BioProject number
    -d, --dir       target directory
    -p, --pl        number of parallel SRA downloads [default: 10]
    -s, --skip      skip download if fastq file already present
    -h, --help      this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    [[ -n $2 ]] && exit $2
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    [[ $? != 0 ]] && error "Command $1 is needed. Exiting..." 1
}


# Progress bar
## Usage: ProgressBar $mystep $myend
function ProgressBar {
    if [[ -t 1 ]]
    then
        # Process data
        let _progress=(${1}*100/${2}*100)/100
        let _done=(${_progress}*4)/10
        let _left=40-$_done
        # Build progressbar string lengths
        _fill=$(printf "%${_done}s")
        _empty=$(printf "%${_left}s")

        # Build progressbar strings and print the ProgressBar line
        # Output example:
        # Progress : [########################################] 100%
        #printf "\rProgress : [${_fill// /=}${_empty// / }] ${_progress}%%"
        printf "\r\e[32mProgress:\e[00m [${_fill// /=}${_empty// / }] ${_progress}%%"
        
        [[ ${_progress} == 100 ]] && echo ""
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 0
}

function clean_up_kill {
    echo
    [[ -n "$run" ]] && rm -rf runinfo "$ldir/$fln/$run"*
    [[ -n $myjobs ]] && kill ${myjobs[@]} 2> /dev/null
    exit 1
}


# Prefetch
function myprefetch {
    #Download
    retry=0
    while [[ $retry -lt 2 ]]
    do
        # Download SRA file
        prefetch -X 1T -C yes -r yes -O "$1/$2/" $3 &>> "$4"
        [[ $(find "$1/$2/$3" 2> /dev/null) ]] && break || ((retry++))  #|| break
    done

    # If max download trials reached, issue message and move to the next
    if [[ $retry -eq 2 ]]
    then
        echo "$run: dowloading error" >> "$4"
        return 1
    fi

    # Convert sra into fastq
    fasterq-dump -O "$1/$2/" -f --split-files "$1/$2/$3" &>> "$4"
    rm -R "$1/$2/$3"

    # Rename file with more meaningful name
    if [[ $(find "$1/$2/" -name "*.fastq" | wc -l) -eq 2 ]]
    then
        mv "$1/$2/${3}_1.fastq" "$1/$2/${2}_R1.fastq"
        mv "$1/$2/${3}_2.fastq" "$1/$2/${2}_R2.fastq"
    else
        mv "$1/$2/${3}.fastq" "$1/$2/${2}_R1.fastq"
    fi

    # Compress
    find "$1/$2/" -name "*.fastq" -exec pigz {} +
}



#==============#
# Dependencies #
#==============#

test_dep wget 
test_dep prefetch
test_dep fasterq-dump



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -b|--bp     ) bioproject="$2" ; shift 2 ;;
        -d|--dir    ) ldir="$2"       ; shift 2 ;;
        -p|--pl     ) pl="$2"         ; shift 2 ;;
        -s|--skip   ) skip="skip"     ; shift   ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
[[ -z "$bioproject" ]] && error "BioProject is required. Exiting...\n$(usage)" 1

# Default values
[[ -z "$ldir" ]] && ldir="."
[[ -z $pl ]] && pl=10

log="$ldir/log"



#============#
# Processing #
#============#

# Handling status and exit
set -o pipefail
trap 'clean_up_kill' SIGINT SIGTERM
trap 'clean_up runinfo $log' EXIT

# Download related information to data project
runinfo=$(wget -q -O - "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=${bioproject}")
[[ $(sed "/^$/d" <<< "$runinfo" | wc -l) -eq 0 ]] && error "Bioproject does not seem to exists. Exiting..." 0

# Field of interest (library name and weblink)
fdn=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "SampleName" | cut -d ":" -f 1)
fdr=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "Run" | cut -d ":" -f 1)
flk=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "download_path" | cut -d ":" -f 1)

# Download fastq files
myend=$(tail -n +2 <<< "$runinfo" | sed "/^$/d" | wc -l)
while read line
do
    # Counter
    ((mystep++))

    # Filename, run and download link
    fln=$(cut -d "," -f $fdn <<<$line)
    run=$(cut -d "," -f $fdr <<<$line)
    
    # Folder
    echo "$fln" >> "$log"
    [[ ! -d "$ldir/$fln/" ]] && mkdir -p "$ldir/$fln/"
    
    # Check for fastq if requested
    if [[ "$skip" ]]
    then
        ls_fl=$(find "$ldir/$fln/" -name "${fln}_R*.fastq*" | wc -l)
        [[ $ls_fl -gt 0 ]] && continue
    fi
    
    # Pause when number of running job is expected
    myjobs=($(jobs -p))
    while [[ ${#myjobs[@]} -ge $pl ]]
    do
        sleep 30s
        myjobs=($(jobs -p))
    done

    # Fetch SRA file in the background (except if last)
    if [[ $mystep -lt $myend ]]
    then
        (myprefetch "$ldir" "$fln" "$run" "$log") &
    else
        myprefetch "$ldir" "$fln" "$run" "$log"
    fi

    # Progress bar
    ProgressBar $mystep $myend

done < <(tail -n +2 <<< "$runinfo" | sed "/^$/d")

# Wait until every background job is done
wait $(jobs -p)

# Check for errors
[[ $(grep -E -i "error|warning" "$log") ]] && warning "Errors or waning present in $log." && exit 1

exit 0
