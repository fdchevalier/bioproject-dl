#!/bin/bash
# Title: bioproject-dl.sh
# Version: 0.6
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2022-04-09
# Modified in: 2025-01-17
# License: GPL v3



#======#
# Aims #
#======#

aim="Download fastq files from a given BioProject."



#==========#
# Versions #
#==========#

# v0.6 - 2025-01-17: check runinfo header
# v0.5 - 2023-09-07: remove samples which may be erroneously associated with BioProject
# v0.4 - 2022-09-12: add merge option / create rename function / improve traps
# v0.3 - 2022-09-11: replace wget with esearch/efetch to download runinfo
# v0.2 - 2022-04-13: add option to load runinfo file
# v0.1 - 2022-04-11: add parallel downloads / improve trap / improve runinfo handling
# v0.0 - 2022-04-09: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -b|--bp value -d|--dir path --r|--ri path -p|--pl number -m|--merge -s|--skip -h|--help

Aim: $aim

Version: $version

Options:
    -b, --bp        BioProject number
    -d, --dir       target directory
    -r, --ri        provide a runinfo file instead of the automatic download
    -p, --pl        number of parallel SRA downloads [default: 10]
    -m, --merge     merge fastq files belonging to the same sample
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

        # Header
        [[ -z "$3" ]] && hdr="Progress" || hdr="$3"

        # Build progressbar strings and print the ProgressBar line
        # Output example:
        # Progress : [########################################] 100%
        #printf "\rProgress : [${_fill// /=}${_empty// / }] ${_progress}%%"
        printf "\r\e[32m${hdr}:\e[00m [${_fill// /=}${_empty// / }] ${_progress}%%"

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
    # Stop jobs
    [[ -n $myjobs ]] && kill ${myjobs[@]} 2> /dev/null
    wait

    # Remove files
    a=("$@")
    for i in ${a[@]}
    do
        rm -rf "$i"*
    done

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
}


# Rename fastq
function rename_fastq {
    # Rename if only one read
    if [[ $(find "$1/$2/" \( -name "*.fastq" ! -name "*_1.fastq" ! -name "*_2.fastq" \) | wc -l) -eq 1 ]]
    then
        mv "$1/$2/"*.fastq "$1/$2/${2}_R1.fastq"
    elif [[ $(find "$1/$2/" \( -name "*.fastq" ! -name "*_1.fastq" ! -name "*_2.fastq" \) | wc -l) -gt 1 ]]
    then
        cat "$1/$2/"*.fastq > "$1/$2/${2}_R1.fastq" && find "$1/$2/" \( -name "*.fastq" ! -name "*_R1.fastq" \) -exec rm {} +
    fi

    # Rename read 1 file with more meaningful name
    if [[ $(find "$1/$2/" -name "*_1.fastq" | wc -l) -gt 1 ]]
    then
        cat "$1/$2/"*_1.fastq > "$1/$2/${2}_R1.fastq" && rm "$1/$2/"*_1.fastq
    elif [[ $(find "$1/$2/" -name "*_1.fastq" | wc -l) -eq 1 ]]
    then
        mv "$1/$2/"*_1.fastq "$1/$2/${2}_R1.fastq"
    fi

    # Rename read 2 file with more meaningful name
    if [[ $(find "$1/$2/" -name "*_2.fastq" | wc -l) -gt 1 ]]
    then
        cat "$1/$2/"*_2.fastq > "$1/$2/${2}_R2.fastq" && rm "$1/$2/"*_2.fastq
    elif [[ $(find "$1/$2/" -name "*_2.fastq" | wc -l) -eq 1 ]]
    then
        mv "$1/$2/"*_2.fastq "$1/$2/${2}_R2.fastq"
    fi

    # Compress
    find "$1/$2/" -name "*.fastq" -exec pigz {} +
}



#==============#
# Dependencies #
#==============#

test_dep esearch
test_dep efetch
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
        -r|--ri     ) runinfo="$2"    ; shift 2 ;;
        -p|--pl     ) pl="$2"         ; shift 2 ;;
        -m|--merge  ) merge="merge"   ; shift   ;;
        -s|--skip   ) skip="skip"     ; shift   ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
[[ -z "$bioproject" && -z "$runinfo" ]] && error "BioProject is required. Exiting...\n$(usage)" 1

# Check runinfo
if [[ -n "$runinfo" ]]
then
    [[ ! -e "$runinfo" ]] && error "runinfo file does not exist. Exiting..." 1
    runinfo=$(sed "/^$/d" "$runinfo")
    [[ -z "$runinfo" ]] && error "runinfo file is empty. Exiting..." 1
    [[ ! $(head -1 <<< "$runinfo" | grep -w -n "SampleName") ]] && error "header cannot be detected from the runinfo. Exiting..." 1
fi

# Check folder
[[ -n "$ldir" && -f "$ldir" ]] && error "Cannot create folder $ldir. Exiting..." 1
[[ -n "$ldir" && ! -d "$ldir" ]] && mkdir -p "$ldir"


# Default values
[[ -z "$ldir" ]] && ldir="."
[[ -z $pl ]] && pl=10

log="$ldir/log"



#============#
# Processing #
#============#

# Handling status and exit
set -o pipefail
trap 'clean_up_kill ${to_del[@]}' SIGINT SIGTERM
trap '[[ $? -eq 0 ]] && clean_up $log' EXIT

#---------#
# Runinfo #
#---------#

# Download related information to data project
[[ -z "$runinfo" ]] && runinfo=$(esearch -db sra -query ${bioproject} | efetch -format runinfo | sed "/^$/d")
[[ -z "$runinfo" ]] && error "Bioproject does not seem to exist. Exiting..." 0

# Field of interest (library name and weblink)
fdn=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "SampleName" | cut -d ":" -f 1)
fdr=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "Run" | cut -d ":" -f 1)
flk=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "download_path" | cut -d ":" -f 1)
fbp=$(head -n 1 <<< "$runinfo" | tr "," "\n" | grep -w -n "BioProject" | cut -d ":" -f 1)


#------------------#
# BioProject check #
#------------------#

# Check if non expected BioProject exist
list_bp=$(sed "1d" <<< "$runinfo" | cut -d "," -f "$fbp" | sort -u | grep -w -v "${bioproject}")

if [[ -n "$list_bp" ]]
then
    runinfo=$(grep -E -w -v "$(sed "s/ /|/g" <<< ${list_bp})" <<< "$runinfo")
    warning "Unexpected BioProject(s) detected and removed: $list_bp"
fi


#----------------#
# Fastq Download #
#----------------#

# Download fastq files
myend=$(tail -n +2 <<< "$runinfo" | sed "/^$/d" | wc -l)
while read line
do
    # Progress bar
    ((mystep++))
    [[ $mystep -lt $myend ]] && ProgressBar $mystep $myend "Download"

    # Filename, run and download link
    fln=$(awk -v fdn=$fdn -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '{print $fdn}' <<<$line)
    run=$(cut -d "," -f $fdr <<<$line)

    # Folder
    echo "$fln" >> "$log"
    [[ ! -d "$ldir/$fln/" ]] && mkdir -p "$ldir/$fln/"

    # Check for fastq if requested
    if [[ "$skip" ]]
    then
        ls_rn=$(find "$ldir/$fln/" -name "${run}_*.fastq*" | wc -l)
        [[ $ls_rn -gt 0 ]] && continue

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

    # Fetch SRA file in the background
    (myprefetch "$ldir" "$fln" "$run" "$log") &

    # Store path
    to_del+=("$ldir/$fln/$run")

done < <(tail -n +2 <<< "$runinfo" | sed "/^$/d")

# Wait until every background job is done
wait $(jobs -p)
ProgressBar $mystep $myend "Download"


#-----------#
# Log check #
#-----------#

# Check for errors
[[ $(grep -E -i "error|warning" "$log") ]] && error "Errors or warnings present in $log. Exiting..." 1


#----------------#
# Fastq renaming #
#----------------#

# Sanity check
flns=$(tail -n +2 <<< "$runinfo" | awk -v fdn=$fdn -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '{print $fdn}' | sort -u)
flns_nb=$(wc -l <<< "$flns")

[[ "$flns_nb" -ne "$myend" && -z "$merge" ]] && error "Several fastq files correspond to a single samples but merge was not used. Exiting..." 1

# Rename fastq files
unset to_del
mystep=0
myend=$flns_nb
while read fln
do
    # Progress bar
    ((mystep++))
    [[ $mystep -lt $myend ]] && ProgressBar $mystep $myend "Renaming"

    # Sample being processed
    echo "$fln" >> "$log"

    # Pause when number of running job is expected
    myjobs=($(jobs -p))
    while [[ ${#myjobs[@]} -ge $pl ]]
    do
        sleep 30s
        myjobs=($(jobs -p))
    done

    # Fetch SRA file in the background
    (rename_fastq "$ldir" "$fln" "$log") &

    # Store path
    to_del+=("$ldir/$fln/$fln")

done <<< "$flns"

# Wait until every background job is done
wait $(jobs -p)
ProgressBar $mystep $myend "Renaming"

exit 0
