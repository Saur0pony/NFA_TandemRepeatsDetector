#!/bin/bash


# Parse command-line options

function show_usage (){
    echo "Usage $0 [options]"
    echo ""
    echo "Options:"
    echo "-h, display this help"
    echo "-f, input fasta file"
    echo "-p, the pattern to search in sequence"
    echo "-k, maximum error allowed"
    echo "-g, gap length between 2 motifs"
    echo "-c, cutpat"
    echo "-o, output file"
return 0
}

if [[ $# -eq 0 ]];then # If no argument given, display the help
   show_usage
   exit 1
fi

# To check mandatory arguments are provided
f_flag=false # file
p_flag=false # pattern
o_flag=false # output file


options=':f:p:k:g:co:h'
while getopts $options option
do
    case "$option" in
       h  )
         show_usage; exit 0;;
       f  )
         f_flag=true; file_arg=$OPTARG;;
       p  )
         p_flag=true; pattern_arg=$OPTARG;;
       o  )
         o_flag=true; output_arg=$OPTARG;;
       k  )
         k_arg=$OPTARG;;
       g  )
         g_arg=$OPTARG;;
       c  )
         c_arg=true;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done


if ! $f_flag
then
    echo "You must provide a fasta file as input (-f)" >&2
    exit 1
fi
if ! $p_flag
then
    echo "You must provide a pattern (-p)" >&2
    exit 1
fi
if ! $o_flag
then
    echo "You must provide an output filename (-o)" >&2
    exit 1
fi


# Compute python command to run script
command="pypy ./NFA_TandermRepeat_detector_v5.py --file GENOME_FILE --pattern $pattern_arg --out OUTPUT"

if [ $k_arg ]
then
  add_to_command=" --k $k_arg"
  tmp_command=$command$add_to_command
  command=$tmp_command
fi
if [ $g_arg ]
then
  add_to_command=" --gap $g_arg"
  tmp_command=$command$add_to_command
  command=$tmp_command
fi
if [ $c_arg ]
then
  add_to_command=" --cutpat"
  tmp_command=$command$add_to_command
  command=$tmp_command
fi


# Split fasta file into chromosomes
csplit -s -n 4 -f tmp_genome_part_ -z $file_arg '/>/' '{*}'

for gen_part in ./tmp_genome_part_*
do
  basename_gen_part=$(basename $gen_part)
  output=${output_arg}_${basename_gen_part}
  new_command=$command
  new_command=${new_command/GENOME_FILE/$gen_part}
  new_command=${new_command/OUTPUT/$output}
  tmp_script=./tmp_script_${basename_gen_part}.sh
  echo '#!/bin/bash' > $tmp_script
  echo " " >> $tmp_script
  echo $new_command >> $tmp_script
  echo " " >> $tmp_script
  echo "rm $basename_gen_part" >> $tmp_script
  #echo rm $tmp_script >> $tmp_script

  nb_running_jobs=$(squeue | grep "gprunier" | wc -l)
  waiting_time=20
  echo $nb_running_jobs
  while [ $nb_running_jobs -ge 25 ]
  do
    echo "More than 20 running scripts. Waiting $waiting_time seconds"
    sleep $waiting_time
    if [ $waiting_time -le 600 ]; then
      waiting_time=$((waiting_time*2))
    fi
    nb_running_jobs=$(squeue | grep "gprunier" | wc -l)
  done

  sbatch --cpus-per-task=6 --mem=10G $tmp_script
done


