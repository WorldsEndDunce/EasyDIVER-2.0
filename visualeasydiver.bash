#! /bin/bash

clear

function useEasyDiver {
  {
    options=$(whiptail --title "Select EasyDIVER Options" --checklist "Choose EasyDIVER options by hitting the space key. Press enter to confirm:" 20 60 10 \
        "Input Directory" "" ON \
        "Output Directory" "" OFF \
        "Forward Primer Sequence" "" OFF \
        "Reverse Primer Sequence" "" OFF \
        "Retain Individual Lane Outputs" "" OFF \
        "Number of Threads" "" OFF \
        "Translate to Amino Acids" "" OFF \
        "Extra Flags for PANDASeq" "" OFF \
        3>&1 1>&2 2>&3)
#    echo "options: $options"
    if [ $? -eq 0 ]; then
        selected_options="Selected options:\n$options"
        whiptail --title "Options Selected" --msgbox "$selected_options" 20 60

        command="bash easydiver.sh"
        if [[ $options == *"Input Directory"* ]]; then
            input_dir=$(whiptail --inputbox "Enter input directory filepath:" 10 60 3>&1 1>&2 2>&3)
            command+=" -i \"$input_dir\""
        fi

        if [[ $options == *"Output Directory"* ]]; then
            output_dir=$(whiptail --inputbox "Enter output directory filepath:" 10 60 3>&1 1>&2 2>&3)
            command+=" -n \"$output_dir\""
        fi

        if [[ $options == *"Forward Primer Sequence"* ]]; then
            forward_primer=$(whiptail --inputbox "Enter forward primer sequence:" 10 60 3>&1 1>&2 2>&3)
            command+=" -p \"$forward_primer\""
        fi

        if [[ $options == *"Reverse Primer Sequence"* ]]; then
            reverse_primer=$(whiptail --inputbox "Enter reverse primer sequence:" 10 60 3>&1 1>&2 2>&3)
            command+=" -q \"$reverse_primer\""
        fi

        if [[ $options == *"Retain Individual Lane Outputs"* ]]; then
            command+=" -r"
        fi

        if [[ $options == *"Number of Threads"* ]]; then
            num_threads=$(whiptail --inputbox "Enter number of threads:" 10 60 3>&1 1>&2 2>&3)
            command+=" -T \"$num_threads\""
        fi

        if [[ $options == *"Translate to Amino Acids"* ]]; then
            command+=" -a"
        fi

        if [[ $options == *"Extra Flags for PANDASeq"* ]]; then
            extra_flags=$(whiptail --inputbox "Enter extra flags for PANDASeq:" 10 60 3>&1 1>&2 2>&3)
            command+=" -e \"$extra_flags\""
        fi
        # Run the constructed command
        eval "$command" | whiptail --gauge "Running the command..." 6 60 0
    else
        whiptail --title "Options Selected" --msgbox "No options selected" 10 40
    fi
    }
}

function findEnrichments {
  {
    type=""
    if (whiptail --title "Find Enrichments" --yesno "Calculate enrichment statistics for amino acid counts? (Yes- AA, No- Nucleotide)" 8 78); then
      type="counts.aa"
    else
      type="counts"
    fi
#        command="bash scripts_enrichments/modified_counts_bash.sh"
        dir=$(whiptail --inputbox "Enter the filepath for the EasyDIVER output directory:" 10 60 3>&1 1>&2 2>&3)
#        command+=" \"$dir\""
      #!/bin/bash
    {
      echo type
      outdir=$dir
      counts_dir="$outdir/$type" # $outdir/counts or $outdir/counts.aa

      # Get the maximum round
      max_round=0
      for file in "$counts_dir"/*-out_$type.txt; do
          filename=$(basename "$file")
          round=${filename%-out_$type.txt}

          if [ $round -gt $max_round ]; then
              max_round=$round
          fi
      done
          progress=$((5))
          echo $progress
      # Check if there are any files matching the format "*-in_counts.aa.txt"
      if [ ! -n "$(find "$counts_dir" -name "*-in_${type}.txt" -print -quit)" ]; then
          # Cases 1A and 1B: Loop up to max_round - 1
          for ((i = 1; i < max_round; i++)); do
              # Run the modified_counts_bash.sh script with the appropriate arguments
              if [ ! -n "$(find "$counts_dir" -name "*-neg_${type}.txt" -print -quit)" ]; then
                python3 ./scripts_enrichments/modified_counts.py -out "$counts_dir/$i-out_$type.txt" -res "./scripts_enrichments/$i-res.txt"
              else
                python3 ./scripts_enrichments/modified_counts.py -out "$counts_dir/$i-out_$type.txt" -neg "$counts_dir/$(($i + 1))-neg_$type.txt" -res "./scripts_enrichments/$i-res.txt"
              fi

              # Calculate progress
              progress=$((i * 100 / (max_round - 1)))
              echo $progress
          done
      else
          # Case 2A and 2B: Loop up to max_round
          for ((i = 1; i <= max_round; i++)); do
              # Run the modified_counts_bash.sh script with the appropriate arguments
              if [ ! -n "$(find "$counts_dir" -name "*-neg_${type}.txt" -print -quit)" ]; then
                python3 ./scripts_enrichments/modified_counts.py -in "$counts_dir/$i-in_$type.txt" -out "$counts_dir/$i-out_$type.txt" -res "./scripts_enrichments/$i-res.txt"
              else
                python3 ./scripts_enrichments/modified_counts.py -in "$counts_dir/$i-in_$type.txt" -out "$counts_dir/$i-out_$type.txt" -neg "$counts_dir/$i-neg_$type.txt" -res "./scripts_enrichments/$i-res.txt"
              fi

              # Calculate progress
              progress=$((i * 100 / max_round))
              echo $progress
          done
      fi
    }| whiptail --gauge "Finding enrichments..." 6 60 0
    }
}

ascii_art=(
$'
  ______                _____ _______      ________ _____    ___    ___
 |  ____|              |  __ \_   _\ \    / /  ____|  __ \  |__ \  / _ \
 | |__   __ _ ___ _   _| |  | || |  \ \  / /| |__  | |__) |    ) || | | |
 |  __| / _` / __| | | | |  | || |   \ \/ / |  __| |  _  /    / / | | | |
 | |___| (_| \__ \ |_| | |__| || |_   \  /  | |____| | \ \   / /_ | |_| |
 |______\__,_|___/\__, |_____/_____|   \/   |______|_|  \_\ |____(_)___/
                   __/ |
                  |___/
                  '
)
whiptail --title "Welcome" --msgbox "${ascii_art[0]}
Welcome to the pipeline for Easy pre-processing and Dereplication of In Vitro Evolution Reads!\n\n
Press OK to continue." 30 85

while [ 1 ];
do
CHOICE=$(
whiptail --title "Main Menu" --menu "Choose an option:"  $LINES $COLUMNS $(( $LINES - 8 )) \
	"1)" "Use EasyDIVER"   \
	"2)" "Calculate enrichment statistics"  \
	"3)" "Graph figures (TODO)" \
	"4)" "Help"  \
	"5)" "End script"  3>&2 2>&1 1>&3
)

result=""
case $CHOICE in
	"1)")
	  useEasyDiver
	;;
	"2)")   
	  findEnrichments
	;;
	"3)")
	whiptail --msgbox "$result" 20 78
  ;;
	"4)")
		result="EasyDIVER 2.0 is a Bash program that takes input files of selex round reads in fastq(.gz) format and analyzes them, providing
sequence type and length distribution, joined reads, and/or enrichment information.
v2.0 Author: Allison Tee and Celia Blanco
contact: ateecup@stanford.edu or cblanco@chem.ucsb.edu

	REQUIRED
	Input directory filepath

	OPTIONAL MODIFIERS
	Output directory filepath
	Forward primer sequence for extraction
	Reverse primer sequence for extraction
	Translating to amino acids
	Retaining individual lane outputs
	# of threads
	extra flags for PANDASeq (use quotes, e.g. \"-L 50\")"
	whiptail --msgbox "$result" $LINES $COLUMNS $(( $LINES - 8 ))
  ;;
	"5)")
		    exit
  ;;
esac
done
exit