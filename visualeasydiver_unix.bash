#! /bin/bash

# This program is just an example of how to make a whiptail menu and some basic commands.
# Copyright (C) 2016  Baljit Singh Sarai

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

clear

function useEasyDiver {
  {
    options=$(whiptail --title "Select Options" --checklist "Choose EasyDIVER options by hitting the space key. Press enter to confirm:" 20 60 10 \
        "Input Directory" "" ON \
        "Output Directory" "" OFF \
        "Forward Primer Sequence" "" OFF \
        "Reverse Primer Sequence" "" OFF \
        "Retain Individual Lane Outputs" "" OFF \
        "Number of Threads" "" OFF \
        "Translate to Amino Acids" "" OFF \
        "Extra Flags for PANDASeq" "" OFF \
        3>&1 1>&2 2>&3)
    echo "options: $options"
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
#eval `resize`
CHOICE=$(
whiptail --title "Main Menu" --menu "Make your choice" $LINES $COLUMNS $(( $LINES - 8 )) \
	"1)" "Use EasyDIVER"   \
	"2)" "Calculate enrichment statistics"  \
	"3)" "Help"  \
	"4)" "End script"  3>&2 2>&1 1>&3
)

result=""
case $CHOICE in
	"1)")
	  useEasyDiver
	;;
	"2)")   
	        OP=$(uptime | awk '{print $3;}')
		result="This system has been up $OP minutes"
		whiptail --msgbox "$result" 20 78
	;;

	"3)")
		result="easydiver.sh is a Bash script that takes input files in fastq format and performs various processing
steps on them. It accepts command-line arguments for specifying input and output directories, forward and
reverse primer sequences, number of threads, extra flags for PANDASeq, translation into amino acids, and retaining
individual line outputs.
v2.0 Author: Allison Tee and Celia Blanco
contact: ateecup@stanford.edu or cblanco@chem.ucsb.edu

USAGE: bash easydiver.sh -i [-n -o -p -q -r -T -h -a -e -s]
where:
	REQUIRED
	-i input directory filepath

	OPTIONAL
	-o output directory filepath
	-p forward primer sequence for extraction
	-q reverse primer sequence for extraction
	-a translate to amino acids
	-r retain individual lane outputs
	-T # of threads
	-e extra flags for PANDASeq (use quotes, e.g. \"-L 50\")
	-s include enrichment stats on output (Only if you have in/neg/out files)"
	whiptail --msgbox "$result" 20 78
  ;;
	"4)")   
#	    		whiptail --msgbox "Thank you for using EasyDIVER!" 20 78
		        exit
        ;;
esac
done
exit