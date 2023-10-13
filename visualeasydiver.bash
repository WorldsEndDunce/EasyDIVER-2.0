#! /bin/bash

# Define colors for the terminal interface
export NEWT_COLORS='
root=,blue
checkbox=,brightblue
entry=,brightblue
title=brightgreen,gray
textbox=,lightgray
window=,gray
sellistbox=,green
actsellistbox=,green
button=,green
actbutton=,green
actcheckbox=,green
border=,gray
'

# Clear the terminal
clear

# Function to use EasyDIVER
function useEasyDiver {
  {
    # Use whiptail to create a checklist for selecting EasyDIVER options
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

    if [ ! -z "$options" ]; then
        selected_options="Selected options:\n$options"
        whiptail --title "Options Selected" --msgbox "$selected_options" 20 60

        # Construct the EasyDIVER command based on selected options
        command="bash easydiver.sh"

        if [[ $options == *"Input Directory"* ]]; then
            input_dir=$(whiptail --inputbox "Enter input directory filepath:" 10 60 3>&1 1>&2 2>&3)
            command+=" -i \"$input_dir\""
        fi
        # ... Other options ...

        # Run the constructed command and display a progress gauge
        eval "$command" | whiptail --gauge "Running EasyDIVER..." 6 60 0
    else
        whiptail --title "Options Selected" --msgbox "No options selected" 10 40
    fi
  }
}

# Function to find enrichments
function findEnrichments {
  {
    type=""
    if (whiptail --title "Find Enrichments" --yesno "Calculate enrichment statistics for amino acid counts? (Yes- AA, No- Nucleotide)" 8 78); then
      type="counts.aa"
    else
      type="counts"
    fi

    # Get the EasyDIVER output directory path
    dir=$(whiptail --inputbox "Enter the filepath for the EasyDIVER output directory:" 10 60 3>&1 1>&2 2>&3)

    # ... Rest of the code to calculate enrichments ...

  }| whiptail --gauge "Finding enrichments..." 6 60 0
}

# Function to graph figures
function graphFigures {
  {
    choice=$(whiptail --title "Graphs" --radiolist \
    "Choose an option:" $LINES $COLUMNS $(( $LINES - 8 )) \
    "Histogram" "Visualize sequence length distribution for a given output/histos folder" ON \
    "Enrichment Scatterplot" "Plot negative control vs output enrichment of each sequence for a given round." OFF \
    "AA Count Line Graph" "Plot unique and total amino acid counts across all rounds" OFF \
    "Txt to Excel" "Convert a given enrichment result output (.txt), convert it to .xslx" OFF \
    3>&1 1>&2 2>&3)

    filepath=""
    graph=""

    if [ $? -eq 0 ]; then
        if [ "$choice" == "Histogram" ]; then
            filepath=$(whiptail --inputbox "Enter EasyDIVER output directory filepath:" 10 60 3>&1 1>&2 2>&3)
            graph="2"
            python3 ./graphs.py "$filepath" "$graph" | whiptail --gauge "Graphing..." 6 60 0
            whiptail --title "Finished" --msgbox "Look in the figures folder to find your output(s)!" 10 40
        elif [ "$choice" == "Enrichment Scatterplot" ]; then
            filepath=$(whiptail --inputbox "Enter Enrichment result .txt filepath:" 10 60 3>&1 1>&2 2>&3)
            graph="1"
            python3 ./graphs.py "$filepath" "$graph" | whiptail --gauge "Graphing..." 6 60 0
            whiptail --title "Finished" --msgbox "Look in the figures folder to find your output(s)!" 10 40
        elif [ "$choice" == "AA Count Line Graph" ]; then
            filepath=$(whiptail --inputbox "Enter EasyDIVER output directory filepath:" 10 60 3>&1 1>&2 2>&3)
            graph="3"
            python3 ./graphs.py "$filepath" "$graph" | whiptail --gauge "Graphing..." 6 60 0
            whiptail --title "Finished" --msgbox "Look in the figures folder to find your output(s)!" 10 40
        elif [ "$choice" == "Txt to Excel" ]; then
            filepath=$(whiptail --inputbox "Enter Enrichment result .txt filepath:" 10 60 3>&1 1>&2 2>&3)
            python3 ./txt_to_xslx.py "$filepath" | whiptail --gauge "Working on it..." 6 60 0
            whiptail --title "Finished" --msgbox "Look in the figures folder to find your output(s)!" 10 40
        fi
    else
        echo "Operation canceled."
    fi
  }
}

# Array containing ASCII art for the welcome message
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
# Display a welcome message with ASCII art
whiptail --title "Welcome" --msgbox "${ascii_art[0]}
Welcome to the pipeline for Easy pre-processing and Dereplication of In Vitro Evolution Reads!\n\n
Press OK to continue." 30 85

# Main menu loop
while [ 1 ];
do
CHOICE=$(
whiptail --title "Main Menu" --menu "Choose an option:"  $LINES $COLUMNS $(( $LINES - 8 )) \
	"1)" "Use EasyDIVER"   \
	"2)" "Calculate enrichment statistics"  \
	"3)" "Figures and Misc." \
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
	  graphFigures
  ;;
	"4)")
		result="EasyDIVER 2.0 is a Bash program that takes input files of selex round reads in fastq(.gz) format and analyzes them, providing
sequence type and length distribution, joined reads, and/or enrichment information. There are also options to visualize sequence distribution,
enrichment, and diversity data.
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