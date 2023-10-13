#!/bin/bash

clear

while [ 1 ];
  do
  # Display the main menu
  CHOICE=$(whiptail --title "Main Menu" --menu "Make your choice" $LINES $COLUMNS $(( $LINES - 8 )) \
    "1)" "Run Calculations"   \
    "2)" "Exit" 3>&2 2>&1 1>&3
  )

  # Process user choice
  case $CHOICE in
    "1)")
      # Display a message while running calculations
      TERM=ansi whiptail --infobox "Currently running calculations..." 10 60
      sleep 2  # Simulating some calculations

      # Display a message while reading files
      TERM=ansi whiptail --infobox "Reading files..." 10 60
      sleep 2  # Simulating file reading

      # Display a message when finished
      TERM=ansi whiptail --infobox "Finished processing." 10 60
      sleep 2  # Pause before returning to the main menu
    ;;
    "2")
      # Exit the script
      whiptail --msgbox "Exiting the script." 10 40
      exit
    ;;
  esac
done
