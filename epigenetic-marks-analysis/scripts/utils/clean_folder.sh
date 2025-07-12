#!/bin/bash

# --- Configuration ---
# Bash version check (associative arrays require Bash 4.0+)
if (( BASH_VERSINFO[0] < 4 )); then
  echo "Error: This script requires Bash version 4.0 or later." >&2
  exit 1
fi

# --- Variables ---
TARGET_DIR=""
URL_LIST_FILE="" # Renamed variable for clarity
DRY_RUN=0
FORCE_YES=0

# --- Functions ---
usage() {
  echo "Usage: $0 [--dry-run] [--yes] <directory_to_check> <url_list_file.txt>"
  echo "  Deletes files in <directory_to_check> whose names do NOT match the"
  echo "  filename part of the URLs listed in <url_list_file.txt>."
  echo
  echo "  Options:"
  echo "    --dry-run : Print files that would be deleted without actually deleting them."
  echo "    --yes     : Delete files without prompting for confirmation (use with caution!)."
  echo
  echo "  Arguments:"
  echo "    <directory_to_check> : The folder containing files to examine."
  echo "    <url_list_file.txt>  : A text file listing URLs (one per line)."
  echo "                           The script extracts filenames from these URLs."
  exit 1
}

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --dry-run)
      DRY_RUN=1
      shift # past argument
      ;;
    --yes)
      FORCE_YES=1
      shift # past argument
      ;;
    *)
      # Assume positional arguments if not options
      if [[ -z "$TARGET_DIR" ]]; then
        TARGET_DIR="$1"
      elif [[ -z "$URL_LIST_FILE" ]]; then
        URL_LIST_FILE="$1"
      else
        echo "Error: Unexpected argument '$1'" >&2
        usage
      fi
      shift # past argument
      ;;
  esac
done

# --- Input Validation ---
if [[ -z "$TARGET_DIR" ]] || [[ -z "$URL_LIST_FILE" ]]; then
  echo "Error: Missing required arguments." >&2
  usage
fi

if [[ ! -d "$TARGET_DIR" ]]; then
  echo "Error: Target directory '$TARGET_DIR' not found or is not a directory." >&2
  exit 1
fi

if [[ ! -f "$URL_LIST_FILE" ]]; then
  echo "Error: URL list file '$URL_LIST_FILE' not found or is not a file." >&2
  exit 1
fi

if [[ ! -r "$URL_LIST_FILE" ]]; then
  echo "Error: Cannot read URL list file '$URL_LIST_FILE'." >&2
  exit 1
fi


# --- Core Logic ---

# Read URLs and extract allowed filenames into an associative array
declare -A allowed_files
echo "Reading URLs from '$URL_LIST_FILE' and extracting filenames..."
while IFS= read -r url || [[ -n "$url" ]]; do
  # Skip empty lines or lines starting with # (comments)
  [[ -z "$url" ]] || [[ "$url" =~ ^\s*# ]] && continue

  # Extract filename from the URL using basename
  allowed_filename=$(basename "$url")

  # Add the extracted filename to the associative array if not empty
  if [[ -n "$allowed_filename" ]]; then
      allowed_files["$allowed_filename"]=1
      # You can uncomment the next line for debugging to see what's being added
      # echo "Debug: Allowing filename '$allowed_filename' from URL '$url'"
  else
      echo "Warning: Could not extract filename from line: $url" >&2
  fi
done < "$URL_LIST_FILE"
echo "Found ${#allowed_files[@]} unique allowed filenames from URLs."

# Find files to delete
declare -a files_to_delete
echo "Scanning directory '$TARGET_DIR'..."
shopt -s nullglob # Prevent loop from running if no files match
for filepath in "$TARGET_DIR"/*; do
  # Ensure it's a file, not a directory or something else
  if [[ -f "$filepath" ]]; then
    filename=$(basename "$filepath")
    # Check if the filename exists as a key in the allowed list
    if [[ -z "${allowed_files[$filename]}" ]]; then
      # The filename from the directory was NOT found in the list extracted from URLs
      files_to_delete+=("$filepath")
    fi
  fi
done
shopt -u nullglob # Turn off nullglob

# --- Execution/Reporting ---

if [[ ${#files_to_delete[@]} -eq 0 ]]; then
  echo "No files found in '$TARGET_DIR' that do not correspond to filenames in '$URL_LIST_FILE'. No action taken."
  exit 0
fi

echo "--------------------------------------------------"
if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY RUN: The following files are in '$TARGET_DIR' but their names were NOT found in the list derived from '$URL_LIST_FILE' and would be deleted:"
  printf "  %s\n" "${files_to_delete[@]}"
  echo "--------------------------------------------------"
  echo "DRY RUN complete. No files were deleted."
  exit 0
fi

echo "The following files are in '$TARGET_DIR' but their names were NOT found in the list derived from '$URL_LIST_FILE':"
printf "  %s\n" "${files_to_delete[@]}"
echo "--------------------------------------------------"

# Confirmation Prompt (unless --yes is used)
confirm="n"
if [[ $FORCE_YES -eq 1 ]]; then
  confirm="y"
  echo "Proceeding with deletion (--yes flag specified)."
else
  # Updated prompt message to reflect NZ locale and time (less critical for script function, but per instructions)
  current_nz_time=$(date) # Get current system time, assumes server locale is NZST or similar
  echo "Current time: $current_nz_time"
  read -p "Are you sure you want to delete these ${#files_to_delete[@]} files? (y/N): " confirm
fi

# Perform Deletion
if [[ "$confirm" =~ ^[Yy]$ ]]; then
  echo "Deleting files..."
  deleted_count=0
  error_count=0
  for file_to_del in "${files_to_delete[@]}"; do
    if rm "$file_to_del"; then
      echo "  Deleted '$file_to_del'"
      ((deleted_count++))
    else
      echo "  Error deleting '$file_to_del'" >&2
      ((error_count++))
    fi
  done
  echo "--------------------------------------------------"
  echo "Deletion complete. $deleted_count file(s) deleted."
  if [[ $error_count -gt 0 ]]; then
      echo "$error_count error(s) occurred during deletion." >&2
      exit 1 # Indicate partial failure
  fi
else
  echo "Deletion cancelled by user."
fi

exit 0