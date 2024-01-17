#!/bin/bash

# Usage: krona_formatting.sh <your_input_file.csv> <your_sample_prefix_ID>

# Specify the input CSV file
input_file=$1

# Specify the string to search for in column names
search_string=$2

# Use awk to extract entire columns starting with the specified string
awk -v search="$search_string" -F ',' '{
  for (i = 1; i <= NF; i++) {
    if (index($i, search) == 1) {
      col[i] = 1
    }
  }
} FNR == 1 { # Process the header
  for (i = 1; i <= NF; i++) {
    if (col[i]) {
      printf "%s,", $i
    }
  }
  printf "\n"
} FNR > 1 { # Process the rest of the lines
  for (i = 1; i <= NF; i++) {
    if (col[i]) {
      printf "%s,", $i
    }
  }
  printf "\n"
}' "$input_file" | tail -n +2 | \
awk -F ',' '{
  sum = 0
  for (i = 1; i <= NF; i++) {
    sum += $i
  }
  print sum
}' > sum.tmp

# Extract taxonomy
columns="domain,phylum,class,order,family,genus,species"

# Use awk to cut columns by name
awk -F',' -v c="Seq" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i+1; break}; next} {for (i=p; i<=NF; i++) printf "%s%s", $i, (i==NF) ? "\n" : ","}' "$input_file" | sed 's/,/\t/g' > tax.tmp

# Combine sums and tax
paste sum.tmp tax.tmp > krona_formatted_file.tsv

# Clean up
rm tax.tmp sum.tmp
