#!/bin/bash

# Directories
POS_DIR="positive"
NEG_DIR="negative"
TRAIN_POS="train/positive"
TRAIN_NEG="train/negative"
VAL_POS="val/positive"
VAL_NEG="val/negative"

mkdir -p "$TRAIN_POS" "$TRAIN_NEG" "$VAL_POS" "$VAL_NEG"

# Get shuffled file lists
POS_FILES=($(printf '%s\n' "$POS_DIR"/* | sort -R))
NEG_FILES=($(printf '%s\n' "$NEG_DIR"/* | sort -R))

TOTAL_POS=${#POS_FILES[@]}
TOTAL_NEG=${#NEG_FILES[@]}

# 80% of negatives sets the cap for balanced training
TRAIN_COUNT=$(( TOTAL_NEG * 80 / 100 ))

echo "Total positive: $TOTAL_POS"
echo "Total negative: $TOTAL_NEG"
echo "Balanced training size per class: $TRAIN_COUNT"
echo "Positive val/excess: $(( TOTAL_POS - TRAIN_COUNT ))"
echo "Negative val: $(( TOTAL_NEG - TRAIN_COUNT ))"

# Copy positives: first TRAIN_COUNT to train, rest to val
for i in "${!POS_FILES[@]}"; do
  if (( i < TRAIN_COUNT )); then
    cp "${POS_FILES[$i]}" "$TRAIN_POS/"
  else
    cp "${POS_FILES[$i]}" "$VAL_POS/"
  fi
done

# Copy negatives: first TRAIN_COUNT to train, rest to val
for i in "${!NEG_FILES[@]}"; do
  if (( i < TRAIN_COUNT )); then
    cp "${NEG_FILES[$i]}" "$TRAIN_NEG/"
  else
    cp "${NEG_FILES[$i]}" "$VAL_NEG/"
  fi
done

echo "Done!"

