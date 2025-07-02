#!/bin/bash


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# generate rename map
rename_map_file="${SCRIPT_DIR}/prototroch_rename.sed"

if [[ ! -f "$rename_map_file" ]]
then
    for i in {1..12}
    do
        for ap in {A,P}
        do
            newnum=$(( $i + 1))
            echo "s/prototroch_${i}clock_$ap/prototroch_${newnum}$ap/g" >> "$rename_map_file"
        done
    done

    sed -i 's/11A/11-12A/g' "$rename_map_file"
    sed -i 's/11clock_P/11clock/g' "$rename_map_file"
    sed -i 's/13/1/g' "$rename_map_file"
    sed -i '/11clock_A/d' "$rename_map_file"

else
    echo "$rename_map_file exists"
fi


# rename prototroch cells in text files
while read -r filename
do
    sed -f "$rename_map_file" -i "$filename"
done < <(grep -rl  --exclude='prototroch_rename.s*' "clock_" "${SCRIPT_DIR}/..")
