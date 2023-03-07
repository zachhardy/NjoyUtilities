#
# Simple cross-section for hydrogen using a custom neutron
# group structure
# FIX ME

if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

CWD=$PWD
cd ../njoy_utilities || exit

#================================= Set properties here
neutron_file="n-001_H_001.endf"
sab_file="tsl-HinH2O.endf"

output_directory="../output/ENDF-B-VIII-0/custom7/"
output_file_prefix="H_water"

#================================= Run NJOY

python3 generate_njoy_mgxs.py \
--njoy-executable=njoy21 \
--path-to-neutron-endf="${ENDF_ROOT}/neutrons/$neutron_file" \
--path_to_sab="${ENDF_ROOT}/thermal_scatt/$sab_file" \
--inelastic-thermal-number=222 \
--inelastic-thermal-num-atoms=2 \
--temperature=293.6 \
--neutron-group-structure=1 \
--custom-neutron-gs-file="$output_directory/custom7g.txt" \
--output-directory=$output_directory \
--output-filename=$output_file_prefix.njoy

# --path-to-neutron-endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
# --inelastic-thermal-number=229 \
# --inelastic-thermal-num-atoms=1 \
# --elastic-thermal-number=230 \
# --path-to-gamma-endf= \
# --temperature=296.0 \
# --neutron-group-structure=22 \
# --neutron-weight-function=8 \
# --output-directory=$output_directory \
# --output-filename=$output_file_prefix.njoy \
# --gamma-group-structure=0 \
# --gamma-weight-function=2 \
# --custom-neutron-gs-file="" \
# --custom-gamma-gs-file="" \
# --custom-neutron-wt-file="" \
# --custom-gamma-wt-file="" \

echo "********** DONE with GENERATION"

python3 njoy_processor.py \
--output-directory=$output_directory \
--njoy-output-filename=$output_file_prefix.njoy \
--xs-filename=$output_file_prefix.xs \
--plot

echo "********** DONE with PROCESSING"

cd "$CWD" || exit
