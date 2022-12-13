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
--njoy_executable=njoy21 \
--path_to_neutron_endf="${ENDF_ROOT}/neutrons/$neutron_file" \
--path_to_sab="${ENDF_ROOT}/thermal_scatt/$sab_file" \
--inelastic_thermal_number=222 \
--inelastic_thermal_num_atoms=2 \
--temperature=293.6 \
--neutron_group_structure=1 \
--custom_neutron_gs_file="$output_directory/custom7g.txt" \
--output_directory=$output_directory \
--output_filename=$output_file_prefix.njoy

# --path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
# --inelastic_thermal_number=229 \
# --inelastic_thermal_num_atoms=1 \
# --elastic_thermal_number=230 \
# --path_to_gamma_endf= \
# --temperature=296.0 \
# --neutron_group_structure=22 \
# --neutron_weight_function=8 \
# --output_directory=$output_directory \
# --output_filename=$output_file_prefix.njoy \
# --gamma_group_structure=0 \
# --gamma_weight_function=2 \
# --custom_neutron_gs_file="" \
# --custom_gamma_gs_file="" \
# --custom_neutron_wt_file="" \
# --custom_gamma_wt_file="" \

echo "********** DONE with GENERATION"

python3 njoy_processor.py \
--output_directory=$output_directory \
--njoy_output_filename=$output_file_prefix.njoy \
--xs_filename=$output_file_prefix.xs \
--plot

echo "********** DONE with PROCESSING"

cd "$CWD" || exit
