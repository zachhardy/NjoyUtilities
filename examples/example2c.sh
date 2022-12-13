#
# Simple cross-section for O16 using mostly defaults and 
# the lanl-187 neutron group structure and the gamma lanl-48 group structure
#

if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

CWD=$PWD
cd ../njoy_utilities || exit

#================================= Set properties here
neutron_file="n-008_O_016.endf"
gamma_file="photoat-008_O_000.endf"
# sab_file="tsl-graphite.endf"

output_directory="../output/ENDF-B-VIII-0/LANL187_LANL48/"
output_file_prefix="O16_n187g48"

#================================= Run NJOY

python3 generate_njoy_mgxs.py \
--njoy_executable=njoy21 \
--path_to_neutron_endf="${ENDF_ROOT}/neutrons/$neutron_file" \
--path_to_gamma_endf="${ENDF_ROOT}/photoat/$gamma_file" \
--temperature=293.6 \
--neutron_group_structure=10 \
--neutron_weight_function=11 \
--gamma_group_structure=6 \
--gamma_weight_function=3 \
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
