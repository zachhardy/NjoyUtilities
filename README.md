# NjoyUtilities
Prepares NJOY inputs and converts NJOY output to multi-group transport 
cross-section formats. This utility is meant to be expanded upon as needed to support 
a number of different multi-group transport cross-section formats.

## How to use this converter

### Step-0: Python requirements

Runs with python3.

### Step 1: Download and set up your flavor of ENDF

- The entire ENDF8.0 library (Preferred) 
  [ENDFVIII.0](https://www.nndc.bnl.gov/endf/b8.0/zips/ENDF-B-VIII.0.zip) 
  (~500Mb on download, ~2Gb extracted)
- The entire ENDF7.1 library 
  [ENDFVII.1](https://ndclx4.bnl.gov/gf/download/frsrelease/138/2242/ENDF-B-VII.1.tar.gz)
- Older versions. 
  Go to Brookhaven National Laboratory's site 
  [National Nuclear Data Center](https://www.nndc.bnl.gov/exfor/endf00.jsp) 
  and download what you need. 

Extract this library in a folder of your choice which we shall just call `ENDF_DIR`.

For example, on Linux machines, you can do this using: 
- ```shell 
  curl -O https://www.nndc.bnl.gov/endf/b8.0/zips/ENDF-B-VIII.0.zip
  ```
- unzip using the ```unzip``` command
- set the envrionement variable `ENDF_ROOT` to point to `ENDF_DIR`. 
  In bash, this is
  ```shell
  export ENDF_ROOT=<path-to-ENDF_DIR>
  ```
  **Note:** this is the absolute path.
  
### Step 2: Download and install NJOY2016 or NJOY21

This converter processes NJOY2016 or NJOY21 output files. 
Note that rigorous testing for each version has not been performed. 
NJOY2016 is hosted on GitHub at 
[https://github.com/njoy/NJOY2016](https://github.com/njoy/NJOY2016), 
and NJOY21 at [https://github.com/njoy/NJOY21](https://github.com/njoy/NJOY21). 
Installation instructions are located on each site and are summarized below for 
installing NJOY2016

Make sure you meet the prerequisites:

- C++17 or higher (If you have GCC or clang then you should be good)
- Fortran 2003 or higher
- Python 3.4+
- CMake 3.2+

After this is good to go, do the following. Go to a folder where you want to install NJOY2016

```shell
# Download the source code
git clone https://github.com/njoy/NJOY2016

# Get the desired version of NJOY21 (1.1.0 in this example)
cd NJOY2016
wget https://raw.githubusercontent.com/njoy/signatures/master/NJOY21/1.1.0-NJOY21.json
./metaconfigure/fetch_subprojects.py 1.1.0-NJOY21.json
```
Yes, we know the file points to NJOY21. We checked and it works for NJOY2016.

```shell
# Configure the build process
mkdir bin
cd bin
cmake -D CMAKE_BUILD_TYPE=Release ..

# Build NJOY21
make

# Test NJOY21
make test
```

Next add `njoy` to your environment variables. 
This may be different depending on your system. 
For `bash` users, in order to add the current `bin`-directory to your `PATH`-environment variable, 
type `pwd` and copy the path (let us use the place-holder `<path-to-njoy>` for this path). 
Next add the following to your `~/.bashrc` or `~/.bash_profile` file:

```shell
export PATH=<path-to-njoy>:$PATH
```

### Step 3: Clone the converter to a folder of your choice

```shell
cd <folder-of-choice>
git clone https://github.com/<username>/NjoyConverter
cd NjoyConverter
```

### Step 4a: Use the automation tool to generate the desired cross-sections

The automation script ```automate.py``` automatically writes input scripts of the same 
format as those found in ```examples```. 
The aim of the script is to provide a streamlined input for generating cross-sections.
Specification instructions can be found by calling ```python automate.py --help``` from
the command line.

### Step 4b: Adapt one of the example scripts to your specific application

The example scripts are a collection of neutron only examples:
- `example1a.sh`, U235 with a xmas172 neutron group structure
- `example1b.sh`, same as 1a but using graphite, showcasing how to include S(a,b) thermal scattering
- `example1c.sh`, same as 1b but now with a custom weighting spectrum
- `example1d.sh`, same as 1b but with completely custom group structure

And a collection of neutron-gamma examples:
- `example2a`, O16 with the lanl30 neutron group structure and the lanl12 gamma group structure
- `example2b`, O16 with the xmas172 neutron group structure and the lanl48 gamma group structure


The basic input for each example is:
```shell
if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

CWD=$PWD
cd njoy_utilities ||exit

#================================= Set properties here

export ENDF_ROOT=/Users/janv4/Desktop/Projects/ENDF/ENDF-B-VII.1
neutron_file="n-092_U_235.endf"

output_directory="output/ENDF-B-VII-1/xmas172/"
output_file_prefix="U235"

#================================= Run NJOY

python njoy_runner.py \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
--temperature=293.6 \
--neutron_group_structure=22 \
--output_directory=../$output_directory \
--output_filename=$output_file_prefix.njoy

cd "$CWD" || exit

echo "********** DONE with GENERATION"

python3 njoy_processor.py \
--output_directory=$output_directory \
--njoy_output_filename=$output_file_prefix.njoy \
--xs_filename=$output_file_prefix.xs \
--plot

echo "********** DONE with PROCESSING"

cd "$CWD" || exit
```

## FAQs:

### FAQ-1: General format of input examples

The two python scripts are `njoy_runner.py` and `njoy_processor.py` located 
in the ```njoy_utilities``` directory.

The `njoy_runner.py` script basically runs NJOY and has the following 
inputs (most of which are optional):
```
# --path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_gamma_endf=$ENDF_ROOT/gammas/$gamma_file \
# --path_to_photoat_endf=$ENDF_ROOT/photoat/$photoat_file \
# --path_to_sab_endf=$ENDF_ROOT/thermal_scatt/$sab_file \
# --neutron_group_structure=22 \
# --custom_neutron_gs_file="" \
# --neutron_weight_function=8 \
# --custom_neutron_wt_file="" \
# --gamma_group_structure=0 \
# --custom_gamma_gs_file="" \
# --gamma_weight_function=2 \
# --custom_gamma_wt_file="" \
# --temperature=296.0 \
# --inelastic_thermal_number=229 \
# --inelastic_thermal_num_atoms=1 \
# --elastic_thermal_number=230 \
# --output_directory=$output_directory \
# --output_filename=../$output_file_prefix.njoy \

```

The `njoy_processor.py` script converts NJOY output to the desired multi-group
tranport cross-section format. 
It only has three required inputs:
```
#--output_directory=$output_directory \
#--njoy_output_filename=$output_file_prefix.njoy \
#--xs_filename=$output_file_prefix.xs
```
Optionally a ```--plot``` flag can be used to visualize the data.

### FAQ-2: How to specify S(a,b) thermal scattering

Simply run `njoy_runner.py` with `--path_to_sab_endf`, 
`--inelastic_thermal_number=`, `--inelastic_thermal_num_atoms`, 
and `--elastic_thermal_number` specified.
The latter is only required for some molecules.

There aren't a lot of materials that have S(a,b) elastic treatment, 
so it is worth doing some homework on them and verifying a semi-infinite 
medium spectrum like done in the `tests` folder.

### FAQ-3: How to produce neutron-gamma cross-sections?
The moment you supply the option `--path_to_gamma_endf` or 
`--path_to_photoat_endf`, then the script will know to run with gamma 
production and make `--gamma_group_structure` required.

### FAQ-4: Format of a custom weighting spectrum file
The file is in ENDF TAB1-record format which can be confusing. 
An example spectrum file is supplied in 
`Docs/SampleWeightingSpectrum/spectrum_file.txt` and is the same for 
custom neutron AND gamma spectra.

**Note:** remember the `/` terminator and the blank line at the end of the file.

### FAQ-5: Format of a custom group structure file
The first line of a group structure file is the total number of groups `G`. 
Then followed by the lowest energy cutoff then a total of `G` upper bin boundaries 
(all in eV):

```
# Number of groups G
6
# G+1 boundaries, lower bound first then all upper bounds (eV)
1.0e-5
5.0e-2
5.0e-1
1.0e2
1.0e5
1.0e6
20.0e6/
```

**Note:** remember the `/` terminator and the blank line at the end of the file.
