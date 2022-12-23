"""Prepares input/output files for NJOY2016 and executes NJOY2016/2021"""
import argparse
import os
import sys
import textwrap
import warnings


def parse_isotope_info(f):
    iso = f.readline(11)
    element = iso[4:6].strip()
    Z, A = int(iso[:3]), int(iso[7:10])
    m = iso[10:].strip()
    return Z, element, A, m


def parse_material_number(f):
    lne = f.readline()
    while "MATERIAL" not in lne:
        lne = f.readline()
    return int(lne.split()[2])


def check_isotope_info(f):
    iso = f.readline(11)
    if int(iso[:3]) != atomic_num or \
            iso[4:6].strip() != symbol or \
            int(iso[7:10]) != mass_num or \
            iso[10:].strip() != metastable:
        raise AssertionError("Differing isotope encountered.")


def check_element_info(f):
    iso = f.readline(11)
    if int(iso[:3]) != atomic_num or \
            iso[4:6].strip() != symbol:
        raise AssertionError("Different element encountered.")


############################################################
# Check python version
############################################################

if sys.version_info[0] < 3:
    print(f"\nError: This script requires python3 but was executed "
          f"with version:\n\n{sys.version}\n")
    sys.exit(1)

############################################################
# Setup parser
############################################################

ign = [*range(1, 35)]  # Allowed numeric options for --neutron_group_structure
niwt = [*range(1, 13)]  # Allowed numeric options for --neutron_weight_function
igg = [*range(0, 11)]  # Allowed numeric options for --gamma_group_structure
giwt = [*range(1, 4)]  # Allowed numeric options for --gamma_weight_function
argparser = argparse.ArgumentParser(
    description="Prepares input/output files "
                "for NJOY2016 and executes NJOY2016",
    formatter_class=argparse.RawTextHelpFormatter,
    epilog='''Additional information:
  Neutron group structure options (ign):
   1   arbitrary structure (read in)
   2   csewg 239-group structure
   3   lanl 30-group structure
   4   anl 27-group structure
   5   rrd 50-group structure
   6   gam-i 68-group structure
   7   gam-ii 100-group structure
   8   laser-thermos 35-group structure
   9   epri-cpm 69-group structure
   10  lanl 187-group structure
   11  lanl 70-group structure
   12  sand-ii 620-group structure
   13  lanl 80-group structure
   14  eurlib 100-group structure
   15  sand-iia 640-group structure
   16  vitamin-e 174-group structure
   17  vitamin-j 175-group structure
   18  xmas nea-lanl
   19  ecco 33-group structure
   20  ecco 1968-group structure
   21  tripoli 315-group structure
   22  xmas lwpc 172-group structure
   23  vit-j lwpc 175-group structure
   24  shem cea 281-group structure
   25  shem epm 295-group structure
   26  shem cea/epm 361-group structure
   27  shem epm 315-group structure
   28  rahab aecl 89-group structure
   29  ccfe 660-group structure (30 MeV)
   30  ukaea 1025-group structure (30 MeV)
   31  ukaea 1067-group structure (200 MeV)
   32  ukaea 1102-group structure (1 GeV)
   33  ukaea 142-group structure (200 MeV)
   34  lanl 618-group structure
   
   Neutron group structure weighting options (iwt):
   1   read in smooth weight function
   2   constant
   3   1/e
   4   1/e + fission spectrum + thermal maxwellian
   5   epri-cell lwr
   6   (thermal) -- (1/e) -- (fission + fusion)
   7   same with t-dep thermal part
   8   thermal--1/e--fast reactor--fission + fusion
   9   claw weight function
   10  claw with t-dependent thermal part
   11  vitamin-e weight function (ornl-5505)
   12  vit-e with t-dep thermal part
   
  Gamma group structure options (igg):
   0   none
   1   arbitrary structure (read in)
   2   csewg 94-group structure
   3   lanl 12-group structure
   4   steiner 21-group gamma-ray structure
   5   straker 22-group  structure
   6   lanl 48-group structure
   7   lanl 24-group structure
   8   vitamin-c 36-group structure
   9   vitamin-e 38-group structure
   10  vitamin-j 42-group structure
   
   Gamma group structure weighting options (giwt):
   1   read in smooth weight function
   2   constant
   3   1/e + rolloffs
   
Example Custom group structure file:
------
# Number of groups
31/
# Group boundaries
1.0000e-05 4.6589e+05 9.3178e+05 1.3977e+06 1.8636e+06 2.3295e+06
2.7953e+06 3.2612e+06 3.7271e+06 4.1930e+06 4.6589e+06 5.1248e+06
5.5907e+06 6.0566e+06 6.5225e+06 6.9884e+06 7.4775e+06 7.9434e+06
8.4093e+06 8.8752e+06 9.3411e+06 9.8070e+06 1.0273e+07 1.0739e+07
1.1205e+07 1.1671e+07 1.2136e+07 1.2602e+07 1.3068e+07 1.3534e+07
1.3999e+07 1.4000e+07/
------

Example Custom smooth weight function:
See spectrum_file.txt
   
''')

argparser.add_argument(
    '--njoy_executable', type=str, default='njoy21', metavar='',
    help="The NJOY executable"
)

argparser.add_argument(
    '--path_to_neutron_endf', type=str, metavar='',
    help="Path to an incident neutron ENDF file."
)

argparser.add_argument(
    '--path_to_gamma_endf', type=str, metavar='',
    help="Path the incident gamma ENDF file."
)

argparser.add_argument(
    '--path_to_photoat_endf', type=str, metavar='',
    help="Path to photo-atomic ENDF file."
)

argparser.add_argument(
    '--path_to_sab_endf', type=str, metavar='',
    help="Path to S(alpha, beta) scattering ENDF file."
)

argparser.add_argument(
    '--neutron_group_structure', type=int,
    choices=ign, default=3, metavar='',
    help=textwrap.dedent('''\
    The neutron group structure.
    If 1, --custom_neutron_gs_file is required.
    Default 3, LANL 30-group structure.''')
)

argparser.add_argument(
    '--custom_neutron_gs_file', type=str, metavar='',
    help=textwrap.dedent('''\
    The path to the custom neutron group structure file.
    If --neutron_group_structure is not 1, this is ignored.'''),
)

argparser.add_argument(
    '--neutron_weight_function', type=int,
    choices=niwt, default=6, metavar='',
    help=textwrap.dedent('''\
    The neutron weight function.
    If 1, --custom_neutron_wt_file is required.
    Default 8, thermal--1/e--fast reactor--fission + fusion.''')
)

argparser.add_argument(
    '--custom_neutron_wt_file', type=str, metavar='',
    help=textwrap.dedent('''\
    The path to the custom neutron weight function file.
    If --neutron_weight_function is not 1, this is ignored.'''),
)

argparser.add_argument(
    '--gamma_group_structure', type=int,
    choices=igg, default=0, metavar='',
    help=textwrap.dedent('''\
    The gamma group structure.
    If 1, --custom_gamma_gs_file is required.
    Default 0, no gamma group structure.''')
)

argparser.add_argument(
    '--custom_gamma_gs_file', type=str, metavar='',
    help=textwrap.dedent('''\
    The path to the custom gamma group structure file.
    If --gamma_group_structure is not 1, this is ignored.'''),
)

argparser.add_argument(
    '--gamma_weight_function', type=int,
    choices=giwt, default=3, metavar='',
    help=textwrap.dedent('''\
    The gamma weight function.
    If 1, --custom_gamma_wt_file is required.
    Default 3, 1/e + rolloffs.''')
)

argparser.add_argument(
    '--custom_gamma_wt_file', type=str, metavar='',
    help=textwrap.dedent('''\
    The path to the custom gamma weight function file.
    If --gamma_weight_function is not 1, this is ignored.'''),
)

argparser.add_argument(
    '--temperature', type=float, default=296.0, metavar='',
    help="The material temperature."
)

argparser.add_argument(
    '--inelastic_thermal_number', type=int, metavar='',
    help="MT number to use for incoherent inelastic scattering"
)

argparser.add_argument(
    '--elastic_thermal_number', type=int, metavar='',
    help="MT number to use for coherent/incoherent elastic scattering"
)

argparser.add_argument(
    '--inelastic_thermal_num_atoms', type=int, default=1, metavar='',
    help="MT number to use for incoherent inelastic scattering",
)

argparser.add_argument(
    '--no_thermal', action='store_true', default=False,
    help="A flag for excluding any thermal scattering."
)

argparser.add_argument(
    '--fissile', action='store_true', default=False,
    help="A flag for fissile materials."
)

argparser.add_argument(
    '--output_directory', type=str, default=os.getcwd(), metavar='',
    help=textwrap.dedent('''\
    A directory to save the output to.
    If unspecified, output are saved to the current directory.''')
)

argparser.add_argument(
    '--output_filename', type=str, metavar='',
    help=textwrap.dedent('''\
    A filename to save the output to.
    If unspecified, output are given a generic name.''')
)

argv = argparser.parse_args()

############################################################
# Check arguments
############################################################

if argv.path_to_neutron_endf:
    if not os.path.isfile(argv.path_to_neutron_endf):
        raise FileNotFoundError(
            f"Value supplied to --path_to_neutron_endf does not "
            f"point to an existing file."
        )

if argv.path_to_gamma_endf:
    if not os.path.isfile(argv.path_to_gamma_endf):
        raise FileNotFoundError(
            "Value supplied to --path_to_gamma_endf does not point "
            "to an existing file."
        )

if argv.path_to_photoat_endf:
    if not os.path.isfile(argv.path_to_photoat_endf):
        raise FileNotFoundError(
            "Value supplied to --path_to_photoat_endf does not point "
            "to an existing file."
        )

if argv.path_to_sab_endf:
    if not os.path.isfile(argv.path_to_sab_endf):
        print(argv.path_to_sab_endf)
        raise FileNotFoundError(
            "Value supplied to --path_to_sab_endf does not point to "
            "an existing file."
        )

if argv.neutron_group_structure == 1:
    if not argv.custom_neutron_gs_file:
        raise ValueError(
            "When supplying 1 for --neutron_group_structure, "
            "then --custom_neutron_gs_file must be supplied."
        )
    elif not os.path.isfile(argv.custom_neutron_gs_file):
        raise FileNotFoundError(
            "Value supplied to --custom_neutron_gs_file does not "
            "point to an existing file."
        )

if argv.neutron_weight_function == 1:
    if not argv.custom_neutron_wt_file:
        raise ValueError(
            "When supplying 1 for --neutron_weight_function, "
            "then --custom_neutron_wt_file must be supplied."
        )
    elif not os.path.isfile(argv.custom_neutron_wt_file):
        raise FileNotFoundError(
            "Value supplied to --custom_neutron_wt_file does not "
            "point to an existing file."
        )

if argv.path_to_gamma_endf or argv.path_to_photoat_endf:
    if argv.gamma_group_structure == 1:
        if not argv.custom_gamma_gs_file:
            raise ValueError(
                "When supplying 1 for --gamma_group_structure, "
                "then --custom_gamma_gs_file must be supplied."
            )
        elif not os.path.isfile(argv.custom_gamma_gs_file):
            raise FileNotFoundError(
                "Value supplied to --custom_gamma_gs_file does not "
                "point to an existing file."
            )

    if argv.gamma_weight_function == 1:
        if not argv.custom_gamma_wt_file:
            raise ValueError(
                "When supplying 1 for --gamma_weight_function, "
                "then --custom_gamma_wt_file must be supplied."
            )
        elif not os.path.isfile(argv.custom_gamma_wt_file):
            raise FileNotFoundError(
                "Value supplied to --custom_gamma_wt_file does not "
                "point to an existing file."
            )

if argv.path_to_sab_endf:
    if not argv.inelastic_thermal_number:
        raise ValueError(
            "When --path_to_sab_endf is supplied, "
            "--inelastic_thermal_number must also be supplied."
        )

output_directory = os.path.abspath(argv.output_directory)
output_directory = os.path.join(output_directory, "njoy")
if not os.path.isdir(output_directory):
    warnings.warn(
        "Value supplied to --output_directory does not point to an "
        f"existing directory. Creating one at {output_directory}..."
    )
    os.makedirs(output_directory)
output_filename = argv.output_filename

with_neutron = argv.path_to_neutron_endf is not None
with_gamma = argv.path_to_gamma_endf is not None
with_photoat = argv.path_to_photoat_endf is not None
with_sab = argv.path_to_sab_endf is not None
with_thermal = not argv.no_thermal

############################################################
# Tape reservations
############################################################

# The following tape numbers are used:
# neutron endf         20
# neutron sab endf     50
# gamma endf           60
# photoat endf         70
# moder neutron  input 20    output 21
# moder gamma    input 60    output 61
# moder photoat  input 70    output 71
# reconr neutron input 21    output 22
# reconr gamma   input 61    output 62
# reconr photoat input 71    output 72
# broadr neutron input 21,22 output 23
# broadr gamma   input 61 62 output 63
# unresr         input 21,23 output 24
# heatr          input 21,24 output 25
# thermr freeg   input 25    output 26
# thermr sab     input 50,26 output 27
# groupr neutron input 21,27 output 28
# groupr gamma   input 61,63 output 64
# gaminr         input 71,72 output 73
# moder print    input 28    output 29
# moder print    input 64    output 65
# moder print    input 73    output 74

#############################################################
# Display type of problem
############################################################

if argv.path_to_neutron_endf:
    if argv.path_to_gamma_endf:
        print("Neutron + Gamma Problem")
    else:
        print("Neutron only Problem")
elif argv.path_to_gamma_endf:
    print("Gamma only Problem")
else:
    raise AssertionError(
        "A neutron or gamma ENDF file must be supplied."
    )

############################################################
# Create input tapes
############################################################

if argv.path_to_neutron_endf:
    os.system(f"ln -fs {argv.path_to_neutron_endf} tape20")

if argv.path_to_sab_endf:
    os.system(f"ln -fs {argv.path_to_sab_endf} tape50")

if argv.path_to_gamma_endf:
    os.system(f"ln -fs {argv.path_to_gamma_endf} tape60")

if argv.path_to_photoat_endf:
    os.system(f"ln -fs {argv.path_to_photoat_endf} tape70")

############################################################
# Extract general information from ENDF files
############################################################

atomic_num = None
symbol = None
mass_num = None
metastable = None
found_iso_info = False

neutron_material_number = None
gamma_material_number = None
photoat_material_number = None
sab_material_number = None

sab_material_name = None

# Neutron ENDF
if argv.path_to_neutron_endf:
    with open("tape20", "r") as endf:
        for _ in range(5):
            endf.readline()

        if not found_iso_info:
            isotope_info = parse_isotope_info(endf)
            atomic_num, symbol, mass_num, metastable = isotope_info
            found_iso_info = True
        else:
            check_isotope_info(endf)

        neutron_material_number = parse_material_number(endf)
        print(f"neutron material number read: "
              f"{neutron_material_number}")

# Gamma ENDF
if argv.path_to_gamma_endf:
    with open("tape60", "r") as endf:
        for _ in range(5):
            endf.readline()

        if not found_iso_info:
            isotope_info = parse_isotope_info(endf)
            atomic_num, symbol, mass_num, metastable = isotope_info
            found_iso_info = True
        else:
            check_isotope_info(endf)

        gamma_material_number = parse_material_number(endf)
        print(f"gamma material number read: "
              f"{gamma_material_number}")

# Photo-atomic ENDF
if argv.path_to_photoat_endf:
    with open("tape70", "r") as endf:
        for _ in range(5):
            endf.readline()

        check_element_info(endf)
        photoat_material_number = parse_material_number(endf)
        print(f"photo-atomic material number read:"
              f"{photoat_material_number}")

# S(alpha, beta) ENDF
if with_thermal and argv.path_to_sab_endf:
    with open("tape50", "r") as endf:
        for _ in range(5):
            endf.readline()

        sab_material_name = endf.readline().split()[0]
        sab_material_number = parse_material_number(endf)
        print(f"s(a, b) material {sab_material_name} read: "
              f"{sab_material_number}")

############################################################
# Define output filename
############################################################

if not output_filename:
    output_filename = f"{symbol}{mass_num}{metastable}"
    if with_thermal and not with_sab:
        output_filename = f"{output_filename}_freegas"
    output_filename = f"{output_filename}.txt"

############################################################
# Start writing NJOY input deck
############################################################

with open("NJOY_INPUT.txt", "w") as njoy_input:
    njoy_input.write("-- Processing ENDF to PENDF\n")

    # ============================== MODER inputs
    njoy_input.write("moder\n")
    njoy_input.write("20 -21/\n")
    if with_gamma:
        njoy_input.write("moder\n")
        njoy_input.write("60 -61/\n")
    if with_photoat:
        njoy_input.write("moder\n")
        njoy_input.write("70 -71\n")

    # ============================== RECONR inputs
    njoy_input.write("reconr\n")
    njoy_input.write("-21 -22/\n")
    njoy_input.write("'pendf neutron tape'/\n")
    njoy_input.write(f"{neutron_material_number} 0/\n")
    njoy_input.write("0.001/\n")
    njoy_input.write("0/\n")

    if with_gamma:
        njoy_input.write("reconr\n")
        njoy_input.write("-61 -62/\n")
        njoy_input.write("'pendf gamma tape'/\n")
        njoy_input.write(f"{gamma_material_number} 0/\n")
        njoy_input.write("0.001/\n")
        njoy_input.write("0/\n")

    if with_photoat:
        njoy_input.write("reconr\n")
        njoy_input.write("-71 -72/\n")
        njoy_input.write("'pendf photo-atomic tape'/\n")
        njoy_input.write(f"{photoat_material_number} 0/\n")
        njoy_input.write("0.001/\n")
        njoy_input.write("0/\n")

    # ============================== BROADR inputs
    njoy_input.write("broadr\n")
    njoy_input.write("-21 -22 -23/\n")
    njoy_input.write(f"{neutron_material_number} 1 0 0 0/\n")
    njoy_input.write("0.001/\n")
    njoy_input.write(f"{argv.temperature}/\n")
    njoy_input.write("0/\n")

    if with_gamma:
        njoy_input.write("broadr\n")
        njoy_input.write("-61 -62 -63/\n")
        njoy_input.write(f"{gamma_material_number} 1 0 0 0/\n")
        njoy_input.write("0.001/\n")
        njoy_input.write(f"{argv.temperature}/\n")
        njoy_input.write("0/\n")

    # ============================== UNRESR inputs
    njoy_input.write("unresr\n")
    njoy_input.write("-21 -23 -24/\n")
    njoy_input.write(f"{neutron_material_number} 1 1 0/\n")
    njoy_input.write(f"{argv.temperature}/\n")
    njoy_input.write("0.0/\n")
    njoy_input.write("0/\n")

    # ============================== HEATR inputs
    njoy_input.write("heatr\n")
    njoy_input.write("-21 -24 -25/\n")
    njoy_input.write(f"{neutron_material_number} 0/\n")

    # ============================== THERMR inputs
    if with_thermal:
        # free-gas
        njoy_input.write("thermr\n")
        njoy_input.write("0 -25 -26/\n")

        njoy_input.write("0 ")
        njoy_input.write(f"{neutron_material_number} ")
        njoy_input.write("16 ")  # number of angle bins
        njoy_input.write("1 ")  # number of temperatures
        # inelastic option -- 0=none, 1=free-gas, 2=s(a,b)
        njoy_input.write("1 ")
        # elastic option -- 0=none, 1=ENDF6 format
        njoy_input.write("0 ")
        njoy_input.write("0 ")  # output format
        njoy_input.write("1 ")  # number of principal atoms
        njoy_input.write("221 ")  # MT for inelastic reactions
        njoy_input.write("1/\n")  # print option -- 0=min, 1=max

        njoy_input.write(f"{argv.temperature}/\n")  # temperatures

        njoy_input.write("0.005 ")  # tolerance
        njoy_input.write("5.0/\n")  # maximum energy for thermal

        # S(alpha, beta)
        if with_sab:
            njoy_input.write("thermr\n")
            njoy_input.write("50 -26 -27/\n")

            njoy_input.write(f"{sab_material_number} ")
            njoy_input.write(f"{neutron_material_number} ")
            njoy_input.write("16 ")  # number of angle bins
            njoy_input.write("1 ")  # number of temperaturs
            # inelastic option -- 0=none, 1=free-gas, 2=s(a,b)
            njoy_input.write("2 ")
            # elastic option -- 0=none, 1=ENDF6 format
            njoy_input.write("0 ")
            njoy_input.write("0 ")  # output format
            njoy_input.write(f"{argv.inelastic_thermal_num_atoms} ")
            njoy_input.write(f"{argv.inelastic_thermal_number} ")
            njoy_input.write(f"1/\n")  # prinnt option -- 0=min, 1=max

            njoy_input.write(f"{argv.temperature}/\n")  # temperatures

            njoy_input.write("0.005 ")  # tolerance
            njoy_input.write("5.0/\n")  # maximum energy for thermal

    # ============================== GROUPR inputs
    njoy_input.write("groupr\n")
    njoy_input.write("-21 -27 0 -28/\n"
                     if with_thermal and with_sab else
                     "-21 -26 0 -28/\n"
                     if with_thermal and not with_sab else
                     "-21 -25 0 -28/\n")

    njoy_input.write(f"{neutron_material_number} ")
    njoy_input.write(f"{argv.neutron_group_structure} ")
    njoy_input.write(f"{argv.gamma_group_structure} ")
    njoy_input.write(f"{argv.neutron_weight_function} ")
    njoy_input.write("7 ")  # legendre order
    njoy_input.write("1 ")  # number of temperatures
    njoy_input.write("1 ")  # number of sigma zeros
    njoy_input.write("1 ")  # long print option -- 0=min, 1=max
    njoy_input.write("1/\n")  # smoothing operator -- 0=min, 1=max

    title = f"{symbol}{mass_num}{metastable}"
    njoy_input.write(f"'{title}'/\n")  # title

    njoy_input.write(f"{argv.temperature}/\n")  # temperatures
    njoy_input.write("0.0/\n")  # sigma zeros value

    # Custom neutron group structure
    if argv.neutron_group_structure == 1:
        with open(argv.custom_neutron_gs_file, "r") as wt_file:
            for line in wt_file:
                if line[0] != '#':
                    njoy_input.write(line)

    # Custom gamma group structure
    if argv.gamma_group_structure == 1:
        with open(argv.custom_gamma_gs_file, "r") as wt_file:
            for line in wt_file:
                if line[0] != "#":
                    njoy_input.write(line)

    # Custom neutron weight function
    if argv.neutron_weight_function == 1:
        with open(argv.custom_neutron_wt_file, "r") as wt_file:
            for line in wt_file:
                if line[0] != '#':
                    njoy_input.write(line)

    # MF3 cross-sections
    njoy_input.write("3/\n")

    # thermal cross-sections
    if with_thermal:
        njoy_input.write("3 221 'free gas'/\n")
        if with_sab:
            njoy_input.write(f"3 {argv.inelastic_thermal_number} "
                             f"'inelastic s(a,b)'/\n")
            if argv.elastic_thermal_number:
                njoy_input.write(f"3 {argv.elastic_thermal_number} "
                                 f"'elastic_s(a,b)'/\n")

    # fission data
    if argv.fissile:
        njoy_input.write("3 452 'total nubar (neutron)'/\n")

        njoy_input.write("3 456 'prompt nubar (neutron)'/\n")
        njoy_input.write("5 18 'prompt chi (neutron)'/\n")

        njoy_input.write("3 455 'delayed nubar (neutron)'/\n")
        njoy_input.write("5 455 'delayed chi (neutron)'/\n")

    # special cross-sections
    njoy_input.write("3 257 'average energy (neutron)'/\n")
    njoy_input.write("3 259 'inverse velocity'/\n")

    # MF6 transfer matrices
    njoy_input.write("6/\n")

    # thermal transfer matrices
    if with_thermal:
        njoy_input.write("6 221 'free gas neutron matrix'/\n")
        if with_sab:
            njoy_input.write(f"6 {argv.inelastic_thermal_number} ")
            njoy_input.write("'inelastic s(a,b) neutron matrix'/\n")

            if argv.elastic_thermal_number:
                njoy_input.write(f"6 {argv.elastic_thermal_number} "
                                 f"'elastic s(a,b) neutron matrix'/\n")

    # fission matrix
    if argv.fissile:
        njoy_input.write("6 18/ '(n,fission) neutron matrix'/\n")

    # MF13 photon production cross-sections
    njoy_input.write("13/\n" if with_gamma else "")

    # MF16 neutron-gamma matrices
    njoy_input.write("16/\n" if with_gamma else "")

    # End groupr
    njoy_input.write("0/\n")  # terminate reactions
    njoy_input.write("0/\n")  # terminate groupr

    if with_gamma:
        njoy_input.write("groupr\n")
        njoy_input.write("-61 -63 0 -64/\n")

        njoy_input.write(f"{gamma_material_number} ")
        njoy_input.write(f"{argv.neutron_group_structure} ")
        njoy_input.write(f"{argv.gamma_group_structure} ")
        njoy_input.write(f"{argv.neutron_weight_function} ")
        njoy_input.write("7 ")  # legendre order
        njoy_input.write("1 ")  # number of temperatures
        njoy_input.write("1 ")  # number of sigma zeros
        njoy_input.write("1 ")  # long print option -- 0=min, 1=max
        njoy_input.write("1/\n")  # smoothing operator -- 0=min, 1=max

        title = f"{symbol}{mass_num}{metastable}"
        njoy_input.write(f"'{title}'/\n")  # title

        njoy_input.write(f"{argv.temperature}/\n")  # temperatures
        njoy_input.write("0.0/\n")  # sigma zeros value

        # Custom neutron group structure
        if argv.neutron_group_structure == 1:
            with open(argv.custom_neutron_gs_file, "r") as wt_file:
                for line in wt_file:
                    if line[0] != '#':
                        njoy_input.write(line)

        # Custom gamma group structure
        if argv.gamma_group_structure == 1:
            with open(argv.custom_gamma_gs_file, "r") as wt_file:
                for line in wt_file:
                    if line[0] != "#":
                        njoy_input.write(line)

        # Custom neutron weight function
        if argv.neutron_weight_function == 1:
            with open(argv.custom_neutron_wt_file, "r") as wt_file:
                for line in wt_file:
                    if line[0] != '#':
                        njoy_input.write(line)

        # MF3 cross-sections
        njoy_input.write("3/\n")

        njoy_input.write("3 257 'average energy (gamma)'/\n")

        if argv.fissile:
            njoy_input.write("3 452 'total nubar (gamma)'/\n")

            njoy_input.write("3 456 'prompt nubar (gamma)'/\n")
            njoy_input.write("5 18 'prompt chi (gamma)'/\n")

            njoy_input.write("3 455 'delayed nubar (gamma)'/\n")
            njoy_input.write("5 455 'delayed chi (gamma)'/\n")

        njoy_input.write("6/\n")  # MF6 transfer matrices

        if argv.fissile:
            njoy_input.write("6 18/ '(g,fission) neutron matrix'/\n")

        # MF16 neutron-gamma matrices
        njoy_input.write("16/\n" if with_gamma else "")

        # End groupr
        njoy_input.write("0/\n")  # terminate reactions
        njoy_input.write("0/\n")  # terminate groupr

    # ============================== GAMINR inputs
    if with_photoat:
        njoy_input.write("gaminr\n")
        njoy_input.write("-71 -72 0 -73/\n")

        njoy_input.write(f"{photoat_material_number} ")
        njoy_input.write(f"{argv.gamma_group_structure} ")
        njoy_input.write(f"{argv.gamma_weight_function} ")
        njoy_input.write("7 ")  # legendre order
        njoy_input.write("1/\n")  # print options

        njoy_input.write(f"'Photon Interactions'/\n")  # title

        # Custom gamma group structure
        if argv.gamma_group_structure == 1:
            with open(argv.custom_gamma_gs_file) as wt_file:
                for line in wt_file:
                    if line[0] != "#":
                        njoy_input.write(line)

        # Custom gamma weight function
        if argv.gamma_weight_function == 1:
            with open(argv.custom_gamma_wt_file) as wt_file:
                for line in wt_file:
                    if line[0] != "#":
                        njoy_input.write(line)

        njoy_input.write("-1/\n")  # all reaction data
        njoy_input.write("0/\n")  # terminal gaminr
    njoy_input.write("stop\n")

os.system("rm -f out")
if argv.njoy_executable == "njoy21":
    cmd_line = argv.njoy_executable + " -i NJOY_INPUT.txt -o out"
else:
    cmd_line = argv.njoy_executable + " < NJOY_INPUT.txt"
print('command line =', cmd_line)

os.system(cmd_line)
os.system("rm tape* NJOY_INPUT.txt")

output_path = os.path.join(output_directory, output_filename)
print(f"Copying output file to {output_path}")
os.system(f"cp out {output_path} && rm -f out")
