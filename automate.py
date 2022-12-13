import os
import argparse
import textwrap
import warnings

from njoy_utilities import utils


if __name__ == "__main__":

    # ------------------------------------------------------------
    # Setup Command-Line Interface
    # ------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
        A script for processing cross-sections for several isotopes.'''),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        Notes
        -----
        Materials are specified using the atomic number, element symbol,
        mass number, and optionally molecule name separated by dashes.
        Examples:
            Hydrogen in light water:   1-H-1-H2O
            Uranium 235            :   92-U-235
            Hydrogen in ZrH        :   1-H-1-ZrH
            Carbon in graphite     :   6-C-12-graphite
        ''')
    )

    parser.add_argument(
        '-m', '--materials',
        type=str, nargs='*', required=True, metavar='str',
        help=textwrap.dedent('''\
        A list of the desired materials.
        See the notes for specification instructions.
        ''')
    )

    parser.add_argument(
        '-gs', '--group_structures',
        type=str, nargs='*', required=True, metavar='str',
        help=textwrap.dedent('''\
        A list of the desired neutron, gamma, or neutron-gamma 
        group structures. If the latter, the group structures 
        should be separated by a `-`.
        ''')
    )

    parser.add_argument(
        '-t', '--temperatures',
        type=float, nargs='*', default=[293.6], metavar='float',
        help="The desired temperatures."
    )

    parser.add_argument(
        '--no_thermal',
        type=int, nargs='*', default=[], metavar='int',
        help="The indices of the materials to exclude thermal treatment for."
    )

    parser.add_argument(
        '--option',
        type=int, choices=[0, 1, 2], default=0, metavar='',
        help=textwrap.dedent('''\
        The mode to run the script in.
        0 runs NJOY and creates a \'.xs\' file, 1 only runs NJOY, and
        2 only creates a \'.xs\' file from an existing NJOY output.''')
    )

    parser.add_argument(
        "--plot", action='store_true', default=False,
        help="A flag for plotting."
    )

    argv = parser.parse_args()

    # ------------------------------------------------------------
    # Loop Over Materials
    # ------------------------------------------------------------

    for m, material in enumerate(argv.materials):
        material: str = material
        if material.count('-') < 2:
            raise ValueError("Invalid material specification.")

        # ----------------------------------------
        # Parse Material
        # ----------------------------------------

        material_info = utils.get_material_info(material)
        symbol = material_info['symbol']
        isotope = material_info['isotope'].replace('_', '-')
        molecule = material_info['molecule']
        neutron_endf = material_info['neutron_endf']
        gamma_endf = material_info['gamma_endf']
        photoat_endf = material_info['photoat_endf']

        # ----------------------------------------
        # Parse Thermal Options
        # ----------------------------------------

        # check whether thermal scattering is included
        with_thermal = m not in argv.no_thermal

        # get thermal scattering parameters
        sab_info = {}
        if molecule and with_thermal:
            sab_info = utils.get_thermal_info(symbol, molecule)

        # ----------------------------------------
        # Loop Over Group Structures
        # ----------------------------------------

        for gs in argv.group_structures:
            gs_list = gs.split('-')
            gs_info = utils.get_group_structure_info(gs_list)

            gs_outdir = gs_info['outdir']

            # ----------------------------------------
            # Loop over temperatures
            # ----------------------------------------

            for temperature in argv.temperatures:

                if temperature in [296.0, 293.6]:
                    temperature_name = "room"
                    if material_info['molecule'] == "H2O":
                        temperature = 293.6
                    elif material_info['molecule'] == "graphite":
                        temperature = 296.0

                elif temperature >= 0.0:
                    temperature_name = f"{str(temperature).replace('.', '_')}k"
                    if temperature.is_integer():
                        temperature_name = f"{temperature_name.split('_')[0]}k"

                else:
                    raise ValueError("Temperature must be positive.")

                outdir = os.path.join(gs_outdir, temperature_name)

                # ----------------------------------------
                # Prepare Output Location
                # ----------------------------------------

                os.makedirs(outdir, exist_ok=True)
                outdir = os.path.abspath(outdir)

                # ----------------------------------------
                # Define Filename
                # ----------------------------------------

                filename = material_info['symbol']
                filename += str(material_info['mass_number'])
                if with_thermal:
                    suffix = material_info['molecule']
                    suffix = suffix if suffix else "freegas"
                    filename = f"{filename}_{suffix}"

                # ----------------------------------------
                # Print Summary
                # ----------------------------------------

                msg = f"Processing isotope {isotope}, "
                if molecule:
                    msg += f"molecule {molecule}, "
                msg += f"temperature {temperature_name}..."
                print(msg)

                # ----------------------------------------
                # Write the Run and Process Script
                # ----------------------------------------
                with open("tmp.sh", "w") as f:

                    f.write("CWD=\"$PWD\"\n")
                    f.write("cd njoy_utilities || exit\n\n")

                    # -------------------- run njoy block
                    f.write(
                        "if [[ $1 == '0' ]] || [[ $1 == '1' ]]\n"
                        "then\n"
                    )

                    njoy = f"  python njoy_runner.py \\\n" \
                           f"  --njoy_executable=njoy21 \\\n" \
                           f"  --temperature={temperature} \\\n"

                    # neutron data
                    if neutron_endf and gs_info['neutron']:
                        njoy += f"  --path_to_neutron_endf={neutron_endf} \\\n"

                        # add group structure info
                        neutron_gs = gs_info['neutron']['gs_id']
                        njoy += f"  --neutron_group_structure={neutron_gs} \\\n"

                    # add custom group structure file
                    if neutron_gs == 1:
                        neutron_gs_file = gs_info['neutron']['gs_file']
                        njoy += f"  --custom_neutron_gs_file={neutron_gs_file} \\\n"

                    # gamma/photo-atomic data
                    print(gs_info['gamma'])
                    if gs_info['gamma']:
                        if not gamma_endf and not photoat_endf:
                            continue

                        if gamma_endf:
                            njoy += f"  --path_to_gamma_endf={gamma_endf} \\\n"
                        if photoat_endf:
                            njoy += f"  --path_to_photoat_endf={photoat_endf} \\\n"

                        # gamma group structure
                        gamma_gs = gs_info['gamma']['gs_id']
                        njoy += f"  --gamma_group_structure={gamma_gs} \\\n"

                        # add custom group structure file
                        if gamma_gs == 1:
                            gamma_gs_file = gs_info['gamma']['gs_file']
                            njoy += f"  --custom_gamma_gs_file={gamma_gs_file} \\\n"

                    # thermal scattering data
                    if sab_info:
                        sab_endf = sab_info['sab_endf']
                        njoy += f"  --path_to_sab_endf={sab_endf} \\\n"

                        # add inelastic thermal scattering
                        mti = sab_info['mti']
                        njoy += f"  --inelastic_thermal_number={mti} \\\n"

                        # add elastic thermal scattering
                        if "mtc" in sab_info:
                            mtc = sab_info['mtc']
                            njoy += f"  --elastic_thermal_number={mtc} \\\n"

                        # add the number of principal atoms
                        n_atoms = sab_info['n_atoms']
                        njoy += f"  --inelastic_thermal_num_atoms={n_atoms} \\\n"

                    # fission flag
                    if material_info['atomic_number'] >= 90:
                        njoy += "  --fissile \\\n"

                    # output data
                    njoy += f"  --output_directory={outdir} \\\n"
                    njoy += f"  --output_filename={filename}.njoy\n\n"

                    f.write(njoy)
                    f.write("fi\n\n")

                    # -------------------- processor block
                    f.write(
                        "if [[ $1 == '0' ]] || [[ $1 == '2' ]]\n"
                        "then\n"
                    )

                    process = f"  python3 njoy_processor.py \\\n" \
                              f"  --output_directory={outdir} \\\n" \
                              f"  --njoy_output_filename={filename}.njoy \\\n" \
                              f"  --xs_filename={filename}.xs"
                    process += f"  \\\n  --plot\n\n" if argv.plot else "\n\n"

                    f.write(process)
                    f.write("fi\n")

                    f.write("cd \"$CWD\" || exit\n")

                os.system(f"source tmp.sh {argv.option}")
                os.system("rm -f tmp.sh")
