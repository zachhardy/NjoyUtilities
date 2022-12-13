import os
import argparse
import textwrap
import warnings

from njoy_utilities import utils


if __name__ == "__main__":

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
        '-ngs', '--neutron_group_structure',
        type=str, required=True, metavar='str',
        choices=[
            'custom1', 'custom2', 'custom3', 'custom5',
            'custom7', 'custom31', 'lanl30', 'lanl70',
            'lanl80', 'lanl187', 'lanl618'
        ],
        help="The desired neutron group structure."
    )

    parser.add_argument(
        '-ggs', '--gamma_group_structure',
        type=str, metavar='str',
        choices=['lanl12', 'lanl24', 'lanl48'],
        help="The desired gamma group structure."
    )

    parser.add_argument(
        '-t', '--temperature',
        type=float, nargs='*', default=293.6, metavar='float',
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

    args = parser.parse_args()
    with_neutron = args.neutron_group_structure is not None
    with_gamma = args.gamma_group_structure is not None
    with_photoat = with_gamma

    ############################################################
    # Loop over materials
    ############################################################

    for m, material in enumerate(args.materials):
        material: str = material
        if material.count('-') < 2:
            raise ValueError("Invalid material specification.")

        ############################################################
        # Parse material argument
        ############################################################

        # get material data
        material_info = utils.get_material_info(material)
        symbol = material_info['symbol']
        isotope = material_info['isotope'].replace('_', '-')
        molecule = material_info['molecule']
        neutron_endf = material_info['neutron_endf']
        gamma_endf = material_info['gamma_endf']
        photoat_endf = material_info['photoat_endf']

        # check that neutron data exists
        if with_neutron:
            if not os.path.isfile(neutron_endf):
                raise FileNotFoundError(
                    f"{neutron_endf} is not a valid "
                    f"neutron ENDF file."
                )

        # if no gamma file, set gamma flag to False
        if with_gamma and not os.path.isfile(gamma_endf):
            warnings.warn("No gamma ENDF file for this isotope exists.")
            with_gamma = False

        # check photo-atomic data
        if with_photoat:
            if not os.path.isfile(photoat_endf):
                raise FileNotFoundError(
                    f"{photoat_endf} is not a valid photo-atomic ENDF file."
                )

        ############################################################
        # Determine thermalization options
        ############################################################

        # check whether thermal scattering is included
        with_thermal = m not in args.no_thermal

        # get thermal scattering parameters
        sab_info = {}
        if molecule and with_thermal:
            sab_info = utils.get_thermal_info(symbol, molecule)

        ############################################################
        # Determine group structure information
        ############################################################

        gs_info = utils.get_group_structure_info(
            args.neutron_group_structure, args.gamma_group_structure
        )
        output_directory = gs_info['output_directory']

        ############################################################
        # Determine temperatures
        ############################################################

        if args.temperature in [296.0, 293.6]:
            temperature_name = "room"
            if material_info['molecule'] == "H2O":
                args.temperature = 293.6
            elif material_info['molecule'] == "graphite":
                args.temperature = 296.0
        elif args.temperature >= 0.0:
            temperature_name = f"{str(args.temperature).replace('.', '_')}k"
        else:
            raise ValueError("Temperature must be positive.")

        output_directory = f"{output_directory}/{temperature_name}"

        ############################################################
        # Prepare output location
        ############################################################

        # create directory
        os.makedirs(output_directory, exist_ok=True)
        output_directory = os.path.abspath(output_directory)

        # define filename
        filename = material_info['symbol']
        filename += str(material_info['mass_number'])
        if with_thermal:
            suffix = material_info['molecule']
            suffix = suffix if suffix else "freegas"
            filename = f"{filename}_{suffix}"

        ############################################################
        # Write and run NJOY file
        ############################################################

        msg = f"Processing isotope {isotope}, "
        if molecule:
            msg += f"molecule {molecule}, "
        msg += f"temperature {temperature_name}..."
        print(msg)

        ############################################################
        # Write script to run NJOY and create cross-section file
        ############################################################

        with open("tmp.sh", "w") as f:

            f.write("CWD=\"$PWD\"\n")
            f.write("cd njoy_utilities || exit\n\n")

            ##################################################
            # Run NJOY
            ##################################################

            f.write(
                "if [[ $1 == '0' ]] || [[ $1 == '1' ]]\n"
                "then\n"
            )

            njoy = f"  python njoy_runner.py \\\n" \
                   f"  --njoy_executable=njoy21 \\\n" \
                   f"  --temperature={args.temperature} \\\n"

            # neutron data
            if with_neutron:
                njoy += f"  --path_to_neutron_endf={neutron_endf} \\\n"

                # add group structure info
                neutron_gs = gs_info['neutron']['gs_id']
                njoy += f"  --neutron_group_structure={neutron_gs} \\\n"

                # add custom group structure file
                if neutron_gs == 1:
                    neutron_gs_file = gs_info['neutron']['gs_file']
                    njoy += f"  --custom_neutron_gs_file={neutron_gs_file} \\\n"

            # gamma/photo-atomic data
            if with_gamma or with_photoat:
                if with_gamma:
                    njoy += f"  --path_to_gamma_endf={gamma_endf} \\\n"
                if with_photoat:
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

            # fissile
            if material_info['atomic_number'] >= 90:
                njoy += "  --fissile \\\n"

            # outputs
            njoy += f"  --output_directory={output_directory} \\\n"
            njoy += f"  --output_filename={filename}.njoy\n\n"

            f.write(njoy)
            f.write("fi\n\n")

            ##################################################
            # Run NJOY converter
            ###################################################

            f.write(
                "if [[ $1 == '0' ]] || [[ $1 == '2' ]]\n"
                "then\n"
            )

            process = f"  python3 njoy_processor.py \\\n" \
                      f"  --output_directory={output_directory} \\\n" \
                      f"  --njoy_output_filename={filename}.njoy \\\n" \
                      f"  --xs_filename={filename}.xs"
            process += f"  \\\n  --plot\n\n" if args.plot else "\n\n"

            f.write(process)
            f.write("fi\n")

            f.write("cd \"$CWD\" || exit\n")

        os.system(f"source tmp.sh {args.option}")
        os.system("rm -f tmp.sh")
