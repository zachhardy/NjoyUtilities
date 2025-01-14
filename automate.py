"""
An automation script that can run NJOY, process outputs, and write cross-section
data files for combinations of listed materials, group structures, and temperatures.
"""

import os
import argparse
import textwrap
import warnings

from njoy_utilities import utils


################################################################################
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.MetavarTypeHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


if __name__ == "__main__":

    # ----------------------------------------------------------------------
    # Setup Command-Line Interface
    # ----------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
        A script for processing cross-sections for several isotopes.'''),
        formatter_class=CustomFormatter,
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
        type=str,
        nargs='*',
        required=True,
        help=textwrap.dedent('''\
        A list of the desired materials.
        See the notes for specification instructions.
        ''')
    )

    parser.add_argument(
        '-gs', '--group-structures',
        type=str,
        nargs='*',
        required=True,
        help=textwrap.dedent('''\
        A list of the desired neutron, gamma, or neutron-gamma 
        group structures. If the latter, the group structures 
        should be separated by a `-`.
        ''')
    )

    parser.add_argument(
        '-t', '--temperatures',
        type=float,
        nargs='*',
        default=[293.6],
        help="The desired temperatures."
    )

    parser.add_argument(
        '--no-thermal',
        type=int,
        nargs='*',
        default=[],
        help="The indices of the materials to exclude thermal treatment for."
    )

    parser.add_argument(
        '--option',
        type=int,
        choices=[0, 1, 2],
        default=0,
        help=textwrap.dedent('''\
        The mode to run the script in.
        0 runs NJOY and creates a \'.xs\' file, 1 only runs NJOY, and
        2 only creates a \'.xs\' file from an existing NJOY output.''')
    )

    parser.add_argument(
        '--plot',
        action='store_true',
        default=False,
        help="A flag for plotting."
    )

    argv = parser.parse_args()

    # ----------------------------------------------------------------------
    # Loop over specified materials
    # ----------------------------------------------------------------------

    for m, material in enumerate(argv.materials):
        material: str = material
        if material.count('-') < 2:
            raise ValueError("Invalid material specification.")

        # ------------------------------------------------------------
        # Parse material info
        # ------------------------------------------------------------

        material_info = utils.get_material_info(material)
        symbol = material_info['symbol']
        isotope = material_info['isotope'].replace('_', '-')
        molecule = material_info['molecule']
        neutron_endf = material_info['neutron_endf']
        gamma_endf = material_info['gamma_endf']
        photoat_endf = material_info['photoat_endf']

        # ------------------------------------------------------------
        # Parse thermal options
        # ------------------------------------------------------------

        # check whether thermal scattering is included
        with_thermal = m not in argv.no_thermal

        # get thermal scattering parameters
        sab_info = {}
        if molecule and with_thermal:
            sab_info = utils.get_thermal_info(symbol, molecule)

        # ------------------------------------------------------------
        # Loop over specified group structures
        # ------------------------------------------------------------

        for gs in argv.group_structures:
            gs_list = gs.split('-')
            gs_info = utils.get_group_structure_info(gs_list)

            gs_outdir = gs_info['outdir']

            # --------------------------------------------------
            # Loop over specified temperatures
            # --------------------------------------------------

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

                # --------------------------------------------------
                # Prepare output location
                # --------------------------------------------------

                os.makedirs(outdir, exist_ok=True)
                outdir = os.path.abspath(outdir)

                filename = material_info['symbol']
                filename += str(material_info['mass_number'])
                if material_info['molecule']:
                    filename = f"{filename}_{material_info['molecule']}"
                elif not with_thermal:
                    filename = f"{filename}_notherm"

                # --------------------------------------------------
                # Print summary
                # --------------------------------------------------

                print()
                print(f"##################################################")
                print(f"{'Isotope':<20}: {isotope}")
                if molecule:
                    print(f"{'Molecule':<20}: {molecule}")
                print(f"{'Group Structure':<20}: {gs}")
                print(f"{'Temperature':<20}: {temperature}")
                print(f"##################################################")
                print()

                # --------------------------------------------------
                # Write the run and process script
                # --------------------------------------------------
                with open("tmp.sh", "w") as f:

                    f.write("CWD=\"$PWD\"\n")
                    f.write("cd njoy_utilities || exit\n\n")

                    # ----------------------------------------
                    # Run NJOY block
                    # ----------------------------------------

                    f.write(
                        "if [[ $1 == '0' ]] || [[ $1 == '1' ]]\n"
                        "then\n"
                    )

                    njoy = f"  python njoy_runner.py \\\n" \
                           f"  --njoy-executable=njoy21 \\\n" \
                           f"  --temperature={temperature} \\\n"

                    # neutron data
                    if neutron_endf and gs_info['neutron']:
                        njoy += f"  --path-to-neutron-endf={neutron_endf} \\\n"

                        # add group structure info
                        neutron_gs = gs_info['neutron']['gs_id']
                        neutron_n_grps = gs_info['neutron']['n_groups']
                        njoy += f"  --neutron-group-structure={neutron_gs} \\\n"
                        njoy += f"  --neutron-num-groups={neutron_n_grps} \\\n"

                    # add custom group structure file
                    if neutron_gs == 1:
                        neutron_gs_file = gs_info['neutron']['gs_file']
                        njoy += f"  --custom-neutron-gs-file={neutron_gs_file} \\\n"

                    # gamma/photo-atomic data
                    if gs_info['gamma']:
                        if not gamma_endf and not photoat_endf:
                            continue

                        if gamma_endf:
                            njoy += f"  --path-to-gamma-endf={gamma_endf} \\\n"
                        if photoat_endf:
                            njoy += f"  --path-to-photoat-endf={photoat_endf} \\\n"

                        # gamma group structure
                        gamma_gs = gs_info['gamma']['gs_id']
                        gamma_n_grps = gs_info['gamma']['n_groups']
                        njoy += f"  --gamma-group-structure={gamma_gs} \\\n"
                        njoy += f"  --gamma-num-groups={gamma_n_grps} \\\n"

                        # add custom group structure file
                        if gamma_gs == 1:
                            gamma_gs_file = gs_info['gamma']['gs_file']
                            njoy += f"  --custom-gamma-gs-file={gamma_gs_file} \\\n"

                    # thermal scattering data
                    if sab_info:
                        sab_endf = sab_info['sab_endf']
                        njoy += f"  --path-to-sab-endf={sab_endf} \\\n"

                        # add inelastic thermal scattering
                        mti = sab_info['mti']
                        njoy += f"  --inelastic-thermal-number={mti} \\\n"

                        # add elastic thermal scattering
                        if "mtc" in sab_info:
                            mtc = sab_info['mtc']
                            njoy += f"  --elastic-thermal-number={mtc} \\\n"

                        # add the number of principal atoms
                        n_atoms = sab_info['n_atoms']
                        njoy += f"  --inelastic-thermal-num-atoms={n_atoms} \\\n"

                    # fission flag
                    if material_info['atomic_number'] >= 90:
                        njoy += "  --fissile \\\n"

                    # output data
                    njoy += f"  --output-directory={outdir} \\\n"
                    njoy += f"  --output-filename={filename}.njoy\n\n"

                    f.write(njoy)
                    f.write("fi\n\n")

                    # ----------------------------------------
                    # NJOY Processor block
                    # ----------------------------------------

                    f.write(
                        "if [[ $1 == '0' ]] || [[ $1 == '2' ]]\n"
                        "then\n"
                    )

                    process = f"  python3 njoy_processor.py \\\n" \
                              f"  --output-directory={outdir} \\\n" \
                              f"  --njoy-output-filename={filename}.njoy \\\n" \
                              f"  --xs-filename={filename}.xs"
                    process += f"  \\\n  --plot\n\n" if argv.plot else "\n\n"

                    f.write(process)
                    f.write("fi\n")

                    f.write("cd \"$CWD\" || exit\n")

                os.system(f"source tmp.sh {argv.option}")
                os.system("rm -f tmp.sh")
