import os
import sys
import argparse
import textwrap


def get_material_info(material: str) -> dict:
    """
    Return the ENDF paths for a particular material input.

    Parameters
    ----------
    material : str
        The material input formatted via <Z>-<element>-<A>-<molecule>.

    Returns
    -------
    A dictionary containing the ENDF file paths.
    """
    entries = material.split('-')
    symbol = entries[1]
    Z, A = int(entries[0]), int(entries[2])

    isotope = "_".join([str(Z).zfill(3), symbol, str(A).zfill(3)])
    element = "_".join([str(Z).zfill(3), symbol, "000"])
    molecule = None if len(entries) < 4 else entries[3]

    endf_root = os.environ['ENDF_ROOT']
    return {'isotope': isotope, 'molecule': molecule,
            'atomic_number': Z, 'symbol': symbol, 'mass_number': A,
            'neutron_endf': f"{endf_root}/neutrons/n-{isotope}.endf",
            'gamma_endf': f"{endf_root}/gammas/g-{isotope}.endf",
            'photoat_endf': f"{endf_root}/photoat/photoat-{element}.endf"}


def get_thermal_info(element: str, molecule: str) -> dict:
    """
    Return the thermal scattering parameters, including the
    S(alpha, beta) ENDF file path, the inelastic and elastic
    reaction numbers, and the number of principal atoms.

    Parameters
    ----------
    element : str, The element symbol.
    molecule : str, The molecule name.

    Returns
    -------
    dict
    """
    info = {}
    if element == "H" and molecule == "H2O":
        file_prefix = "HinH2O"
        info['mti'] = 222
        info['n_atoms'] = 2
    elif element == "H" and molecule == "CH2":
        file_prefix = "HinCH2"
        info['mti'] = 223
        info['mtc'] = 224
        info['n_atoms'] = 2
    elif element == "C" and molecule == "graphite":
        file_prefix = "crystalline-graphite"
        info['mti'] = 229
        info['mtc'] = 230
        info['n_atoms'] = 1
    else:
        raise ValueError("Unrecognized molecule.")

    sab_root = f"{os.environ['ENDF_ROOT']}/thermal_scatt"
    info['sab_endf'] = f"{sab_root}/tsl-{file_prefix}.endf"
    if not os.path.isfile(info['sab_endf']):
        raise FileNotFoundError(
            f"{info['sab_endf']} is not a valid S(alpha, beta) ENDF file."
        )

    return info


def get_group_structure_info(neutron_gs: str, gamma_gs: str) -> dict:
    """
    Return the NJOY group structure inputs.

    Parameters
    ----------
    neutron_gs : str
    gamma_gs : str

    Returns
    -------
    int, The NJOY neutron group set index.
    int, The NJOY gamma group set index.
    str, A path to a custom neutron group structure file.
    str, A path to a custom gamma group structure file.
    """

    info = {'output_directory': None,
            'neutron': {}, 'gamma': {}}

    # define output directory
    outdir = f"output/ENDF-B-VIII-0"
    if neutron_gs:
        outdir = f"{outdir}/{neutron_gs}n"
    if gamma_gs:
        sep = "/" if not neutron_gs else "_"
        outdir = f"{outdir}{sep}{gamma_gs}g"
    info['output_directory'] = outdir

    ##################################################
    # Neutron group structure
    ##################################################

    if neutron_gs:

        # check custom
        if neutron_gs.startswith("custom"):
            info['neutron']['gs_id'] = 1

            # check the custom file
            ngs_file = f"{outdir}/{neutron_gs}n.txt"
            if not os.path.isfile(ngs_file):
                raise FileNotFoundError(f"{ngs_file} is not a valid file.")
            info['neutron']['gs_file'] = ngs_file

        # check lanl
        elif neutron_gs.startswith("lanl"):

            # check for a valid group structure
            valid_opts = [str(g) for g in [30, 70, 80, 187, 618]]
            if not any(g in neutron_gs for g in valid_opts):
                raise ValueError("Invalid LANL neutron group structure.")

            # define group structure id
            info['neutron']['gs_id'] = \
                3 if "30" in neutron_gs else \
                11 if "70" in neutron_gs else \
                13 if "80" in neutron_gs else \
                10 if "187" in neutron_gs else 34

        else:
            raise ValueError(f"Invalid neutron group structure.")

    ##################################################
    # Gamma group structure
    ##################################################

    if gamma_gs:

        # check custom
        if gamma_gs.startswith("custom"):
            info['gamma']['gs_id'] = 1

            # check the custom file
            ggs_file = f"{outdir}/{gamma_gs}g.txt"
            if not os.path.isfile(ggs_file):
                raise FileNotFoundError(f"{ggs_file} is not a valid file.")
            info['gamma']['gs_file'] = ggs_file

        # check lanl
        elif gamma_gs.startswith("lanl"):

            # check for a valid group structure
            valid_opts = [str(g) for g in [12, 24, 48]]
            if not any(g in gamma_gs for g in valid_opts):
                raise ValueError("Invalid LANL gamma group structure.")

            # define group structure id
            info['gamma']['gs_id'] = \
                3 if "12" in gamma_gs else \
                7 if "24" in gamma_gs else 6

        else:
            raise ValueError("Invalid gamma group structure.")

    return info
