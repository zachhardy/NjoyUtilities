import os
import warnings
import argparse


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.MetavarTypeHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


def get_material_info(material):
    """
    Return the ENDF paths for a particular material input.

    Parameters
    ----------
    material : str
        The material input formatted via <Z>-<element>-<A>-<molecule>.

    Returns
    -------
    dict
        A dictionary containing the ENDF file paths.
    """
    entries = material.split("-")
    symbol = entries[1]
    z, a = int(entries[0]), int(entries[2])

    isotope = "_".join([str(z).zfill(3), symbol, str(a).zfill(3)])
    element = "_".join([str(z).zfill(3), symbol, "000"])
    molecule = None if len(entries) < 4 else entries[3]

    endf_root = os.getenv("ENDF_ROOT")

    neutron_endf = f"{endf_root}/neutrons/n-{isotope}.endf"
    if not os.path.isfile(neutron_endf):
        warnings.warn("No neutron ENDF file found.")
        neutron_endf = None

    gamma_endf = f"{endf_root}/gammas/g-{isotope}.endf"
    if not os.path.isfile(gamma_endf):
        warnings.warn("No gamma ENDF file found.")
        gamma_endf = None

    photoat_endf = f"{endf_root}/photoat/photoat-{element}.endf"
    if not os.path.isfile(photoat_endf):
        warnings.warn("No photo-atomic ENDF file found.")
        photoat_endf = None

    return {"isotope": isotope,
            "molecule": molecule,
            "atomic_number": z,
            "symbol": symbol,
            "mass_number": a,
            "neutron_endf": neutron_endf,
            "gamma_endf": gamma_endf,
            "photoat_endf": photoat_endf}


def get_thermal_info(element: str, molecule: str) -> dict:
    """
    Return the thermal scattering parameters, including the S(alpha, beta) ENDF file path, the
    inelastic and elastic reaction numbers, and the number of principal atoms.

    Parameters
    ----------
    element : str
        The element symbol.
    molecule : str
        The molecule name.

    Returns
    -------
    dict
        Thermal scattering information.
    """
    info = {}
    if element == "H" and molecule == "H2O":
        file_prefix = "HinH2O"
        info["mti"] = 222
        info["n_atoms"] = 2
    elif element == "H" and molecule == "CH2":
        file_prefix = "HinCH2"
        info["mti"] = 223
        info["mtc"] = 224
        info["n_atoms"] = 2
    elif element == "C" and molecule == "graphite":
        file_prefix = "crystalline-graphite"
        info["mti"] = 229
        info["mtc"] = 230
        info["n_atoms"] = 1
    else:
        raise ValueError("Unrecognized molecule.")

    sab_root = f"{os.environ['ENDF_ROOT']}/thermal_scatt"
    info["sab_endf"] = f"{sab_root}/tsl-{file_prefix}.endf"
    if not os.path.isfile(info["sab_endf"]):
        raise FileNotFoundError(
            f"{info['sab_endf']} is not a valid S(alpha, beta) ENDF file."
        )

    return info


def get_group_structure_info(group_structures: list[str]) -> dict:
    """
    Return the NJOY group structure inputs.

    Parameters
    ----------
    group_structures : list of str

    Returns
    -------
    int
        The NJOY neutron group set index.
    int
        The NJOY gamma group set index.
    str
        A path to a custom neutron group structure file.
    str
        A path to a custom gamma group structure file.
    """

    info = {"outdir": "output/ENDF-B-VIII-0",
            "neutron": {}, "gamma": {}}

    # Find present group structures
    ngs, ggs = None, None
    for gs in group_structures:
        if gs.endswith("n"):
            ngs = gs
        elif gs.endswith("g"):
            ggs = gs
        else:
            raise AssertionError(
                "Unrecognized group structure specification."
            )

    # Define the output directory
    if ngs and ggs:
        outdir = f"{ngs}_{ggs}"
    elif ngs and not ggs:
        outdir = f"{ngs}"
    elif ggs and not ngs:
        outdir = f"{ggs}"
    else:
        raise AssertionError("No neutron or gamma group structure found.")
    info["outdir"] = os.path.join(info["outdir"], outdir)

    # Define the neutron NJOY input
    if ngs:

        # Custom group structures
        if ngs.startswith("custom"):
            info["neutron"]["gs_id"] = 1

            # Get number of groups
            loc = ngs.rfind("custom") + 6
            info["neutron"]["n_groups"] = int(ngs[loc:-1])

            # Ensure there is a custom file
            n_gs_file = os.path.join(info["outdir"], f"{ngs}.txt")
            if not os.path.isfile(n_gs_file):
                raise FileNotFoundError(f"{n_gs_file} is not a valid file.")
            info["neutron"]["gs_file"] = n_gs_file

        # LANL group structures
        elif ngs.startswith("lanl"):

            # Check for a valid group structure
            valid_opts = [str(g) for g in [30, 70, 80, 187, 618]]
            if not any(g in ngs for g in valid_opts):
                raise ValueError("Invalid LANL neutron group structure.")

            # Define group structure id
            info["neutron"]["gs_id"] = \
                3 if "30" in ngs else \
                    11 if "70" in ngs else \
                        13 if "80" in ngs else \
                            10 if "187" in ngs else 34

            # get number of groups
            loc = ngs.find("lanl") + 4
            info["neutron"]["n_groups"] = int(ngs[loc:-1])

        else:
            raise ValueError("Invalid neutron group structure.")

    # Define gamma NJOY input
    if ggs:

        # Custom group structures
        if ggs.startswith("custom"):
            info["gamma"]["gs_id"] = 1

            # get number of groups
            loc = ggs.rfind("custom") + 6
            info["gamma"]["n_groups"] = int(ggs[loc:-1])

            # check the custom file
            g_gs_file = os.path.join(info["outdir"], f"{ggs}.txt")
            if not os.path.isfile(g_gs_file):
                raise FileNotFoundError(f"{g_gs_file} is not a valid file.")
            info["gamma"]["gs_file"] = g_gs_file

        # LANL group structures
        elif ggs.startswith("lanl"):

            # check for a valid group structure
            valid_opts = [str(g) for g in [12, 24, 48]]
            if not any(g in ggs for g in valid_opts):
                raise ValueError("Invalid LANL gamma group structure.")

            # define group structure id
            info["gamma"]["gs_id"] = \
                3 if "12" in ggs else \
                    7 if "24" in ggs else 6

            # get number of groups
            loc = ggs.rfind("lanl") + 4
            info["gamma"]["n_groups"] = int(ggs[loc:-1])

        else:
            raise ValueError("Invalid gamma group structure.")

    return info
