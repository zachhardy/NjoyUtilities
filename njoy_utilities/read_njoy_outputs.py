"""
Contains all the necessary functionality to read NJOY
cross-section files.
"""
import sys
from typing import Union

RawXS = list[tuple[int, float]]
RawDelayedChi = list[tuple[int, list[float]]]
RawTransferMatrix = list[tuple[int, int, list[float]]]


######################################################################
def print_line(line_num: int, line: str) -> None:
    """ Prints a file line without the line feed. """
    print("Line: ", line_num, line, end='')


######################################################################
def string_to_float(string: str) -> float:
    """
    Convert a float string to a float.

    Parameters
    ----------
    string : str

    Returns
    -------
    float
    """
    number = string.replace("+", "E+")
    if "-" in number:
        if number.count("-") == 2:
            pos = number.rfind("-")
            number = f"{number[:pos]}E-{number[pos + 1:]}"
        elif number.find("-") != 0:
            number = number.replace("-", "E-")
    return float(number)


######################################################################
def process_group_structure(
        n: int,
        lines: list[str]
) -> list[list]:
    """
    Processes a group structure.

    Parameters
    ----------
    n : int, The current line number.
    lines : list[str], The list of lines in the file.

    Returns
    -------
    list[list]
    """
    # next line until first entry is found
    words = lines[n].split()
    while not words or words[0] != "1":
        n += 1
        words = lines[n].split()

    # Go through group structure
    group_struct = []
    while words:
        # parse the current line
        group = int(words[0]) - 1
        e_low, e_high = float(words[1]), float(words[3])
        group_struct.append([group, e_low, e_high])

        # go to next line
        n += 1
        words = lines[n].split()
    return group_struct


######################################################################
def process_cross_section(
        n: int,
        lines: list[str],
        header_size: int = 5,
        line_incr: int = 1,
        wt_spectrum: bool = False
) -> RawXS:
    """
    Reads 1D neutron cross-sections from the lines.

    Parameters
    ----------
    n : int, The current line number.
    lines : list[str], The list of lines in the file.
    header_size : int
    line_incr : int
    wt_spectrum : bool

    Returns
    -------
    A table containing the group wise xs.
    """
    # go to first line of data
    n += 1 if "particle emission" in lines[n + 1] else 0
    n += header_size
    words = lines[n].split()

    # if no data, return empty
    if len(words) == 0:
        return []

    # parse the cross-section
    xs = []
    while len(words) >= 2:
        value = string_to_float(words[1])
        if wt_spectrum and words[0] == "flx":
            xs.append([value])
        else:
            group = int(words[0]) - 1
            xs.append((group, value))

        n += line_incr
        words = lines[n].split()
    return xs


######################################################################
def process_prompt_chi(
        n: int,
        lines: list[str],
        header_size: int = 4,
        line_incr: int = 1
) -> RawXS:
    """
    Reads 1D prompt fission spectrum from the lines.

    Parameters
    ----------
    n : int, The current line number.
    lines : list[str], The list of lines in the file.
    header_size : int
    line_incr : int
      
    Returns
    -------
    A table containing the group-wise xs.
    """
    # go to first line of data
    n += 2 if "spectrum constant" in lines[n + 2] else 0
    n += header_size
    words = lines[n].split()

    # parse the data
    spec = []
    while len(words) >= 2:

        group = int(words[0]) - 1
        for number_word in words[1:]:
            spec.append([group, string_to_float(number_word)])
            group += 1

        n += line_incr
        words = lines[n].split()
    return spec


######################################################################
def process_decay_constants(
        n: int,
        lines: list[str]
) -> list[float]:
    """
    Reads delayed neutron decay constants from delayed chi.

    Parameters
    ----------
    n : int, File line number before data starts.
    lines : list[str], The list of lines in the file.

    Returns
    -------
    A table containing delayed neutron group decay constants.
    """
    # go to first line of data
    n += 5 if "spectrum constant" in lines[n + 2] else 3
    words = lines[n].split()
    if words[0] != "group":
        raise AssertionError("Unexpected line encountered.")

    # parse the data
    return [StringToFloat(val) for val in words[1:]]


######################################################################
def process_delayed_chi(
        n: int,
        lines: list[str]
) -> RawDelayedChi:
    """
    Reads the delayed neutron spectra from the lines.

    Parameters
    ----------
    n : int, File line number before data starts.
    lines : list[str], The list of lines in the file.

    Returns
    -------
    A table containing the group wise delayed  neutron precursor
        spectrum coefficients.
    """
    # go to first line of data
    n += 7 if "spectrum constant" in lines[n + 2] else 5
    words = lines[n].split()

    # parse the data
    matrix = []
    while len(words) >= 2:
        group = int(words[0]) - 1
        vals = [StringToFloat(word) for word in words[1:]]
        matrix.append((group, vals))

        n += 1
        words = lines[n].split()
    return matrix


######################################################################
def process_transfer_matrix(
        n: int,
        lines: list[str],
        overflow: bool = False
) -> RawTransferMatrix:
    """
    Reads transfer matrices from the lines.

    Parameters
    ----------
    n : int, File line number before data starts.
    lines : list[str], The list of lines in the file.
    overflow : bool, A flag for overflow lines

    Returns
    -------
    A table containing the group wise transfer coefficients.
    """
    # go to the first line of informative data
    skip = 1 if "mf26" in lines[n] else 0
    skip += 1 if "particle emission" in lines[n + 1] else 0
    skip += 4 if "spectrum constant" in lines[n + 2] else 2
    n += skip

    # determine matrix type
    if "legendre" in lines[n]:
        matrix_type = "legendre"
    elif "isotropic" in lines[n]:
        matrix_type = "isotropic"
    else:
        raise ValueError("Unexpected line encountered.")
    n += 3
    words = lines[n].split()

    # parse the data
    matrix = []
    while len(words) >= 2:

        # In fission matrices, there are lines with
        # `spec` and `prod`, and in some gamma matrices there are lines
        # with `xs` and `heat`. The try/except block ensures that these
        # lines are safely passed over without having to explicitly
        # hard-code the words to skip
        try:
            gprime = int(words[0]) - 1
            group = int(words[1]) - 1
            vals = [string_to_float(val) for val in words[2:]]
            if overflow:
                n += 1
                words = lines[n].split()
                vals.extend([string_to_float(val) for val in words])

            if matrix_type == "legendre":
                matrix.append((gprime, group, vals))
            else:
                for g in range(len(vals)):
                    matrix.append((gprime, group + g, [vals[g]]))

            n += 1
            words = lines[n].split()
        except (ValueError, IndexError):
            n += 1
            words = lines[n].split()

    return matrix


######################################################################
def read_njoy_file(
        njoy_filename: str = "output",
        verbose: bool = False
) -> dict:
    """
    Read an NJOY output file.

    Parameters
    ----------
    njoy_filename : str, Name of the NJOY file to process.
    verbose : bool, default False

    Returns
    -------
    Returns a complex dictionary of raw data.
    """

    njoy_raw_data = {}
    group_structures = {}
    cross_sections = {}
    transfer_matrices = {"neutron": {}, "gamma": {}}
    weight_spectrum = {}

    # flag_run_processed             = False
    flag_gamma_structure_processed = False

    # Read the file
    with open(njoy_filename, "r") as njoy_file:

        line_num = -1
        file_lines = njoy_file.readlines()
        while line_num < len(file_lines) - 1:
            if not file_lines[line_num + 1]:
                line_num += 1
                continue

            line_num += 1
            line = file_lines[line_num].strip()
            words = file_lines[line_num].split()
            num_words = len(words)

            ##################################################
            # group-structures
            ##################################################

            if "sigma zeroes" in line:
                group_structures["neutron"] = process_group_structure(
                    line_num, file_lines
                )

            if "gamma group structure......" in line:
                if not flag_gamma_structure_processed:
                    group_structures["gamma"] = process_group_structure(
                        line_num, file_lines
                    )
                    flag_gamma_structure_processed = True

            ##################################################
            # reactions
            ##################################################

            if line.startswith("for mf") and "mt" in line:

                # define file and reaction number
                if words[1] == "mf":
                    mf = int(words[2])
                    mt = int(words[5] if words[4] == "mt" else
                             words[4].strip("mt"))
                else:
                    mf = int(words[1].strip("mf"))
                    mt = int(words[3].strip("mt"))

                ######################################################
                # cross sections
                ######################################################

                if line.endswith("cross section"):

                    # the reaction type
                    rxn_type = words[num_words - 3]

                    # handle total interaction differently
                    if mf == 3 and mt == 1:
                        cross_sections[rxn_type] = process_cross_section(
                            line_num, file_lines, line_incr=2
                        )

                        # this is to get the weighting spectrum
                        weight_spectrum["neutron"] = process_cross_section(
                            line_num, file_lines, wt_spectrum=True
                        )

                    # all other reactions
                    else:
                        cross_sections[rxn_type] = \
                            process_cross_section(line_num, file_lines)

                ######################################################
                # other cross-sections
                ######################################################

                mtnames = ["inverse velocity", "average energy",
                           "free gas", "inelastic s(a,b)", "elastic s(a,b)",
                           "total nubar", "prompt nubar", "delayed nubar"]

                if mf == 3:
                    for mtname in mtnames:
                        if mtname in line:
                            rxn_type = line[line.find(mtname):]
                            cross_sections[rxn_type] = process_cross_section(
                                line_num, file_lines
                            )

                ######################################################
                # prompt/delayed fission data
                ######################################################

                # prompt fission spectrum
                if mf == 5 and "prompt chi" in line:
                    rxn_type = line[line.find("prompt chi"):]
                    cross_sections[rxn_type] = process_prompt_chi(
                        line_num, file_lines
                    )

                if mf == 5 and "delayed_chi" in line:
                    cross_sections["decay constants"] = \
                        process_decay_constants(line_num, file_lines)

                    rxn_type = line[line.find("delayed chi"):]
                    cross_sections[rxn_type] = process_delayed_chi(
                        line_num, file_lines
                    )

                ######################################################
                # neutron transfer matrix data
                ######################################################

                # caveat:
                #   the pp values from transfer(mf26) = 2x the pp
                #   from the xsec(mf23)
                # caveat:
                #   the n,2n values from transfer(mf8/mt16) = 2x the n,2n
                #   from the xsec(mf3/mt16)
                # caveat:
                #   the n,3n values from transfer(mf8/mt17) = 3x the n,3n
                #   from the xsec(mf3/mt17
                if line.endswith("matrix"):
                    particle_type = words[num_words - 2]
                    rxn_type = words[num_words - 3]

                    if "free gas" in line:
                        rxn_type = "free gas"
                    if "s(a,b)" in line:
                        rxn_type = f"{words[num_words - 4]} s(a,b)"

                    transfer_matrices[particle_type][rxn_type] = \
                        process_transfer_matrix(line_num, file_lines)

                ######################################################
                # photo-atomic cross-section data
                ######################################################

                # Total photon interaction
                if mf == 23 and mt == 501:
                    cross_sections["(g,total)"] = process_cross_section(
                        line_num, file_lines, header_size=4
                    )

                # Photon coherent scattering
                if mf == 23 and mt == 502:
                    cross_sections["(g,coherent)"] = process_cross_section(
                        line_num, file_lines, header_size=4
                    )

                # Photon incoherent scattering
                if mf == 23 and mt == 504:
                    cross_sections["(g,incoherent)"] = process_cross_section(
                            line_num, file_lines, header_size=4
                        )

                # 515: Pair production, electron field
                # 517: Pair production, nuclear field
                # 516: Pair production; sum of MT=515, 517.
                if mf == 23 and mt == 516:
                    cross_sections["(g,pair_production)"] = \
                        process_cross_section(
                            line_num, file_lines, header_size=4
                        )

                # Photoelectric absorption
                if mf == 23 and mt == 522:
                    cross_sections["(g,abst)"] = process_cross_section(
                        line_num, file_lines, header_size=4
                    )

                if mf == 23 and mt == 525:
                    cross_sections["(g,heat)"] = process_cross_section(
                        line_num, file_lines, header_size=4
                    )

                ######################################################
                # photo-atomic transfer matrix data
                ######################################################

                if mf == 26 and mt == 502:
                    transfer_matrices["gamma"]["(g,coherent)"] = \
                        process_transfer_matrix(
                            line_num, file_lines, overflow=True
                        )

                # we do not need to save the xsec(g), it is just the sum_k
                # xs(g->k) the incoh heat value is not saved either,
                # the total heat xs (mf23/mt525) is the sum of the heat from
                # inch(mf23/mt504) + abst(mf23/mt522) + pp(mf23/mt516)
                if mf == 26 and mt == 504:
                    transfer_matrices["gamma"]["(g,incoherent)"] = \
                        process_transfer_matrix(
                            line_num, file_lines, overflow=True
                        )

                # caveat:
                #   the pp values from transfer(mf26) = 2x the pp
                #   from the xsec(mf23)
                if mf == 26 and mt == 516:
                    transfer_matrices['gamma']['(g,pair_production)'] = \
                        process_transfer_matrix(line_num, file_lines)

    njoy_raw_data["group_structures"] = group_structures
    njoy_raw_data["cross_sections"] = cross_sections
    njoy_raw_data["transfer_matrices"] = transfer_matrices

    if verbose:
        print("Cross-sections extracted:")
        xss = njoy_raw_data["cross_sections"]
        for k in xss:
            print(f"\t{k}")

        print("Transfer matrices extracted")
        mats = njoy_raw_data["transfer_matrices"]
        for ptype in mats:
            print(f"\t{ptype}")
            for k in mats[ptype]:
                print(f"\t\t{k}")
    return njoy_raw_data
