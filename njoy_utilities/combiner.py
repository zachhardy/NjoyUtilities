"""
Combines NJOY raw data into comprehensible transport data.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


################################################################################
def build_combined_data(
        raw_njoy_data: dict,
        plot: bool = False
) -> tuple[dict, dict]:
    """
    Combines enjoy raw data into a dictionary of vectors and matrices.
    """

    # -------------------------------------------------- raw njoy data

    group_structures: dict = raw_njoy_data["group_structures"]
    cross_sections: dict = raw_njoy_data["cross_sections"]
    transfer_matrices: dict = raw_njoy_data["transfer_matrices"]

    # -------------------------------------------------- number of groups

    n_gs, g_gs = [], []
    if "neutron" in group_structures:
        n_gs = group_structures["neutron"]
    if "gamma" in group_structures:
        g_gs = group_structures["gamma"]

    G_n = len(n_gs)
    G_g = len(g_gs)
    G = G_n + G_g

    # -------------------------------------------------- problem type

    problem_description = {}
    if n_gs:
        if g_gs:
            problem_description["problem_type"] = "Neutron + Gamma"
        elif not g_gs:
            problem_description["problem_type"] = "Neutron only"
    elif g_gs:
        problem_description["problem_type"] = "Gamma only"
    else:
        problem_description["problem_type"] = "Unknown"
        print(problem_description)
        raise Exception("unknown particle type")

    problem_description["G_n"] = G_n
    problem_description["G_g"] = G_g

    # -------------------------------------------------- total xs

    sig_t = np.zeros(G)

    if "(n,total)" in cross_sections:
        data = cross_sections["(n,total)"]
        for entry in data:
            g, v = entry
            sig_t[G_n - g - 1] = v

    if "(g,total)" in cross_sections:
        data = cross_sections["(g,total)"]
        for entry in data:
            g, v = entry
            sig_t[G - g - 1] = v

    # -------------------------------------------------- heating xs

    sig_heat = np.zeros(G)

    if "(n,heat)" in cross_sections:
        data = cross_sections["(n,heat)"]
        for entry in data:
            g, v = entry
            sig_heat[G_n - g - 1] = v

    if "(g,heat)" in cross_sections:
        data = cross_sections["(g,heat)"]
        for entry in data:
            g, v = entry
            sig_heat[G - g - 1] = v

    # -------------------------------------------------- capture xs

    sig_x = np.zeros(G)

    if "(n,g)" in cross_sections:
        data = cross_sections["(n,g)"]
        for entry in data:
            g, v = entry
            sig_x[G_n - g - 1] = v

    if "(g,x)" in cross_sections:
        data = cross_sections["(g,x)"]
        for entry in data:
            g, v = entry
            sig_x[G - g - 1] = v

    # -------------------------------------------------- scattering xs

    sig_el = np.zeros(G)
    if "(n,elastic)" in cross_sections:
        data = cross_sections["(n,elastic)"]
        for entry in data:
            g, v = entry
            sig_el[G_n - g - 1] = v

    sig_inel = np.zeros(G)
    if "(n,inel)" in cross_sections:
        data = cross_sections["(n,inel)"]
        for entry in data:
            g, v = entry
            sig_inel[G_n - g - 1] = v

    sig_therm = np.zeros(G)
    if "free gas" in cross_sections:
        data = cross_sections["free gas"]
        for entry in data:
            g, v = entry
            sig_therm[G_n - g - 1] = v

    elif "inelastic s(a,b)" in cross_sections:
        data = cross_sections["inelastic s(a,b)"]
        for entry in data:
            g, v = entry
            sig_therm[G_n - g - 1] = v

        if "elastic s(a,b)" in cross_sections:
            data = cross_sections["inelastic s(a,b)"]
            for entry in data:
                g, v = entry
                sig_therm[G_n - g - 1] += v

    sig_nxn = np.zeros(G)
    for nn in range(2, 4 + 1):
        rx = f"(n,{nn:01d}n)"
        if rx in cross_sections:
            data = cross_sections[rx]
            for entry in data:
                g, v = entry
                sig_nxn[G_n - g - 1] += v

    sig_nonel = np.zeros(G)
    if "(g,nonel)" in cross_sections:
        data = cross_sections["(g,nonel)"]
        for entry in data:
            g, v = entry
            sig_nonel[G - g - 1] = v

    sig_coht = np.zeros(G)
    if "(g,coherent)" in cross_sections:
        data = cross_sections["(g,coherent)"]
        for entry in data:
            g, v = entry
            sig_coht[G - g - 1] = v

    sig_incoht = np.zeros(G)
    if "(g,incoherent)" in cross_sections:
        data = cross_sections["(g,incoherent)"]
        for entry in data:
            g, v = entry
            sig_incoht[G - g - 1] = v

    sig_pairprod = np.zeros(G)
    if "(g,pair_production)" in cross_sections:
        data = cross_sections["(g,pair_production)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_pairprod[G - g - 1] += v

    # ------------------------------------------------- auxiliary data

    inv_v = np.zeros(G)
    if "inverse velocity" in cross_sections:
        data = cross_sections["inverse velocity"]
        for entry in data:
            g, v = entry
            inv_v[G_n - g - 1] = v

    e_avg = np.zeros(G)
    for particle in ['neutron', 'gamma']:
        ref = G_n if particle == 'neutron' else G_n + G_g
        if f"average energy ({particle})" in cross_sections:
            data = cross_sections[f"average energy ({particle})"]
            for entry in data:
                g, v = entry
                e_avg[ref - g - 1] = v * 1.0e-6

    # if gamma groups are used but average energies were not parsed
    # use the average of the photon group bounds
    if G_g > 0 and np.sum(e_avg[G_n:]) == 0.0:
        e_bounds = group_structures["gamma"]
        for entry in e_bounds:
            g, elow, ehigh = entry
            e_avg[G - g - 1] = 0.5 * (elow + ehigh)

    # -------------------------------------------------- fission data

    sig_f = np.zeros(G)

    if "(n,fission)" in cross_sections:
        data = cross_sections["(n,fission)"]
        for entry in data:
            g, v = entry
            sig_f[G_n - g - 1] = v

    if "(g,fission)" in cross_sections:
        data = cross_sections["(g,fission)"]
        for entry in data:
            g, v = entry
            sig_f[G - g - 1] = v

    nu_total = np.zeros(G)
    nu_prompt = np.zeros(G)
    nu_delayed = np.zeros(G)

    chi_prompt = {"neutron": np.zeros(G_n),
                  "gamma": np.zeros(G_n)}

    decay_const = {"neutron": [], "gamma": []}
    precursor_fraction = {"neutron": [], "gamma": []}
    chi_delayed = {"neutron": [], "gamma": []}

    for particle in ["neutron", "gamma"]:
        ref = G_n if particle == "neutron" else G

        if f"total nubar ({particle})" in cross_sections:
            data = cross_sections[f"total nubar ({particle})"]
            for entry in data:
                g, v = entry
                nu_total[ref - g - 1] = v

        if f"prompt nubar ({particle})" in cross_sections:
            data = cross_sections[f"prompt nubar ({particle})"]
            for entry in data:
                g, v = entry
                nu_prompt[ref - g - 1] = v

        if f"delayed nubar ({particle})" in cross_sections:
            data = cross_sections[f"delayed nubar ({particle})"]
            for entry in data:
                g, v = entry
                nu_delayed[ref - g - 1] = v

        if f"prompt chi ({particle})" in cross_sections:
            data = cross_sections[f"prompt chi ({particle})"]
            for entry in data:
                g, v = entry
                chi_prompt[particle][G_n - g - 1] = v
            chi_prompt[particle] /= sum(chi_prompt[particle])

        if f"decay constants ({particle})" in cross_sections:
            data = cross_sections[f"decay constants ({particle})"]
            lambdas = np.array([entry[1] for entry in data])
            decay_const[particle] = lambdas

            # get chi delayed
            chi_delayed[particle] = np.zeros((G_n, len(lambdas)))
            data = cross_sections[f"delayed chi ({particle})"]
            for entry in data:
                g, v = entry[0], entry[1:]
                chi_delayed[particle][G_n - g - 1] = v

            # compute yield fractions
            frac = np.sum(chi_delayed[particle], axis=0)
            precursor_fraction[particle] = frac

            # normalize chi delayed
            chi_delayed[particle] /= precursor_fraction[particle]

    # -------------------------------------------------- reconciliation

    # compute the total scattering cross-section
    sig_s = np.zeros(G)
    sig_s[:G_n] = sig_el[:G_n] + sig_inel[:G_n]
    sig_s[G_n:] = sig_coht[G_n:] + sig_incoht[G_n:]

    # the absorption is given by the total minus scattering
    sig_a = sig_t - sig_s

    # -------------------------------------------------- number of moments

    M = 0
    for particle, transfers in transfer_matrices.items():
        for rxn_type, data in transfers.items():
            M = max(M, len(data[0][2])) if data else M

    # -------------------------------------------------- find thermal cutoff

    thermal_rxns = []
    for key in transfer_matrices["neutron"].keys():
        if any(s in key for s in ["free gas", "s(a,b)"]):
            thermal_rxns.append(key)
    if len(thermal_rxns) > 1 and "free gas" in thermal_rxns:
        thermal_rxns.pop(thermal_rxns.index("free gas"))
        transfer_matrices["neutron"].pop("free gas")

    n_therm = -1
    for rxn_type in thermal_rxns:
        data = transfer_matrices["neutron"][rxn_type]
        if data:
            n_therm = max(n_therm, max(max(*row[:2]) for row in data))

    # -------------------------------------------------- utility functions

    ############################################################
    # Lambda-ish to add neutron data
    # Always operates on transfer_mats
    def add_transfer_neutron(data_values, dst, cross_particle=False):
        for data_value in data_values:
            G_ref = G_n if not cross_particle else G_n + G_g
            gprime = G_ref - data_value[0] - 1
            group = G_n - data_value[1] - 1
            num_moms = len(data_value[2])
            for moment in range(0, num_moms):
                val = data_value[2][moment]
                dst[moment][gprime][group] += val

    ############################################################
    # Lambda-ish to add gamma data
    # Always operates on transfer_mats
    def add_transfer_gamma(data_values, dst, cross_particle=False):
        for data_value in data_values:
            G_ref = G_n + G_g if not cross_particle else G_n
            gprime = G_ref - data_value[0] - 1
            group = G_n + G_g - data_value[1] - 1

            # Determine number of moments
            num_moms = len(data_value[2])
            for moment in range(0, num_moms):
                val = data_value[2][moment]
                dst[moment][gprime][group] += val

    # -------------------------------------------------- transfer matrices
    # include all scattering except for thermal

    fission_mat = np.zeros((1, G, G))
    scattering_mats = np.zeros((M, G, G))
    for particle, transfers in transfer_matrices.items():
        for rxn_type, data in transfers.items():

            # -------------------- scattering reactions
            if "fission" not in rxn_type:

                # skip thermal reactions
                skip = ["free_gas", "s(a,b)"]
                if any(s in rxn_type for s in skip):
                    continue

                if particle == "neutron":
                    if "n," in rxn_type:
                        add_transfer_neutron(data, scattering_mats)
                    elif "g," in rxn_type:
                        add_transfer_neutron(data, scattering_mats, True)

                elif particle == "gamma":
                    if "g," in rxn_type:
                        add_transfer_gamma(data, scattering_mats)
                    elif "n," in rxn_type:
                        add_transfer_gamma(data, scattering_mats, True)

            # -------------------- fission reactions
            else:
                if particle == "neutron":
                    if "n," in rxn_type:
                        add_transfer_neutron(data, fission_mat)
                    elif "g," in rxn_type:
                        add_transfer_neutron(data, fission_mat, True)
                elif particle == "gamma":
                    if "n," in rxn_type:
                        add_transfer_gamma(data, fission_mat, True)

    # -------------------------------------------------- apply thermal treatment

    if thermal_rxns and n_therm >= -1:

        # shorthand for neutron-neutron transfers
        nn_transfers = transfer_matrices["neutron"]

        # start and end of the thermal groups
        bgn, end = G_n - n_therm, G_n + 1

        # subtract thermal MT2 from the transfer matrix
        elastic = np.zeros((M, G, G))
        add_transfer_neutron(nn_transfers["(n,elastic)"], elastic)
        elastic = elastic[:, bgn:end, bgn:end]
        scattering_mats[:, bgn:end, bgn:end] -= elastic

        # add thermal treatment to the transfer matrix
        thermal = np.zeros((M, G, G))
        for rxn_type in thermal_rxns:
            add_transfer_neutron(nn_transfers[rxn_type], thermal)
        thermal = thermal[:, bgn:end, bgn:end]
        scattering_mats[:, bgn:end, bgn:end] += thermal

        # apply thermal treatment to toal scattering
        sig_s[bgn:end] = sig_therm[bgn:end]

        # recompute the total cross-section in the thermal region
        # from absorption and thermal scattering
        sig_t[bgn:end] = sig_a[bgn:end] + sig_s[bgn:end]

    # -------------------------------------------------- compute sparsity

    scattering_mats_sparsity = []
    for m in range(0, M):
        mat_non_zeros = []
        for gp in range(0, G):
            non_zeros = []
            for g in range(0, G):
                if np.abs(scattering_mats[m][gp, g]) > 1.0e-16:
                    non_zeros.append(g)
            mat_non_zeros.append(non_zeros)
        scattering_mats_sparsity.append(mat_non_zeros)

    # -------------------------------------------------- plotting

    if plot:
        Atest = np.copy(scattering_mats[0])
        nz = np.nonzero(Atest)
        Atest[nz] = np.log10(Atest[nz]) + 10.0

        plt.figure()
        plt.imshow(Atest, cmap=cm.Greys)
        if G < 6:
            ticks = np.arange(0, G, 1)
        elif 6 <= G < 21:
            ticks = np.arange(0, G, 2)
        elif 21 <= G < 41:
            ticks = np.arange(0, G, 4)
        elif 41 <= G < 81:
            ticks = np.arange(0, G, 8)
        elif 81 <= G < 121:
            ticks = np.arange(0, G, 12)
        elif 121 <= G < 201:
            ticks = np.arange(0, G, 20)
        else:
            ticks = np.arange(0, G, 50)

        plt.xticks(ticks, [str(g) for g in ticks])
        plt.yticks(ticks, [str(g) for g in ticks])
        plt.xlabel("Destination energy group")
        plt.ylabel("Source energy group")
        plt.gca().xaxis.set_ticks_position("top")
        plt.gca().xaxis.set_label_position("top")

        plt.tight_layout()
        plt.savefig("TransferMatrix_NJOY.png")

        if n_gs:
            # plot cross-sections
            fig, ax = plt.subplots(ncols=2, figsize=(8, 4))
            for iax in ax:
                iax.set_xlabel("Energy (MeV)")
                iax.set_ylabel("Cross-Section (b)")

            # general
            change = np.max(sig_t[:G_n]) / np.min(sig_t[:G_n])
            plotter = ax[0].semilogx if change < 100.0 else ax[0].loglog
            plotter(e_avg[:G_n], sig_t[:G_n], label=r"$\sigma_t$")
            plotter(e_avg[:G_n], sig_a[:G_n], label=r"$\sigma_a$")
            plotter(e_avg[:G_n], sig_f[:G_n], label=r"$\sigma_f$")
            plotter(e_avg[:G_n], sig_s[:G_n], label=r"$\sigma_s$")
            ax[0].legend()
            ax[0].grid(True)

            # scattering
            change = np.max(sig_s[:G_n]) / np.min(sig_s[:G_n])
            plotter = ax[1].semilogx if change < 100.0 else ax[1].semilogx
            plotter(e_avg[:G_n], sig_s[:G_n], label=r"$\sigma_s$")
            plotter(e_avg[:G_n], sig_el[:G_n], label=r"$\sigma_{el}$")
            plotter(e_avg[:G_n], sig_inel[:G_n], label=r"$\sigma_{inel}$")
            if sum(sig_therm) > 0.0:
                plotter(e_avg[:G_n], sig_therm[:G_n], label=r"$\sigma_{therm}$")
            ax[1].legend()
            ax[1].grid(True)

            plt.tight_layout()
            plt.savefig("CrossSections_n_NJOY.png")

            # plot fission data
            if np.linalg.norm(sig_f) > 1.0e-20:
                fig, ax = plt.subplots(ncols=2, figsize=(8, 4))
                for i, iax in enumerate(ax):
                    iax.set_xlabel("Energy (MeV)")
                    iax.set_ylabel("$\chi$ (MeV$^{-1}$)" if i == 0 else
                                   "Neutrons Per Fission")

                ax[0].semilogx(e_avg[:G_n], chi_prompt["neutron"])
                ax[0].grid(True)

                line0, = ax[1].semilogx(e_avg[:G_n], nu_total[:G_n], "b")
                line1, = ax[1].semilogx(e_avg[:G_n], nu_prompt[:G_n], "r")

                twin_ax = ax[1].twinx()
                line2, = twin_ax.semilogx(e_avg[:G_n], nu_delayed[:G_n], "g")
                twin_ax.set_ylabel("Delayed Neutrons Per Fission")

                lines = [line0, line1, line2]
                labels = ["Total", "Prompt", "Delayed"]
                ax[1].legend(lines, labels)
                ax[1].grid(True)

                plt.tight_layout()
                plt.savefig("FissionData_n_NJOY.png")

        if g_gs:
            # plot cross-sections
            plt.figure()
            plt.xlabel("Energy (MeV)")
            plt.ylabel("Cross-Section (b)")
            plt.semilogx(e_avg[G_n:], sig_t[G_n:], label=r"$\sigma_t$")
            plt.semilogx(e_avg[G_n:], sig_a[G_n:], label=r"$\sigma_a$")
            plt.semilogx(e_avg[G_n:], sig_x[G_n:], label=r"$\sigma_x$")
            plt.semilogx(e_avg[G_n:], sig_s[G_n:], label=r"$\sigma_s$")
            plt.legend()
            plt.grid(True)

            plt.tight_layout()

            plt.savefig("CrossSections_g_NJOY.png")

            # plot fission data
            if np.linalg.norm(sig_f) > 1.0e-20 and n_gs:

                fig, ax = plt.subplots(ncols=2, figsize=(8, 4))
                for i, iax in enumerate(ax):
                    iax.set_xlabel("Energy (MeV)")
                    iax.set_ylabel("$\chi$ (MeV$^{-1}$)" if i == 0 else
                                   r"Neutrons Per Fission")

                ax[0].semilogx(e_avg[:G_n], chi_prompt["gamma"])
                ax[0].grid(True)

                line0, = ax[1].semilogx(e_avg[G_n:], nu_total[G_n:], "b")
                line1, = ax[1].semilogx(e_avg[G_n:], nu_prompt[G_n:], "r")

                twin_ax = ax[1].twinx()
                line2, = twin_ax.semilogx(e_avg[G_n:], nu_delayed[G_n:], "g")
                twin_ax.set_ylabel("Delayed Neutrons Per Fission")

                lines = [line0, line1, line2]
                labels = ["Total", "Prompt", "Delayed"]
                ax[1].legend(lines, labels)
                ax[1].grid(True)

                plt.tight_layout()
                plt.savefig("FissionData_g_NJOY.png")

        plt.show()

    ###########################################################################
    # Build return data
    ###########################################################################

    return_data = {"neutron_gs": n_gs,
                   "gamma_gs": g_gs,
                   "sigma_t": sig_t,
                   "sigma_a": sig_a,
                   "sigma_s": sig_s,
                   "sigma_s_el": sig_el,
                   "sigma_s_inel": sig_inel,
                   "sigma_s_thermal": sig_therm,
                   "sigma_f": sig_f,
                   "sigma_heat": sig_heat,
                   "nu_total": nu_total,
                   "nu_prompt": nu_prompt,
                   "nu_delayed": nu_delayed,
                   "chi_prompt": chi_prompt,
                   "chi_delayed": chi_delayed,
                   "decay_constants": decay_const,
                   "precursor_fraction": precursor_fraction,
                   "inv_velocity": inv_v,
                   "avg_energy": e_avg,
                   "scattering_matrices": scattering_mats,
                   "scattering_matrices_sparsity": scattering_mats_sparsity,
                   "fission_matrices": fission_mat}

    return return_data, problem_description
