"""Writes cross-sections to file"""
import os.path
import warnings

import numpy as np


######################################################################
def write_xs_file(data, full_path, problem_description):
    with open(full_path, "w") as xsf:

        # --------------------------------------------------
        # Write header
        # --------------------------------------------------

        xsf.write("#========== Problem Description ==========" + "\n")

        n_group = problem_description["G_n"]
        g_group = problem_description["G_g"]

        xsf.write(f"# Isotope: {problem_description['isotope']}\n")
        xsf.write(f"# Problem type: {problem_description['problem_type']}\n")

        if n_group > 0:
            xsf.write(f"# Neutron group structure: {n_group} groups\n")
        if g_group > 0:
            xsf.write(f"# Gamma group structure: {g_group} groups\n")
        xsf.write("\n")

        # --------------------------------------------------
        # Get data
        # --------------------------------------------------

        sig_t = data["sigma_t"]
        sig_a = data["sigma_a"]
        sig_s = data["sigma_s"]
        sig_f = data["sigma_f"]
        sig_heat = data["sigma_heat"]
        nu_total = data["nu_total"]
        nu_prompt = data["nu_prompt"]
        nu_delayed = data["nu_delayed"]
        chi_prompt = data["chi_prompt"]
        chi_delayed = data["chi_delayed"]
        decay_const = data["decay_constants"]
        precursor_fraction = data["precursor_fraction"]
        inv_velocity = data["inv_velocity"]
        E_avg = data["avg_energy"]

        scattering_mats = data["scattering_matrices"]
        scattering_mats_nonzeros = data["scattering_matrices_sparsity"]

        fission_mats = data["fission_matrices"]
        fission_mats_nonzeros = data["fission_matrices_sparsity"]

        neutron_gs = data["neutron_gs"]
        gamma_gs = data["gamma_gs"]

        # ------------------------------------------------------------
        # Write general info
        # ------------------------------------------------------------

        G = n_group + g_group
        M = len(scattering_mats)
        J = len(decay_const['neutron'])

        xsf.write("# Output\n")
        xsf.write(f"NUM_GROUPS {G}\n")
        xsf.write(f"NUM_MOMENTS {M}\n")
        xsf.write(f"NUM_PRECURSORS {J}\n")
        xsf.write("\n")

        # ------------------------------------------------------------
        # Write group structures
        # ------------------------------------------------------------

        if neutron_gs:
            # print(n_group, len(neutron_gs))
            # neutron_gs = np.flip(neutron_gs)
            xsf.write("NEUTRON_GS_BEGIN  # MeV\n")
            for n in range(n_group):
                xsf.write(f"{n:<4d} ")
                xsf.write(f"{neutron_gs[n][1]:<12.8g} ")
                xsf.write(f"{neutron_gs[n][2]:<12.8g}\n")
            xsf.write("NEUTRON_GS_END\n\n")

        if gamma_gs:
            # print(g_group, len(gamma_gs))
            # gamma_gs = np.flip(gamma_gs)
            xsf.write("GAMMA_GS_BEGIN  # MeV\n")
            for g in range(g_group):
                xsf.write(f"{g:<4d} ")
                xsf.write(f"{gamma_gs[g][1]:<12.8g} ")
                xsf.write(f"{gamma_gs[g][2]:<12.8g}\n")
            xsf.write("GAMMA_GS_END\n\n")

        xsf.write("AVG_ENERGY_BEGIN  # MeV\n")
        for g in range(G):
            xsf.write(f"{g:<4d} {E_avg[g]:<12.8g}\n")
        xsf.write("AVG_ENERGY_END\n\n")

        # ------------------------------------------------------------
        # Write reaction data
        # ------------------------------------------------------------

        xsf.write("SIGMA_T_BEGIN\n")
        for g in range(G):
            xsf.write(f"{g:<4d} {sig_t[g]:<g}\n")
        xsf.write("SIGMA_T_END\n\n")

        xsf.write("SIGMA_A_BEGIN\n")
        for g in range(G):
            xsf.write(f"{g:<4d} {sig_a[g]:<g}\n")
        xsf.write("SIGMA_A_END\n\n")

        xsf.write("SIGMA_S_BEGIN\n")
        for g in range(G):
            xsf.write(f"{g:<4d} {sig_s[g]:<g}\n")
        xsf.write("SIGMA_S_END\n\n")

        xsf.write("SIGMA_HEAT_BEGIN\n")
        for g in range(G):
            xsf.write(f"{g:<4d} {sig_heat[g]:<g}\n")
        xsf.write("SIGMA_HEAT_END\n\n")

        if np.linalg.norm(inv_velocity) > 1.0e-20:
            xsf.write("VELOCITY_BEGIN  # cm/sh\n")
            for g in range(n_group):
                velocity = 1.0e-6 / inv_velocity[g]
                xsf.write(f"{g:<4d} {velocity:<g}\n")
            xsf.write("VELOCITY_END\n\n")

        # ------------------------------------------------------------
        # Write fission data
        # ------------------------------------------------------------

        if np.linalg.norm(sig_f) > 1.0e-20:
            xsf.write("SIGMA_F_BEGIN\n")
            for g in range(G):
                xsf.write(f"{g:<4d} {sig_f[g]:<g}\n")
            xsf.write("SIGMA_F_END\n\n")

        if np.linalg.norm(nu_total) > 1.0e-20:
            xsf.write("NU_BEGIN\n")
            for g in range(G):
                xsf.write(f"{g:<4d} {nu_total[g]:<g}\n")
            xsf.write("NU_END\n\n")

        if np.linalg.norm(nu_prompt) > 1.0e-20:
            xsf.write("NU_PROMPT_BEGIN\n")
            for g in range(G):
                xsf.write(f"{g:<4d} {nu_prompt[g]:<12g}\n")
            xsf.write("NU_PROMPT_END\n\n")

        if J > 0:
            xsf.write("NU_DELAYED_BEGIN\n")
            for g in range(G):
                xsf.write(f"{g:<4d} {nu_delayed[g]:<g}\n")
            xsf.write("NU_DELAYED_END\n\n")

        # this is the neutron induced spectrum
        if np.linalg.norm(chi_prompt["neutron"]) > 1.0e-20:
            xsf.write("CHI_PROMPT_BEGIN\n")
            for g in range(n_group):
                xsf.write(f"{g:<4d} {chi_prompt['neutron'][g]:<g}\n")
            xsf.write("CHI_PROMPT_END\n\n")

        if J > 0:
            xsf.write("CHI_DELAYED_BEGIN\n")
            for g in range(G):
                for j in range(J):
                    xsf.write("G_PRECURSORJ_VAL ")
                    xsf.write(f"{g:<4d} {j:<3d} "
                              f"{chi_delayed['neutron'][g][j]:<g}\n")
            xsf.write("CHI_DELAYED_END\n\n")

            xsf.write("PRECURSOR_LAMBDA_BEGIN\n")
            for j in range(J):
                xsf.write(f"{j:<4d} "
                          f"{decay_const['neutron'][j]:<12g}\n")
            xsf.write("PRECURSOR_LAMBDA_END\n\n")

            xsf.write("PRECURSOR_YIELD_BEGIN\n")
            for j in range(J):
                xsf.write(f"{j:<4d} "
                          f"{precursor_fraction['neutron'][j]:<g}\n")
            xsf.write("PRECURSOR_YIELD_END\n\n")

        # ------------------------------------------------------------
        # Write transfer matrix data
        # ------------------------------------------------------------

        xsf.write("TRANSFER_MOMENTS_BEGIN\n")
        for m in range(M):
            xsf.write(f"# l = {m}\n")
            for gp in range(G):
                for g in scattering_mats_nonzeros[m][gp]:
                    xsf.write("M_GPRIME_G_VAL ")
                    xsf.write(f"{m:<3d} {gp:<4d} {g:<4d} ")
                    xsf.write(f"{scattering_mats[m][gp][g]:<g}\n")
        xsf.write("TRANSFER_MOMENTS_END\n\n")

        # xsf.write("PRODUCTION_MATRIX_BEGIN\n")
        # for gp in range(G):
        #     for g in fission_mats_nonzeros[0][gp]:
        #         xsf.write("M_GPRIME_G_VAL ")
        #         xsf.write(f"0 {gp:<4d} {g:<4d} ")
        #         xsf.write(f"{fission_mats[0][gp][g]:<g}\n")
        # xsf.write("PRODUCTION_MATRIX_END\n")
