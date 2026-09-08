/*******************************************************************************
 *
 * Test for the split-even telescopic estimator vs standard estimator
 * Copyright (c) 2026 P. Butti
 * All rights reserved.
 *
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *
 * Tests that the split-even estimator output is consistent with two
 * independent standard estimators run at masses m1 and m2.
 *
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

// /* Deterministic spin-diluted point source at (tau,0,0,0), colour 0.
//    Writes exactly 4 spinor fields -> drop-in for the 4-vector pipeline.
//    No RNG, so it is identical across every call by construction. */
// void create_point_source_equal_atau(spinor_field *source, int tau) {
//     for (int i = 0; i < 4; ++i) {
// #ifdef WITH_GPU
//         zero_spinor_field_cpu(&source[i]);
// #endif
//         zero_spinor_field(&source[i]);
//     }
//     if (zerocoord[0] <= tau && tau < zerocoord[0] + T &&
//         COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0) {
//         int ix = ipt(tau - zerocoord[0], 0, 0, 0);
//         for (int beta = 0; beta < 4; ++beta) {
//             _FIELD_AT(&source[beta], ix)->c[beta].c[0] = 1.0;  /* spin beta, colour 0 */
//         }
//     }
// #ifdef WITH_GPU
//     for (int i = 0; i < 4; ++i) { copy_to_gpu(&source[i]); }
// #endif
//     for (int i = 0; i < 4; ++i) {
//         start_sendrecv_spinor_field(&source[i]);
//         complete_sendrecv_spinor_field(&source[i]);
//     }
// }


/* ------------------------------------------------------------------ *
 * All physics parameters are set here directly.                       *
 * Only the lattice geometry and RNG are read from the input file.     *
 * ------------------------------------------------------------------ */

/* Two masses: m1 < m2 */
// static const double mass_m1 = -0.73095783;
// static const double mass_m2 = -0.72;

static const double mass_m1 = 0.4;
static const double mass_m2 = 0.5;

/* Estimator parameters */
static const double precision   = 1e-9;
static const int    nhits       = 20;
static const int    source_type = 6;   /* time + spin dilution + gaussian smearing */
static const int    n_mom       = 1;
static const int    n_smr       = 0;
static const double alpha       = 0.2;

/* No HYP smearing */
#define HYP_WEIGHTS NULL

char input_filename[256] = "input_file_test_spliteven";

int main(int argc, char *argv[]) {

    Timer clock;
    timer_set(&clock);

    /* setup process id and communications */
    setup_process(&argc, &argv);
    setup_gauge_fields();

    init_BCs(NULL);

    /* ---------------------------------------------------------------- *
     * Summary of test parameters                                        *
     * ---------------------------------------------------------------- */
    lprintf("MAIN", 0, "============================================================\n");
    lprintf("MAIN", 0, "  Split-even vs standard estimator test w/ full point sources\n");
    lprintf("MAIN", 0, "============================================================\n");
    lprintf("MAIN", 0, "  m1              = %f\n", mass_m1);
    lprintf("MAIN", 0, "  m2              = %f\n", mass_m2);
    lprintf("MAIN", 0, "  precision       = %e\n", precision);
    lprintf("MAIN", 0, "  nhits           = %d\n", nhits);
    lprintf("MAIN", 0, "  source_type     = %d\n", source_type);
    lprintf("MAIN", 0, "  n_mom           = %d\n", n_mom);
    lprintf("MAIN", 0, "  n_smr           = %d\n", n_smr);
    lprintf("MAIN", 0, "  alpha           = %f\n", alpha);
    lprintf("MAIN", 0, "  HYP smearing    = disabled\n");
    lprintf("MAIN", 0, "============================================================\n");

    /* ---------------------------------------------------------------- *
     * Gauge field: unit configuration                                   *
     * ---------------------------------------------------------------- */
    unit_u(u_gauge);
    represent_gauge_field();

    lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
    full_plaquette();

    /* ---------------------------------------------------------------- *
     * Section 1: Standard estimator at m1                              *
     * ---------------------------------------------------------------- */
    lprintf("MAIN", 0, "------------------------------------------------------------\n");
    lprintf("MAIN", 0, "  STANDARD ESTIMATOR at m1 = %f\n", mass_m1);
    lprintf("MAIN", 0, "------------------------------------------------------------\n");

    {
        double m[1] = { mass_m1 };
        measure_loops_smeared(m, nhits, 1, precision, source_type, n_mom, n_smr,
                              alpha, HYP_WEIGHTS, DONTSTORE, NULL);
    }

    /* ---------------------------------------------------------------- *
     * Section 2: Standard estimator at m2                              *
     * ---------------------------------------------------------------- */
    lprintf("MAIN", 0, "------------------------------------------------------------\n");
    lprintf("MAIN", 0, "  STANDARD ESTIMATOR at m2 = %f\n", mass_m2);
    lprintf("MAIN", 0, "------------------------------------------------------------\n");

    {
        double m[1] = { mass_m2 };
        measure_loops_smeared(m, nhits, 1, precision, source_type, n_mom, n_smr,
                              alpha, HYP_WEIGHTS, DONTSTORE, NULL);
    }

    /* ---------------------------------------------------------------- *
     * Section 3: Split-even telescopic estimator with {m1, m2}        *
     * ---------------------------------------------------------------- */
    lprintf("MAIN", 0, "------------------------------------------------------------\n");
    lprintf("MAIN", 0, "  SPLIT-EVEN ESTIMATOR with masses {%f, %f}\n", mass_m1, mass_m2);
    lprintf("MAIN", 0, "------------------------------------------------------------\n");
    {
        double masses[2] = { mass_m1, mass_m2 };
        measure_loops_spliteven(masses, 2, nhits, precision,
                                n_smr, alpha, HYP_WEIGHTS,
                                DONTSTORE, NULL);
    }

    /* ---------------------------------------------------------------- */
    double elapsed_sec = timer_lap(&clock) * 1.e-6;
    lprintf("TIMING", 0, "Test completed [%lf sec]\n", elapsed_sec);

    /* close communications */
    finalize_process();

    return 0;
}
