/*******************************************************************************
 * Compute disconnected loops with split-even telescopic estimator
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif
#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

typedef struct input_loops {
    char masses_string[256];  /* semicolon-separated, e.g. "-0.5;-0.3;-0.1" */
    double precision;
    int nhits;
    double csw;
    int n_smr;
    double alpha;
    char configlist[256];
    input_record_t read[8];
} input_loops;

#define init_input_loops(varname) {                                                              \
    .read = {                                                                                    \
        { "Masses", "disc:masses = %s", STRING_T, (varname).masses_string },                     \
        { "inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision },        \
        { "number of inversions per cnfg", "disc:nhits = %d", INT_T, &(varname).nhits },         \
        { "csw", "disc:csw = %lf", DOUBLE_T, &(varname).csw },                                  \
        { "n_smr", "disc:n_smr = %d", INT_T, &(varname).n_smr },                                \
        { "alpha", "disc:alpha = %lf", DOUBLE_T, &(varname).alpha },                             \
        { "Configuration list:", "disc:configlist = %s", STRING_T, (varname).configlist },       \
        { NULL, NULL, INT_T, NULL }                                                              \
    }                                                                                            \
}

typedef struct input_HYP {
    double weight[3];
    input_record_t read[4];
} input_HYP;

#define init_input_HYP(varname) {                                                                \
    .read = {                                                                                    \
        { "HYP smearing weight[0]", "HYP:weight0 = %lf", DOUBLE_T, &((varname).weight[0]) },    \
        { "HYP smearing weight[1]", "HYP:weight1 = %lf", DOUBLE_T, &((varname).weight[1]) },    \
        { "HYP smearing weight[2]", "HYP:weight2 = %lf", DOUBLE_T, &((varname).weight[2]) },    \
        { NULL, NULL, INT_T, NULL }                                                              \
    }                                                                                            \
}

char input_filename[256] = "input_file_spliteven";

input_loops disc_var = init_input_loops(disc_var);
input_HYP   HYP_var  = init_input_HYP(HYP_var);

int main(int argc, char *argv[]) {
    FILE *list;
    char list_filename[256] = "";
    char cnfg_filename[256] = "";
    int i;

    Timer clock;
    timer_set(&clock);

    setup_process(&argc, &argv);
    setup_gauge_fields();

    read_input(disc_var.read, get_input_filename());
    HYP_var.weight[0] = HYP_var.weight[1] = HYP_var.weight[2] = 0.;
    read_input(HYP_var.read, get_input_filename());
    double *hyp = (HYP_var.weight[0] == 0. && HYP_var.weight[1] == 0. && HYP_var.weight[2] == 0.)
              ? NULL : HYP_var.weight;

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&disc_var.csw);
#endif

    /* parse semicolon-separated mass string */
    double masses[64];
    int n_masses = 0;
    char masses_copy[256];
    strncpy(masses_copy, disc_var.masses_string, 255);
    char *token = strtok(masses_copy, ";");
    while (token != NULL && n_masses < 64) {
        masses[n_masses++] = atof(token);
        token = strtok(NULL, ";");
    }
    error(n_masses < 2, 1, "main [compute_loops_spliteven.c]",
          "Need at least 2 masses for telescopic split-even");

    lprintf("MAIN", 0, "Number of masses: %d\n", n_masses);
    for (int k = 0; k < n_masses; k++)
        lprintf("MAIN", 0, "  masses[%d] = %f\n", k, masses[k]);
    lprintf("MAIN", 0, "nhits = %d\n", disc_var.nhits);
    lprintf("MAIN", 0, "n_smr = %d, alpha = %f\n", disc_var.n_smr, disc_var.alpha);
    lprintf("MAIN", 0, "HYP weights: %f %f %f\n",
            HYP_var.weight[0], HYP_var.weight[1], HYP_var.weight[2]);

    strcpy(list_filename, disc_var.configlist);
    error((list = fopen(list_filename, "r")) == NULL, 1,
          "main [compute_loops_spliteven.c]", "Failed to open config list\n");

    init_BCs(NULL);

    i = 0;
    while (++i) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) break;
        }

        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);
        represent_gauge_field();

        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
        full_plaquette();

        measure_loops_spliteven(masses, n_masses, disc_var.nhits, disc_var.precision,
                                disc_var.n_smr, disc_var.alpha, hyp,
                                DONTSTORE, NULL);

        if (list == NULL) break;
    }

    if (list != NULL) fclose(list);

    double elapsed = timer_lap(&clock) * 1.e-6;
    lprintf("TIMING", 0, "Done [%lf sec]\n", elapsed);

    finalize_process();
    return 0;
}