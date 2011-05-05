/************************************************************************
 *
 *  File name: main.c
 *
 *  Description: main file
 *
 *  Copyright � 2002-2004 Weizmann Institute of Science,
 *			   76100 Rehovot Israel, All rights reserved
 *
 *************************************************************************/

#include "handle.h"

/******************************* Externs *********************************/

// global variables
extern int DEBUG_LEVEL;


extern Network *G_N;
extern Network *G_WN;
extern Gnrl_st GNRL_ST;
//result table
extern Res_tbl RES_TBL;

extern char *input_network_fname;

//Gneralized motifs res table
extern Res_tbl GMTF_RES_TBL;
extern Res_tbl GMTF_CMPLX_RES_TBL;

extern FILE *out_fp;
extern FILE *log_fp;
extern time_t start_time;
extern time_t end_time;

extern list64 *final_res;
extern list64 *final_res_all;

extern list64 *res_sub_motif;


//for grassberger
extern int died_colonies;
extern int over_populated_colonies;
extern double *weights_arr;

extern Hash *Role_hash;
extern role_st *Role_trans;
//Roles result table
extern Res_vec_tbl Roles_res_tbl;
//final results for role statistics
extern Role_res *Roles_final_res;
//Roles members 2 dim array [Role_id][node]
list **Roles_members;

//Gneralized motifs final res lists
extern list64 *gmtf_final_res;
extern list64 *gmtf_cmplx_final_res;

extern void
init_random_seed();

// return random integer in the interval 1 to max_val
extern int
get_rand(int max_val);
extern double get_rand_double();

/******************************************************************************/

int
main(int argc, char *argv[]) {
    int rc = 0;
    Network *N = NULL;

    time_measure_start(&GNRL_ST.total_time);

    //process input arguments and init global structure
    rc |= process_input_args(argc, argv);

    if (GNRL_ST.quiet_mode == FALSE)
        printf("mfinder Version %.2f\n\n", VERSION);


    if (rc == RC_ERR)
        at_exit(-1);

    //general initialization
    rc |= gnrl_init();
    if (rc == RC_ERR)
        at_exit(-1);

    // load network from input file
    if (GNRL_ST.quiet_mode == FALSE)
        printf("Loading Network\n");
    load_network(&G_N, input_network_fname);
    duplicate_network(G_N, &N, "real_network");

    init_random_seed();

    if (rc == RC_ERR)
        at_exit(-1);
    if (GNRL_ST.quiet_mode == FALSE)
        printf("Searching motifs size %d\nProcessing Real network...\n", GNRL_ST.mtf_sz);

    //search motifs size n in Real network
    if (GNRL_ST.dont_search_real != TRUE)
        rc |= motifs_search_real(N);
    if (rc == RC_ERR)
        at_exit(-1);

/*
    printf("RES TBL real\n");
    for (l_id = list64_get_next(RES_TBL.real, NULL); l_id != NULL;
            l_id = list64_get_next(RES_TBL.real, l_id))
    {
        printf("id: %d\ncount: %.0f\n", (int)((Motif*)l_id->p)->id,
                ((Motif*)l_id->p)->count);
    }
*/


    if (GNRL_ST.quiet_mode == FALSE)
        printf("Processing Random networks\n");
    if (GNRL_ST.rnd_net_num > 0) {
        // create random networks with same single node statisticfs as the input network
        if (GNRL_ST.r_grassberger == FALSE) {
            //use switches or stubs
            rc |= process_rand_networks(&RES_TBL, GNRL_ST.mtf_sz);
        } else {
            //use grassberger alg
            weights_arr = (double*) calloc(GNRL_ST.rnd_net_num + 1, sizeof (double));
            rc |= process_rand_networks_grassberger(&RES_TBL, GNRL_ST.mtf_sz, weights_arr);
        }
        if (rc == RC_ERR)
            at_exit(-1);
    }

    if (GNRL_ST.quiet_mode == FALSE)
        printf("Calculating Results...\n");
    if (GNRL_ST.rnd_net_num >= 0) {
        //calculate final results and dump them to the results file
        if (!GNRL_ST.r_grassberger) {
            calc_final_results(&RES_TBL, &final_res, &final_res_all, GNRL_ST.rnd_net_num);
        } else {
            //Nadav change for GRASS NEW
            calc_final_results_grassberger(&RES_TBL, FALSE, res_sub_motif, &final_res, &final_res_all, GNRL_ST.rnd_net_num, weights_arr);
        }
    }

    //calculate final results
    time_measure_stop(&GNRL_ST.total_time);
    //output results
    rc |= output_results(final_res, final_res_all);

    free_network_mem(G_N);
    free(G_N);

    final_res_free(final_res);
    final_res_free(final_res_all);

    if (GNRL_ST.r_grassberger)
        free(weights_arr);

    res_tbl_mem_free(&RES_TBL);

    if (GNRL_ST.calc_roles == TRUE) {
        free_mem_role_hash();
        free_roles_res_tbl(GNRL_ST.rnd_net_num);
        free_role_members();
    }

    if (rc == RC_ERR)
        at_exit(-1);

    exit(at_exit(rc));
}
