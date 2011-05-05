/************************************************************************
 *
 *  File name: interface.c
 *
 *  Description: main swig interface file
 *
 *  Copyright 2010 Jacobs University Bremen, all rights reserved
 *
 *************************************************************************/

#include "handle.h"

/******************************* Externs *********************************/

// global variables
extern int DEBUG_LEVEL;


//result table
extern Res_tbl RES_TBL;

extern char *input_network_fname;


/******************************************************************************/

int
make_network(Network** N_p, int* edges, int edges_num) {
    int i, j, rc, mat_val;
    unsigned s, t, max = 0;
    int self_edge_number = 0;
    rc = RC_OK;
    Network* N;
    N = (Network*) calloc(1, sizeof (Network));
    N->name = "net";
    N->edges_num = 0;
    N->edges_num = edges_num;
    N->e_arr = (Edge*) calloc(N->edges_num + 1, sizeof (Edge));
    for (i = N->edges_num; i > 0; --i) {
        j = i << 1;
        s = edges[j - 2];
        t = edges[j - 1];
        max = (s > t) ? ((s > max) ? s : max) : ((t > max) ? t : max);
        N->e_arr[i].s = s;
        N->e_arr[i].t = t;
        N->e_arr[i].weight = 1;
        if (s == t) {
            ++self_edge_number;
        }
    }
    if (self_edge_number > 0) {
        return RC_ERR;
    }
    N->vertices_num = max;
    rc |= MatInit(&N->mat, N->vertices_num, SPARSE);
    if (rc == RC_ERR) {
        return rc;
    }

    //assign matrix entries according to edges
    for (i = 1; i <= N->edges_num; i++) {
        if (MatGet(N->mat, N->e_arr[i].s, N->e_arr[i].t) != 0) {
            return RC_ERR;
        } else {
            MatAsgn(N->mat, N->e_arr[i].s, N->e_arr[i].t, 1);
        }
    }
    //allocate and fill arrays of single edges and double edges
    N->e_arr_sin = (Edge*) calloc(N->edges_num + 1, sizeof (Edge));
    N->e_arr_dbl = (Edge*) calloc((N->edges_num + 1), sizeof (Edge));
    N->e_sin_num = 0;
    N->e_dbl_num = 0;
    N->roots_num = 0;
    N->leafs_num = 0;
    //allocate indeg and out deg arrays
    N->indeg = (int*) calloc(N->vertices_num + 2, sizeof (int));
    N->outdeg = (int*) calloc(N->vertices_num + 2, sizeof (int));
    N->doubledeg = (int*) calloc(N->vertices_num + 2, sizeof (int));
    if (GNRL_ST.calc_self_edges == TRUE) {
        N->self_edge = (int*) calloc(N->vertices_num + 2, sizeof (int));
    }
    //actually matrix is sparse anyway now
    if (N->mat->type == SPARSE) {
        for (i = 1; i <= N->vertices_num; i++) {
            for (j = 1; j <= N->vertices_num; j++) {
                if ((mat_val = MatGet(N->mat, i, j))) {
                    //if an edge and is not self edge
                    if (i != j) {
                        //if the twin edge exists
                        if (MatGet(N->mat, j, i)) {
                            //not inserted yet
                            if (j > i) {
                                //double edge- this way always the twin pair has indexes 2x-1,2x
                                N->e_arr_dbl[++N->e_dbl_num].s = i;
                                N->e_arr_dbl[N->e_dbl_num].t = j;
                                N->e_arr_dbl[N->e_dbl_num].weight = mat_val;
                                N->e_arr_dbl[++N->e_dbl_num].s = j;
                                N->e_arr_dbl[N->e_dbl_num].t = i;
                                N->e_arr_dbl[N->e_dbl_num].weight = mat_val;
                            }
                        } else {
                            //single edge
                            N->e_arr_sin[++N->e_sin_num].s = i;
                            N->e_arr_sin[N->e_sin_num].t = j;
                            N->e_arr_sin[N->e_sin_num].weight = mat_val;
                        }
                    } else //self edge
                    {
                        N->e_arr_sin[++N->e_sin_num].s = i;
                        N->e_arr_sin[N->e_sin_num].t = j;
                        N->e_arr_sin[N->e_sin_num].weight = mat_val;
                    }
                }
            }
        }

        //fill in deg and out deg arrays
        for (i = 0; i <= N->vertices_num; i++) {
            N->indeg[i] = 0;
            N->outdeg[i] = 0;
            if (GNRL_ST.calc_self_edges == TRUE) {
                N->self_edge[i] = FALSE;
            }
        }
        for (i = 1; i <= N->vertices_num; i++) {
            if (N->mat->spr->m[i].to == NULL)
                N->outdeg[i] = 0;
            else
                N->outdeg[i] = N->mat->spr->m[i].to->size;
            if (N->mat->spr->m[i].from == NULL)
                N->indeg[i] = 0;
            else
                N->indeg[i] = N->mat->spr->m[i].from->size;
            if ((N->mat->spr->m[i].self_edge == 1) && (GNRL_ST.calc_self_edges == TRUE)) {
                N->self_edge[i] = TRUE;
            }
        }
    }

    //statistics and global info about the network
    N->con_vertices_num = 0;
    N->hub_deg = 0;
    N->hub = 0;
    N->in_hub_deg = 0;
    N->in_hub = 0;
    N->out_hub_deg = 0;
    N->out_hub = 0;
    //calc total num of connected vertices
    //and find hub preferenced
    for (i = 1; i <= N->vertices_num; i++) {
        if ((N->indeg[i] != 0) || (N->outdeg[i] != 0))
            N->con_vertices_num++;
        if (((N->indeg[i] + N->outdeg[i]) > N->hub_deg)) {
            N->hub_deg = N->indeg[i] + N->outdeg[i];
            N->hub = i;
        }
        if (N->indeg[i] > N->in_hub_deg) {
            N->in_hub_deg = N->indeg[i];
            N->in_hub = i;
        }
        if (N->outdeg[i] > N->out_hub_deg) {
            N->out_hub_deg = N->outdeg[i];
            N->out_hub = i;
        }
    }
    *N_p = N;
    return RC_OK;
}

Res_tbl*
subgraphs_interface(int* edges, int edges_num, int mtf_sz, int max) {
    int rc = 0;
    Network *N = NULL;

    input_network_fname = "dummy";

    //init GNRL_ST accoridng to default then overide
    //according to arguments
    GNRL_ST.mtf_sz = mtf_sz;
    GNRL_ST.rnd_net_num = 0;
    GNRL_ST.t_init = T0;
    GNRL_ST.iteration_factor = ITERATION_F;
    GNRL_ST.e_thresh = ETHRESH;
    GNRL_ST.use_stubs_method = FALSE;
    GNRL_ST.long_out_flag = FALSE;
    GNRL_ST.calc_unique_flag = TRUE;
    GNRL_ST.calc_roles = FALSE;
    GNRL_ST.input_net_format = SRC_TRG_FORMAT;
    GNRL_ST.undirected_flag = FALSE;
    GNRL_ST.calc_self_edges = FALSE;
    GNRL_ST.calc_weights = FALSE;
    GNRL_ST.run_prob_app = FALSE;
    GNRL_ST.prob_base_samples_num = 0;
    GNRL_ST.prob_converge_mode = FALSE;
    GNRL_ST.prob_conv_diff = CONVERGNESS_DIFF_CONST;
    GNRL_ST.prob_conv_conc_th = CONC_THRESHOLD;
    GNRL_ST.unique_th = UNIQUE_TH;
    GNRL_ST.force_unique_th = FALSE;
    GNRL_ST.mfactor_th = MFACTOR_TH;
    GNRL_ST.zfactor_th = ZSCORE_TH;
    GNRL_ST.pval_th = PVAL_TH;
    GNRL_ST.max_members_list_sz = max;
    GNRL_ST.top_motifs_list_sz = TOP_MOTIF_LIST_SZ;
    GNRL_ST.out = FALSE;
    GNRL_ST.out_log_flag = FALSE;
    GNRL_ST.out_s_mat_flag = FALSE;
    GNRL_ST.out_l_mat_flag = FALSE;
    GNRL_ST.out_s_c_mat_flag = FALSE;
    GNRL_ST.out_members = FALSE;
    GNRL_ST.out_rand_mat_flag = FALSE;
    GNRL_ST.out_roles = FALSE;
    GNRL_ST.out_clustering = FALSE;
    GNRL_ST.quiet_mode = TRUE;
    GNRL_ST.use_metropolis = FALSE;
    GNRL_ST.use_clustering = FALSE;
    GNRL_ST.dont_search_real = FALSE;
    GNRL_ST.out_intermediate = FALSE;
    GNRL_ST.out_non_dangling_motifs = FALSE;
    GNRL_ST.r_switch_factor = 0.;
    GNRL_ST.r_global_switch_mode = FALSE;
    GNRL_ST.list_members = TRUE;
    GNRL_ST.specific_subgraph_members = 0; //no specific subgraph members
    GNRL_ST.efc_prob = FALSE;
    GNRL_ST.out_all_rand_networks = FALSE;
    GNRL_ST.r_grassberger = FALSE;
    GNRL_ST.r_grass_colony_sz = 0;
    GNRL_ST.r_grass_colony_max_population = MAX_COLONY_SZ_RATIO;
    GNRL_ST.r_dont_conserve_mutuals = FALSE;
    GNRL_ST.dont_die = FALSE;
    GNRL_ST.r_conserve_layers = FALSE;
    GNRL_ST.r_layers_num = 0;

    // initialisation of mfinder specific structs
    init_res_tbl(&RES_TBL);

    rc = gnrl_init();
    if (rc == RC_ERR) {
        printf("general init failed\n");
        return NULL;
    }

    rc = make_network(&N, edges, edges_num);
    if (rc == RC_ERR) {
        printf("network init failed\n");
        return NULL;
    }

    //search motif size n
    count_subgraphs_size_n(N, GNRL_ST.mtf_sz, &RES_TBL.real, REAL_NET, 0);

    //calc result after isomorphism of ids
    join_subgraphs_res(&RES_TBL.real, GNRL_ST.mtf_sz, 0);

    free_network_mem(N);
    free(N);

    return &RES_TBL;
}