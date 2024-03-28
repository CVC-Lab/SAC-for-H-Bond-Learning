#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include _MOL_INCLUDE_


// Creates an atomgroup from user specified paramters (See header file for documentation)
struct atomgrp *create_atomgrp (char *pdbFile, char *pdbFixedFile, char *psfFile,
                                char *mol2File, char *prmFile, char *rtfFile,
                                char *aprmFile, int useHbond)
{

    // Return null if we are finding hbond energy (THIS COULD CHANGE IN THE
    // FUTURE IF WE ADD MORE ENERGY TYPES)
    if (!useHbond){
        return NULL;
    }

    // Read atomic parameter file int new prm struct
    printf( "Processing %s... ", aprmFile ); fflush( stdout );
    struct prm *aprm = read_prm( aprmFile, "0.0.6" );   // read in atomic parameters
    printf( "done\n" );
    // Read pdb file int a new atomgroup struct
    printf( "Processing %s... ", pdbFile ); fflush( stdout );
    struct atomgrp *ag = read_pdb( pdbFile, aprm );     // read in input file and it's forcefield parameters
    printf( "done\n" );

    double start, end, charmmTime;
    // Read charm psf prm and rtf files into the atomgroup struct
    printf( "Processing %s, %s and %s... ", psfFile, prmFile, rtfFile ); fflush( stdout );
    start = clock( );
    read_ff_charmm( psfFile, prmFile, rtfFile, ag );    // read in CHARMM forcefield parameters
    end = clock( );   
    charmmTime = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
    printf( "%.3lf sec\n", charmmTime );

    // Check if we are going to keep track on hydrogen bonds
    if ( useHbond )
    {   
        // Read in the hybridization states from the mol2 file
        printf( "Processing %s... ", mol2File ); fflush( stdout );
        if ( !read_hybridization_states_from_mol2( mol2File, pdbFile, ag ) ) return 1;  // read in hybridization states
        printf( "done\n" );
        
        printf( "Fixing acceptor bases... " ); fflush( stdout );
        fix_acceptor_bases( ag, aprm );  // find base1 and base2 of each acceptor
        printf( "done\n" );
    }

    int nfix, *fix = ( int * ) mymalloc( ag->natoms * sizeof( int ) );
    // Read Fixed pdb file
    if ( pdbFixedFile != NULL )
    {
        printf( "Processing %s... ", pdbFixedFile ); fflush( stdout );
    }
    // FOR NOW pdbFixFile IS ALWAYS NULL
    read_fix( pdbFixedFile, &nfix, &fix );  // read in fixed part of the molecule
    fixed_init( ag );
    fixed_update( ag, nfix, fix );
    if ( pdbFixedFile != NULL ) printf( "done\n" );


    zero_grads( ag );  	// finish structure initialization
    fill_ingrp( ag );
    assign_combined_residue_sequence( ag );
    return ag;
}

// Runs hydrogen bonding minimization with the Octree
void hydro_min (struct atomgrp * ag, char* pdbFile, char *outnFile, char *outoFile, 
    int maxIter, int useHbond, double *dcHbond, double *dcVdw, double *dcElec) {
    printf ("Pointer of atomgroup: %lu\n", ag);
    // Storing initial active list of atoms before minimization in array 
    int ndim = ag->nactives * 3;
    double *startag = mymalloc( ndim * sizeof( double ) );
    ag2array( startag, ag );

    double nblist_cons = 0, octree_cons = 0, start = 0, end = 0;
    int snblist = 0, soctree = 0;

    // Initialize the nonbonded list and agsetup
    struct agsetup *ags = mymalloc( sizeof( struct agsetup ) );
    init_nblst( ag, ags );

    // Set nblist params
    ags->nblst->nbcof = 0;
    // Manually deactivate vdw and elec energies
    if ( 0 && ( *dcVdw > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcVdw;
    if ( 0 && ( *dcElec > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcElec;
    if ( useHbond && ( *dcHbond > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcHbond;

    ags->nblst->nbcut = ags->nblst->nbcof + 1;

    OCTREE octree;
    printf( "Building octree... " ); fflush( stdout );
    // Build the octree 
    start = clock( );
    if ( !build_octree( &octree, 60, 6, 1.0, ag ) )
        {
        print_error( "Failed to build octree!" );
        return 1;
        }
    end = clock( );
    octree_cons = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
    printf( "%.3lf sec\n", octree_cons );
    // Set octree size
    soctree = get_octree_size( &octree );

    printf("\n");

    // Set Octree params
    OCTREE_PARAMS eng_params;
    eng_params.ags = ags;
    eng_params.eps = 1.0;
    eng_params.octree_static = &octree;
    eng_params.octree_moving = &octree;
    eng_params.hdist_cutoff = 3.0;
    eng_params.fixed_cull = 1;
    eng_params.trans = NULL;
    eng_params.engcat = NULL;
    eng_params.proc_func_params = &eng_params;
    // Set minimization params
    struct my_par mupa;
    mupa.nbupdate = 1;
    mupa.ag = ag;
    mupa.ags = ags;
    mupa.efacv = 1.0;
    mupa.eface = 1.0;
    mupa.useVdw = 0;
    mupa.useElec = 0;
    mupa.useHbond = useHbond;
    mupa.dcVdw = 12.0;
    mupa.dcElec = 15.0;
    mupa.dcHbond = *dcHbond;
    mupa.aprxdcVdw = 1.0 * mupa.dcVdw;
    mupa.aprxdcElec = 1.0 * mupa.dcElec;

    double nblist_init_E = 0, octree_init_E = 0, nblist_E = 0, octree_E = 0, nblist_min = 0, octree_min = 0;
    const char *min_method[ ] = { "LBFGS" };
    mupa.octpar = &eng_params;

    // Run the minimization 
    my_en_grad( 0, NULL, ( void * )( &mupa ), &octree_init_E, NULL );
    printf( "\nApplying %s with OCTREES ( initial energy = %lf kcal/mol )...\n\n", min_method[ 0 ], octree_init_E ); fflush( stdout );
    start = clock( );
    minimize_ag( 0, maxIter, 1E-5, ag, ( void * )( &mupa ), my_en_grad );
    end = clock( );
    octree_min = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
    my_en_grad( 0, NULL, ( void * )( &mupa ), &octree_E, NULL );
    printf("\ndone ( time = %.3f sec, final energy = %f kcal/mol )\n\n", octree_min, octree_E );
    
    // Write to the out file
    printf( "Writing %s... ", outoFile ); fflush( stdout );
    write_pdb_nopar( ag, pdbFile, outoFile );
    printf( "done\n" );
    // Reset the atomic coordinates of the atomgroup struct
    array2ag( startag, ag );    
    zero_grads( ag );                       
    // Free the octree
    destroy_octree( &octree );

    // Print the results
    printf( "OCTREE results:\n" );     
    printf( "\tconstruction time = %.3lf sec\n", octree_cons );
    printf( "\trunning time = %.3lf sec\n", octree_min );                   
    printf( "\tsize = %d bytes ( %.2lf KB, %.2lf MB )\n",
            soctree, soctree / 1024.0, soctree / ( 1024.0 * 1024.0 ) );
    printf( "\tinitial energy = %lf kcal/mol\n", octree_init_E );                                         
    printf( "\tfinal energy = %lf kcal/mol\n", octree_E );                         
    printf( "\n" );                                    
    printf( "### END: RUN ###\n\n");         
    // Free the atom setup struct
    free_agsetup( ags );
}

// Creates a my_par struct from parameters and builds the octree
struct my_par *create_par(struct atomgrp *ag, int useVdw, int useElec, 
                         int useHbond ,double *dcVdw, double *dcElec, 
                         double *dcHbond, double aprxVdw, double aprxElec) {
    
    // Initialize the nonbonded list and agsetup
    struct agsetup *ags = mymalloc( sizeof( struct agsetup ) );
    init_nblst( ag, ags );

    // Set nblist params
    ags->nblst->nbcof = 0;
    // Manually deactivate vdw and elec energies
    if ( 0 && ( *dcVdw > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcVdw;
    if ( 0 && ( *dcElec > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcElec;
    if ( useHbond && ( *dcHbond > ags->nblst->nbcof ) ) ags->nblst->nbcof = *dcHbond;
    ags->nblst->nbcut = ags->nblst->nbcof + 1;
    
    // Build the octree 
    printf( "Building octree... " ); fflush( stdout );
    OCTREE *octree = malloc (sizeof (OCTREE));
    if ( !build_octree( octree, 60, 6, 1.25, ag ) )
        {
        print_error( "Failed to build octree!" );
        return NULL;
        }
        
    // Get octree size
    //int soctree = get_octree_size( octree );

    printf("\n");
    // Set Octree params
    OCTREE_PARAMS *eng_params = malloc (sizeof (struct OPAR));
    eng_params->ags = ags;
    eng_params->eps = 1.0;
    eng_params->octree_static = octree;
    eng_params->octree_moving = octree;
    eng_params->hdist_cutoff = 3.0;
    eng_params->fixed_cull = 1;
    eng_params->trans = NULL;
    eng_params->engcat = NULL;
    eng_params->proc_func_params = eng_params;

    // Set my_par parameters
    struct my_par *parm = malloc (sizeof(struct my_par));
    parm->nbupdate = 1;
    parm->ag = ag;
    parm->ags = ags;
    parm->efacv = 1.0;
    parm->eface = 1.0;
    parm->useVdw = useVdw;
    parm->useElec = useElec;
    parm->useHbond = useHbond;
    parm->dcVdw = *dcVdw;
    parm->dcElec = *dcElec;
    parm->dcHbond = *dcHbond;
    parm->aprxdcVdw = aprxVdw * parm->dcVdw;
    parm->aprxdcElec = aprxElec * parm->dcElec;
    parm->octpar = eng_params;

    return parm;
}

// Frees the memory of the my_par struct
void free_parm (struct my_par *parm) {
    struct agsetup *ags = parm->ags;
    OCTREE_PARAMS *op = parm->octpar;
    // Free octree
    destroy_octree( op->octree_moving);
    // Free parm and octree parm
    freeMem (parm);
    freeMem (op);
    // Free agsetup 
    free_agsetup (ags);
}

// Set all num_neighbors to 0
void clear_neighbors (struct atomgrp *ag) {
  int i;
  for (i = 0; i < ag->natoms; i++) {
    ag->atoms[i].num_neighbors = 0;
  }
  return;
}

// Finds the neighborhood of each atom
// Must be sent into the accumulate function
void find_neighbors (OCTREE_PARAMS *octpar, double *energy) {
    OCTREE *octree_static = octpar->octree_static;
    OCTREE *octree_moving = octpar->octree_moving;   

    
    double dist_cutoff = octpar->dist_cutoff;
    double *trans_mat = octpar->trans;

    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) octpar->proc_func_params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ octpar->node_static ] );
    OCTREE_NODE *mnode = &( octree_moving->nodes[ octpar->node_moving ] );

    double rc = dist_cutoff;
    double rc2 = rc * rc;

    int nf = snode->nfixed;

    double *engcat = prms->engcat;
    
    *energy = 0;

    // Case where the moving node has more atoms than the static node
    if ( ( trans_mat != NULL ) || ( mnode->n - mnode->nfixed <= snode->n ) )
      {
        // Loop through all atoms in the moving node
        for ( int i = mnode->nfixed; i < mnode->n; i++ )
          {
            // Grab coords of the atom in the moving node
            int ai = mnode->indices[ i ];
            mol_atom *atom_i = &( octree_moving->atoms[ ai ] );
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
            /*if (ai == 8) {
              for ( int k = 0; k < snode->n; k++ ) {
                printf ("Checking node with 8 compared to node that has: %d\n", snode->indices[ k]);
              }
            }*/
            
            if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h
      
            double d2 = min_pt2bx_dist2( snode->lx, snode->ly, snode->lz, snode->dim, x, y, z );
            // Check if moving node is within cutoff from "ai"
            
            // Check if the static node is within the cutoff from "ai"
            if ( d2 < rc2 )
              {
                // Check if each atom is within the cutoff
                for ( int j = 0; j < snode->n; j++ )
                    {
                      int aj = snode->indices[ j ];
                      //if (ai == 0) {
                      //   printf("saj: %d\n", aj);
                      //}
                      // if (ai == 2 && aj == 8) {
                      //   printf ("sNodes Addrs for ai==2: %p\n Nodes Addrs for aj==8: %p\n", snode, mnode);
                      // }
                      
                      // Skip the the atom if its index number is less than "ai"
                      if ( ( j >= nf ) && ( aj <= ai ) ) continue;
                      //if ( (aj == ai ) ) continue;
                                          
                      mol_atom *atom_j = &( octree_static->atoms[ aj ] );

                      // Calculate the dist between atoms          
                      double dx = atom_j->X - x, 
                      dy = atom_j->Y - y,
                      dz = atom_j->Z - z;
                      d2 = dx * dx + dy * dy + dz * dz;
                      
                      // Skip the atom if not within dist cutoff 
                      if ( d2 >= rc2 ) continue; 

                      // Uncomment these lines to exclude atoms in the exclusion list (i.e. bonded atoms, atoms that are too close)
                      //int k;
                      //if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                      //else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                      //if ( k > 0 ) continue;
                      add_neighbor (atom_i, aj);
                      add_neighbor (atom_j, ai);
                    }
              }
          }
      }
    // Case where the static node has more atoms than the moving node
    else
      {
        // Loop through all atoms the static node
        for ( int i = 0; i < snode->n; i++ )
          {
            // Grab coords of atom in static node
            int ai = snode->indices[ i ];
            mol_atom *atom_i = &( octree_static->atoms[ ai ] );
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
            /*if (ai == 8) {
              for ( int k = 0; k < mnode->n; k++ ) {
                printf ("Checking node with 8 compared to node that has: %d\n", mnode->indices[ k]);
              }
            }*/
      
            double d2 = min_pt2bx_dist2( mnode->lx, mnode->ly, mnode->lz, mnode->dim, x, y, z );
            // Check if moving node is within cutoff from "ai"
            if ( d2 < rc2 )
              {
                // Loop through all the atoms in the moving node
                for ( int j = mnode->nfixed; j < mnode->n; j++ )
                    {
                      int aj = mnode->indices[ j ];
                      
                      if ( ( i >= nf ) && ( aj >= ai ) ) continue;
                      //if ( (aj == ai ) ) continue;
                      //if (ai == 0) {
                      //  printf ("maj: %d\n", aj);
                      //}
                      mol_atom *atom_j = &( octree_moving->atoms[ aj ] );

                      // Check if the atoms are within the dist cutoff             
                      double dx = atom_j->X - x, 
                      dy = atom_j->Y - y,
                      dz = atom_j->Z - z;
        
                      d2 = dx * dx + dy * dy + dz * dz;
                      
                      if ( d2 >= rc2 ) continue; 

                      // Uncomment these lines to exclude atoms in the exclusion list (i.e. bonded atoms, atoms that are too close)
                      //int k;
                      //if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                      //else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                      //if ( k > 0 ) continue;
  
                      add_neighbor (atom_i, aj);
                      add_neighbor (atom_j, ai);
                    }
              }
          }      
      }
      return energy; 
}
// Grabs the bonds for a certain atom
int *get_bondls (struct atomgrp *ag, struct atom *a) {
    if (a->nbondis == 0) {
       printf("WARNING: There are no bonds for this atom");
       return 0;
    }
    else {
      // initalize bond list
      int *result = (int*) malloc (sizeof(int) * a->nbondis);
      mol_bond *bonds = ag->bonds;
      int i;
      // Loop through all bonds
      for (i = 0; i < a->nbondis; i++) {
        int id1 = bonds[a->bondis[i]].ai;
        int id2 = bonds[a->bondis[i]].aj;
        // Add the atom that is not the current one
        if (id1 == a->ingrp) {
          result[i] = id2;  
        }
        else {
          result[i] = id1;
        }
      }
      return result;
    }
}

// Adds the neighbor index to the atom's neighborhood list
void add_neighbor(struct atom *a, int neighbor_index) {
    // Resize if there are more neighbors than the capacity
    if (a->num_neighbors >= a->neighbors_cap)
    {
      a->neighbors_cap *= 2;
      a->neighbors = (int *) _mol_realloc (a->neighbors, sizeof (int) * a->neighbors_cap);
    }
    // Add the atom to the neighborhood list
    a->neighbors[a->num_neighbors] = neighbor_index;
    a->num_neighbors++;
}

// Updates the positions of each atom in the octree
// Calculates the gradients 
void my_en_grad( int n, double *inp, void *prms, double *en, double *grad)
{
    static int count = 0;
    struct my_par *prm = ( struct my_par * ) prms;
    struct atomgrp *mag = prm->ag;
    int update_neighbors = 1;

    /*printf("\n NEW INPUT \n");
    // Debugging: Print inp
    for (int index = 0; index < n; index++) {
      printf (" %f ", inp[index]);
    }*/

    if ( inp != NULL)
      {
        if ( n != mag->nactives * 3 )
          {
            print_error( "Mismatch in vector length ( my_en_grad )!");
            exit(0);
          }
        array2ag(inp, mag);
      }

    int mqcheck = ( ( struct my_par * ) prms )->nbupdate;
    struct agsetup *mags = ( ( struct my_par * ) prms )->ags;

    OCTREE_PARAMS *octpar = ( ( struct my_par * ) prms )->octpar;

    if ( octpar != NULL ) reorganize_octree( octpar->octree_static, 1 );
    else
      {
        if ( mqcheck ) check_clusterupdate( mag, mags );
      }

    // insert energy terms here
    *en = 0;
    zero_grads( mag );

    if ( octpar == NULL )
      {
        if ( prm->useVdw ) vdweng( mag,  en, mags->nblst );

        if ( prm->useElec ) eleng( mag, prm->eface, en, mags->nblst );

        if ( prm->useHbond )
          {
            double en1 = 0;
            hbondeng( mag, &en1, mags->nblst );
            *en += en1;
          }
      }
    else
      { 
        // Update the neighbors of each atom
        if (update_neighbors) {
            // Clear the neighbors list
            clear_neighbors (mag);
            // Find the neighbors of each atom based on the hbond distance cutoff and updates the neighbor lists of each atom
            octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcHbond, prm->dcHbond, octpar->fixed_cull, octpar->trans,
                                              octpar->proc_func_params, find_neighbors);
        }
        // Calculate vdw energies
        if ( prm->useVdw )
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcVdw, prm->aprxdcVdw, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, vdweng_octree_single_mol );
        // Calculate electrostatic energies
        if ( prm->useElec )
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcElec, prm->aprxdcElec, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, eleng_octree_single_mol );
        // Calcualte hbond energies
        if ( prm->useHbond ) {
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcHbond, prm->dcHbond, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, hbondeng_octree_single_mol );
        }
      }
    // Calculate bond, angle, and improper angle energy
    beng(mag,  en);
    aeng(mag,  en);
    ieng(mag,  en);
//  teng(mag,  en);

    count++;

    int i;
    if ( grad != NULL )
      {
        for ( i = 0; i < n / 3 ; i++ )
          {
            grad[ 3 * i ] = -1 * mag->atoms[ mag->activelist[ i ] ].GX;
            grad[ 3 * i + 1 ] = -1 * mag->atoms[ mag->activelist[ i ] ].GY;
            grad[ 3 * i + 2 ] = -1 * mag->atoms[ mag->activelist[ i ] ].GZ;
          }
      }
}

void only_find_neighbors(int n, double *inp, void *prms)
{
    struct my_par *prm = ( struct my_par * ) prms;
    struct atomgrp *mag = prm->ag;

    if ( inp != NULL)
      {
        if ( n != mag->nactives * 3 )
          {
            print_error( "Mismatch in vector length ( my_en_grad )!");
            exit(0);
          }
        array2ag(inp, mag);
      }

    int mqcheck = ( ( struct my_par* ) prms )->nbupdate;
    struct agsetup *mags = (( struct my_par* ) prms )->ags;

    OCTREE_PARAMS *octpar = ( ( struct my_par* ) prms )->octpar;

    if ( octpar != NULL ) reorganize_octree( octpar->octree_static, 1 );
    else
      {
        if ( mqcheck ) check_clusterupdate( mag, mags );
      }

    // Clear the neighbors list
    clear_neighbors (mag);
    // Find the neighbors of each atom based on the hbond distance cutoff and updates the neighbor lists of each atom
    octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcHbond, prm->dcHbond, octpar->fixed_cull, octpar->trans,
                                      octpar->proc_func_params, find_neighbors);
}

void read_fix( char *ffile, int *nfix, int **fix )
{
   int linesz = 91;
   char *buffer = mymalloc( sizeof( char ) * linesz );

   *nfix = 0;

   if ( ffile != NULL )
     {
       FILE* fp = fopen( ffile, "r" );

       if ( fp == NULL )
         {
           print_error( "Failed to open fixed PDB %s ( use the --nofixed option if the molecule has no fixed part )!", ffile );
           exit( 0 );
         }

       while ( fgets( buffer, linesz - 1, fp ) != NULL )
         {
           if ( !strncmp( buffer, "ATOM", 4 ) ) ( *nfix )++;
         }
       fclose( fp );
     }

   *fix = mymalloc( ( *nfix ) * sizeof( int ) );

   if ( ffile != NULL )
     {
       FILE *fp = fopen( ffile, "r" );
       int na = 0;
       while ( fgets( buffer, linesz - 1, fp ) != NULL )
         {
           if ( !strncmp( buffer, "ATOM", 4 ) ) ( *fix )[ na++ ] = atoi( buffer + 4 ) - 1;
         }
       fclose( fp );
     }

   free( buffer );
}

