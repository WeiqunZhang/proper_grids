
{
    RealBox rb({AMREX_D_DECL(0., ymin, 0.)}, {AMREX_D_DECL(xmax, ymax, 1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());

    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cellx - 1, n_celly - 1, 1)});

    // defining geometries
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain,&rb, 0, is_periodic.data());
        domain.refine(ref_ratio);
    }

    // defining grid at level 0
    grids[0].define(domain);
    grids[0].maxSize(max_grid_size);

    // defining grids for the remaining levels
    for (int ilev = 1; ilev < nlevels; ++ilev) {        

        // We place the blocks/Boxes of level ilev in a BoxList dlist
        BoxList dlist;

        // for each level, we loop over the blocks and build their grids
        for (int ii = 0; ii < Nb; ++ii) {

            domain = domain0;

            if (bll[ii]>=ilev) {          
                
                // sizetoN is a function to transform the physical
                // location/size into number of points

                // resize and shift the block/Box to the given location and size in the x-direction
                Nsize = sizetoN(blg[4*ii + 1]);
                domain.grow(0, -n_cellx/2 + Nsize/2);
                              
                Nshift=sizetoN(blg[4*ii]);
                domain.shift(0, -n_cellx/2 + Nsize/2 + Nshift);

                // resize and shift the block/Box to the given location and size in the y-direction
                Nsize = sizetoN(blg[4*ii + 3]); 
                domain.grow(1, -n_celly/2 + Nsize/2);
                                
                Nshift=sizetoN(blg[4*ii + 2]); 
                domain.shift(1, -n_celly_min - n_celly/2 + Nsize/2 + Nshift);

                // refine the mesh in the block
                for (int jj = 1; jj <= ilev; ++jj)
                    domain.refine(ref_ratio);

                // fill the fill the block/Box in dlist
                dlist.push_back(domain);
            }
        }
        
        // Define the grid in terms of the BoxList dlist
        grids[ilev].define(dlist);
        grids[ilev].maxSize(max_grid_size);

    }
}

