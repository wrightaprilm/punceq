
taxa <- readTaxonData("data/Simulated_Fossil_Intervals_Thu_Sep_08_11-22-27_2022.tsv")
morpho_f <- readDiscreteCharacterData("data/Simulated_Character_Data_Thu_Sep_08_11-22-27_2022.nex")

moves = VectorMoves()
monitors = VectorMonitors()

n_taxa <- taxa.size()
num_branches <- 2 * n_taxa - 2

# Diversification Rates based on Echinodermata
      speciation_rate ~ dnExponential(1.0);
      # NOTE: If it gets stuck in this script, then set origination & extinction to 1.0
      moves.append(mvScale(speciation_rate, lambda=0.01, weight=5));
      moves.append(mvScale(speciation_rate, lambda=0.10, weight=3));
      moves.append(mvScale(speciation_rate, lambda=1.00, weight=1));

#      turnover ~ dnUnif(0.9, 1.05);
#      moves.append(mvSlide(turnover, delta=0.01, weight=5));
#      moves.append(mvSlide(turnover, delta=0.10, weight=3));
#      moves.append(mvSlide(turnover, delta=1.00, weight=1));
#      extinction_rate := turnover*speciation_rate;
        extinction_rate ~  dnExponential(1.0);

      diversification := speciation_rate - extinction_rate;

      # Fossil Sampling Rates based on collection occupied by Echinodermata
      psi ~ dnExponential(1.0);
      completeness := psi/(extinction_rate+psi);
      moves.append(mvScale(psi, lambda=0.01, weight=5));
      moves.append(mvScale(psi, lambda=0.10, weight=3));
      moves.append(mvScale(psi, lambda=1.00, weight=1));

      # Proportional Taxon Sampling of Youngest Time Slice
      rho <- 0.43;	# 'extant' sampling.

     # Establish Basal Divergence Time
     origin_time ~ dnUnif(5, 8);
     moves.append(mvSlide(origin_time, delta=0.01, weight=5));
     moves.append(mvSlide(origin_time, delta=0.10, weight=3));
     moves.append(mvSlide(origin_time, delta=1.00, weight=1));


     fbd_dist = dnBDSTP(originAge=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, taxa=taxa, condition="sampling")


     outgroup = clade("A0", "P0");

     constraints = v(outgroup)


     fbd_tree ~ dnConstrainedTopology(fbd_dist, constraints=constraints)
     fbd_tree

     moves.append(mvFNPR(fbd_tree, weight=15.0))
     moves.append(mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6.0))
     moves.append(mvNodeTimeSlideUniform(fbd_tree, weight=40.0))
     moves.append(mvRootTimeSlideUniform(fbd_tree, origin_time, weight=5.0))


 # Setup the fossil tip sampling #

 # Use a for loop to create a uniform distribution on the occurence time for each fossil #

 # The boundaries of the uniform distribution are specified in the tsv file #

 fossils = fbd_tree.getFossils()
 for(i in 1:fossils.size())
 {
     t[i] := tmrca(fbd_tree, clade(fossils[i]))
#
     a_i = fossils[i].getMinAge()
     b_i = fossils[i].getMaxAge()
#
     F[i] ~ dnUniform(b_i, a_i)
     F[i].clamp( 0 )
 }

 # Add a move to sample the fossil times #
 #moves.append( mvFossilTimeSlideUniform(fbd_tree, origin_time, weight=5.0) )


   num_samp_anc := fbd_tree.numSampledAncestors()

#n_branches <- 2 * n_taxa - 2
#ucln_mean ~ dnExponential(2.0)
#ucln_sigma ~ dnExponential(3.0)
#ucln_var := ucln_sigma * ucln_sigma
#ucln_mu := ln(ucln_mean) - (ucln_var * 0.5)
#moves.append( mvScale(ucln_mean, lambda=1.0, tune=true, weight=4.0))
#moves.append( mvScale(ucln_sigma, lambda=0.5, tune=true, weight=4.0))
#for(i in 1:n_branches){
#   branch_rates[i] ~ dnLnorm(ucln_mu, ucln_sigma)
#   moves.append( mvScale(branch_rates[i], lambda=1, tune=true, weight=2.))
#}

clock_morpho ~ dnExponential(1.0)

moves.append( mvScale(clock_morpho, weight=4.0) )


     alpha_morpho ~ dnUniform( 0, 1E6 )

     rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )

     #Moves on the parameters of the Gamma distribution.

     moves.append(mvScale(alpha_morpho, lambda=1, weight=2.0))

     rate_pr := fbd_tree.treeLength() / 10
     num_rates = 3 * (3-1)
     for ( i in 1:num_rates ) {
         rate[i] ~ dnLognormal(1, 1)
         scaled_rate[i] := rate[i]*speciation_rate
         moves.append( mvScale( rate[i], weight=2 ) )
     }


     ##########################
     # Set up the rate matrix #
     ##########################

  #   Q_morpho := fnFreeK( scaled_rate, rescale=false )
      Q_morpho := fnJC(3)
     # create model of evolution for the character block
      m_morph ~ dnPhyloCTMC( tree=fbd_tree,
                                         Q=Q_morpho,
                                         siteRates=rates_morpho,
                                         branchRates=clock_morpho,
                                         type="Standard")
     	    m_morph.clamp(morpho_f)




     mymodel = model(fbd_tree)



     monitors.append(mnModel(filename="output/cinc6_dated.log", printgen=10))



     monitors.append(mnFile(filename="output/cinc6_dated.trees", printgen=10, fbd_tree))



     monitors.append(mnScreen(printgen=10, num_samp_anc, origin_time))



     mymcmc = mcmc(mymodel, monitors, moves)


     mymcmc.run(generations=1000000)



     q()
