/**** This is just the header output...got tired of looking at it ****/


if (verbose > 0)
  {
  printf("\n\n----------------------------------------\n\n\n");
  printf("LINEAGE RATE/TIME ANALYSIS FOR TREE %s\n\n",TreeName);
  switch (method)
	{
	case PENLIKET:printf("Method = Penalized Likelihood(T)\n");printf("Smoothing factor = %f\n",rbp->smoothing);break;
	case PENLIKE:printf("Method = Penalized Likelihood\n");printf("Smoothing factor = %f\n",rbp->smoothing);
		switch (rbp->NeighborPenalty)
			{
			case 0:printf("Penalty function = Ancestor-Descendant\n");break;
			case 1:printf("Penalty function = Neighbor variance\n");break;
			}
		switch (rbp->PenaltyType)
			{
			case 0:printf("Scale for rate penalty = ADDITIVE\n");break;
			case 1:printf("Scale for rate penalty = LOG\n");break;
			}
		printf("Minimum allowed rate = %f of initial average rate estimate\n",rbp->minRateFactor);
		printf("Minimum allowed duration on 0-length terminal branches = %f of root's age\n",rbp->minDurFactor);
		break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case GAMMA:printf("Method = Gamma-Negative-Binomial\n");break;
	case HMM:printf("Method = Hidden Markov\n");break;
	case NP:
		{
		printf("Method = Non-parametric (exp=%f)\n",rbp->npexp);
		switch (rbp->PenaltyType)
			{
			case 0:printf("Scale for rate penalty = ADDITIVE\n");break;
			case 1:printf("Scale for rate penalty = LOG\n");break;
			}
		if (gVarMinFlag)
		    printf("(Minimizing variance of local rates)\n");
		else
		    printf("(Minimizing local transformations (NPRS))\n");
		}
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");break;
	case TN:     printf("Optimization via Truncated-Newton (TN) method with bound constraints\n");break;
	}

  printf("\n----------------------------------------\n");

    printf("Substitution Model\n");
    if (rbp->RatesAreGamma)
    	{
	printf("\tRates are gamma distributed across sites\n");
	printf("\tShape parameter of gamma distribution = %6.2f\n",rbp->alpha);
	}
    else
	printf("\tRates are equal across sites\n");
	

    printf("Global/Local Search Parameters\n");
    printf("\tNumber of searches from random starts = %i\n",NUM_TIME_GUESSES);
    printf("\tNumber of restarts after each search = %i\n",rbp->num_restarts);
    printf("\tLocal perturbation on restarts = %4.3g\n", rbp->perturb_factor);
    printf("\tLocal fractional tolerance after restarts = %4.3g\n", rbp->local_factor);
	

	if (algorithm==TN)
		printf("Optimization parameters set automatically by TN routine\n");
	else
		{
	  	printf("Optimization parameters\n");
		printf("\tMaximum number of iterations  = %i\n", rbp->maxIter);
		printf("\tFunction Tolerance = %4.3g\n", rbp->ftol);
		printf("\tlinminOffset = %4.3g\n",rbp->linminOffset);
		printf("\tContract Factor = %4.3g\n", rbp->contractFactor);
		printf("\tMax number of contract iterations  = %i\n",rbp->maxContractIter);
    	}
  if (method==NP)
	{
  	if (gClampRoot)
		printf("Root rates are CLAMPED\n");
  	else
        	printf("Root rates are FREE\n");
	}
  if (gisConstrained && algorithm != TN)
    {
    printf("Time constraints are enforced using barrier optimization\n");
    printf("\tBarrier Tol = %4.3g\n", rbp->barrierTol);
    printf("\tMax Barrier Iter = %6i\n", rbp->maxBarrierIter);
    printf("\tInit Barrier Factor = %4.3g\n", rbp->initBarrierFactor);
    printf("\tBarrier Multiplier = %4.3g\n", rbp->barrierMultiplier);
    }
  printf("Length of tree input =  %g\n",treeLength(root));
  printf("Number of taxa  =  %i\n",numdesc(root));
  printf("Number of sites in sequences =  %li\n",rbp->numSites);
  printf("\n----------------------------------------\n");
  }
