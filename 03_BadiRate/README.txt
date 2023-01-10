
# BADIRATE ANALYSIS

# FILES:
	
 # OG_size_counts.txt -> Orthogroup/Family size data table from OrthoFinder analysis.
 # OG_tables_badirate -> Orthogroup size data (same as above but one file per OG) for running BadiRate.
 # Run_badirate.txt -> Commands used to run the Local Rates (LR) and Global Rates (GR) BadiRate models.
 # post_times_tree.nwk -> Dated phylogeny for the analysis.
 # Results_LRT.txt -> Estimated likelihoods for each model and OG, the likelihood ratio tests (LRT) and FDR adjusted p values. Note that a few OGs failed to run.
 # pars_significantFDR_LR.txt -> Parameter estimates for the OGs that had the LR model as a better fit (G0, L0, D0 are the LR model parameters).

# "Reconstructions" folder:
	# Reconstructions obtained from BadiRate output.
	# Note that BadiRate internal labeling of nodes is different from the labeling in APE's "phylo" object.
	# See "Mapping_nodelabels_Badirate_Phylo.txt" for correspondence.
	
	# All_OGs_GainLoss_dynamics_per_branch.txt -> Summary of reconstructed gains/losses across all branches.
	# All_OGs_netGains_dynamics_per_branch.txt -> Summary of net Gains/Losses per branch.
	# All_OGs_PresenceAbsence_at_nodes.txt -> Whether an OG existed at that branch or not.
	# All_OGs_size_at_nodes.txt -> Reconstructed OG size at given branch.
	# Gains/Loss_table_all.txt -> Minimum number of reconstructed gains/losses per OG per branch.

	# "Nodes" folder: (follows phylo object IDs)

		# OGs_NodeXXX.txt -> OGs present at that node (i.e. reconstructed gene repertoire or pseudogenome).

