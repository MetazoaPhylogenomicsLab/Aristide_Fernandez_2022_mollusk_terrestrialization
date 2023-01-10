          seed = -1
       seqfile = matrix50genes_10partitions_rates.phy
      treefile = species_tree_calibrations.nwk
       outfile = res_usedata2.txt

         ndata = 10
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<7'  * safe constraint on root age, used if no fossil for root.
       

         model = 0   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

         alpha = 1    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0   * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 2 2      * gamma prior for alpha

   rgene_gamma = 1 10    * gammaDir prior for rate for genes
  sigma2_gamma = 1 10   * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 60000
      sampfreq = 60
       nsample = 10000
