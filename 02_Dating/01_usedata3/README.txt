For multigene protein alignments, I followed MCMCtree tutorial 4

First run with usedata=3 (mcmctree_usedata3.ctl; results in res_1stpart)
Problem is, MCMCTree calls CODEML with model 0 (poisson), not with WAG+F, so out.BV (containing the hessian, branch lengths, etc for each partition)
is not good as it is.
To overcome this, you need after running with usedata=3 to rerun manually CODEML for each partition. We use the files generated (.ctl,.tree,.txt), but discard the out.bv.
We now modify each .ctl to use model 2 (USE SCRIPT 01_mod_ctl.sh)

model = 2 * 2: Empirical
aaRatefile = wag.dat
fix_alpha = 0
alpha = .5
ncatG = 4


Then we can run for each partition CODEML (use script 02_run_codeml.sh). This will generate an rst2 file for each one (run in different folders).
Then, concatenate (in the appropriate order, use script 03_catr_st2.sh) all the rst2 files into out.BV file (result in 01_usedata3/10partitions/out.BV). Use this one to run MCMCTree with
usedata = 2 in the final step.
