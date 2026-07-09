# For each sampled complex, the following commands were used to explore the stoichiometric space:

# Initialization
multimer_mapper database.fasta --out_path mm_1

# Dimers
multimer_mapper database.fasta --AF_2mers 2mers/ --skip_plots --first_plot --skip_traj --auto_domains --out_path mm_2

# Nmers
multimer_mapper database.fasta --AF_2mers 2mers/ --AF_Nmers Nmers/ --skip_plots --first_plot --skip_traj --auto_domains --out_path mm_N

# Nmers with maximum combinatorial increments of two
multimer_mapper database.fasta --AF_2mers 2mers/ --AF_Nmers Nmers/ --skip_plots --first_plot --skip_traj --auto_domains --max_comb_order 2 --out_path mm_N 

# Once convergence has been reached, use the full method
multimer_mapper database.fasta --AF_2mers 2mers/ --AF_Nmers Nmers/ --skip_plots --first_plot --auto_domains --out_path mm_N
