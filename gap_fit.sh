#!/bin/bash

gap_fit at_file=./Final_Training_Dataset.xyz \
gap={{distance_Nb order=2 cutoff=3.7 covariance_type=ARD_SE theta_uniform=1 n_sparse=15 delta=0.5 add_species=F sparse_method=FILE sparse_file=1.input Z={29 29}}: \
	 {distance_Nb order=2 cutoff=3.7 covariance_type=ARD_SE theta_uniform=1 n_sparse=15 delta=0.5 add_species=F sparse_method=FILE sparse_file=2.input Z={6 29}}: \
	 {distance_Nb order=2 cutoff=3.7 covariance_type=ARD_SE theta_uniform=1 n_sparse=15 delta=0.5 add_species=F sparse_method=FILE sparse_file=3.input Z={6 6}}: \
	 {distance_Nb order=3 cutoff=2.5 covariance_type=ARD_SE theta_uniform=1.0 delta=0.05 n_sparse=200 add_species=F sparse_method=FILE sparse_file=4.input Z={29 29 29 }}: \
	 {distance_Nb order=3 cutoff=2.5 covariance_type=ARD_SE theta_uniform=1.0 delta=0.05 n_sparse=200 add_species=F sparse_method=FILE sparse_file=5.input Z={6 29 29 }}: \
	 {distance_Nb order=3 cutoff=2.5 covariance_type=ARD_SE theta_uniform=1.0 delta=0.05 n_sparse=200 add_species=F sparse_method=FILE sparse_file=6.input Z={6 6 29 }}: \
	 {distance_Nb order=3 cutoff=2.5 covariance_type=ARD_SE theta_uniform=1.0 delta=0.05 n_sparse=200 add_species=F sparse_method=FILE sparse_file=7.input Z={6 6 6}}:\
	 {soap cutoff=3.7 covariance_type=dot_product cutoff_transition_width=1.0 central_weight=1.0 n_sparse=9000 zeta=4 delta=0.2 atom_sigma=0.5 l_max=4 n_max=12 add_species=F sparse_method=FILE sparse_file=8.input n_species=2 Z=29 species_Z={29 6 }}:\
     {soap cutoff=3.7 covariance_type=dot_product cutoff_transition_width=1.0 central_weight=1.0 n_sparse=9000 zeta=4 delta=0.2 atom_sigma=0.5 l_max=4 n_max=12 add_species=F sparse_method=FILE sparse_file=9.input n_species=2 Z=6 species_Z={29 6 }}} \
e0={C:0:Cu:0} \
default_sigma={0.001 0.01 0.05 0.0} \
config_type_sigma={Liquid:0.050:0.5:0.5:0.0: \
                   Liquid_Interface:0.050:0.5:0.5:0.0:\
				   Amorphous_Bulk:0.005:0.2:0.2:0.0: \
				   Amorphous_Surfaces:0.005:0.2:0.2:0.0: \
				   amoc-Cu:0.005:0.2:0.2:0.0: \
				   Cu_diffusion:0.002:0.1:0.2:0.0: \
				   metal_surface:0.002:0.1:0.2:0.0: \
				   CuC_crystal:0.001:0.01:0.05:0.0: \
				   Surfaces:0.002:0.1:0.2:0.0: \
				   Dimer:0.002:0.1:0.2:0.0: \
				   Fullerenes:0.002:0.1:0.2:0.0: \
				   Defects:0.001:0.01:0.05:0.0: \
				   Crystalline_Bulk:0.001:0.01:0.05:0.0: \
				   Nanotubes:0.001:0.01:0.05:0.0: \
				   Graphite:0.001:0.01:0.05:0.0: \
				   Diamond:0.001:0.01:0.05:0.0: \
				   Graphene:0.001:0.01:0.05:0.0: \
				   Graphite_Layer_Sep:0.001:0.01:0.05:0.0: \
				   Single_Atom:0.0001:0.001:0.05:0.0: \
				   cluster: 0.0001:0.001:0.05:0.0} \
energy_parameter_name=energy force_parameter_name=force virial_parameter_name=virial \
sparse_jitter=1.0e-8 \
do_copy_at_file=F sparse_separate_file=T \
mpi_blocksize=101 \
gp_file=./C-Cu.xml