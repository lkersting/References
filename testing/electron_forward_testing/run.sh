#!/bin/sh
# This file is named submit-script.sh
#SBATCH --partition=univ2                # Use the newer infrastructure
#SBATCH --time=0-12:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=4000

# ------------------------------- COMMANDS ------------------------------------

# -- DAGMC Tests --

cd /home/ecmartin3/tests/hpc_scaling/dagmc

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_1kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_1.txt 2>&1

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_10kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_2.txt 2>&1

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_100kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_4.txt 2>&1


# -- ROOT Tests --

cd /home/ecmartin3/tests/hpc_scaling/root

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_1kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_24.txt 2>&1

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_10kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_32.txt 2>&1

mpiexec -n 25 ./facemc-mpi --sim_info="sim_info.xml" --geom_def="h_sphere_geom.xml" --src_def="h_sphere_source.xml" --resp_def="h_sphere_rsp_fn.xml" --est_def="h_sphere_est_100kev.xml" --mat_def="h_sphere_mat.xml" --cross_sec_dir="/home/software/mcnpdata" > MPI_Strong_48.txt 2>&1
