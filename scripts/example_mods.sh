export n_jobs=1
export python=3.7.0
export gcc=9.1.0
export cmake=3.15.3
export cuda=10.1.243

module purge
module load python/$python
module load gcc/$gcc
module load cmake/$cmake
module load cuda/$cuda
