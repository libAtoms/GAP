GAP_F90_sources = [
  'clustering.F90',
  'descriptors.F90',
  'descriptors_wrapper.F90',
  'find_water_triplets_noncommercial.F90',
  'gp_fit.F90',
  'gp_predict.F90',
  'make_permutations_noncommercial_v2.F90',
  'soap_turbo.f90',
  'soap_turbo_angular.f90',
  'soap_turbo_compress.f90',
  'soap_turbo_functions.f90',
  'soap_turbo_radial.f90',
]

GAP = library('GAP',
  GAP_F90_sources,
  dependencies: [
    blas_dep,
    mpi_dep,
    ],
  link_with : [libAtoms,fox],
  link_args: ['-lgomp'],
  )
