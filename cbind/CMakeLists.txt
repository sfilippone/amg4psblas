set(AMG_amgcbind_source_files
  amgprec/amg_zprec_cbind_mod.F90
  amgprec/amg_prec_cbind_mod.F90
  amgprec/amg_dprec_cbind_mod.F90
)
foreach(file IN LISTS AMG_amgcbind_source_files)
  list(APPEND amgcbind_source_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

set(AMG_amgcbind_source_C_files
  amgprec/amg_c_zprec.c
  amgprec/amg_c_dprec.c
)


foreach(file IN LISTS AMG_amgcbind_source_C_files)
  list(APPEND amgcbind_source_C_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

set(AMG_amgcbind_header_C_files
#  amgprec/amg_const.h
  amgprec/amg_c_dprec.h
  amgprec/amg_cbind.h
  amgprec/amg_c_zprec.h
)


foreach(file IN LISTS AMG_amgcbind_header_C_files)
  list(APPEND amgcbind_header_C_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()
