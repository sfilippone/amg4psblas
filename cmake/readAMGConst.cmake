
macro(_amg_read_const)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/amgprec/amg_base_prec_type.F90" _amg_const_mod)
string(REGEX MATCH "amg_version_major_[ \t]+=[ \t]+([0-9]+)" _amg_version_major_match "${_amg_const_mod}")
set(AMGMAJOR  "${CMAKE_MATCH_1}")
string(REGEX MATCH "amg_version_minor_[ \t]+=[ \t]+([0-9]+)" _amg_version_minor_match "${_amg_const_mod}")
set(AMGMINOR  "${CMAKE_MATCH_1}")
string(REGEX MATCH "amg_patchlevel_[ \t]+=[ \t]+([0-9]+)" _amg_patchlevel_match "${_amg_const_mod}")
set(AMGPATCH  "${CMAKE_MATCH_1}")
string(REGEX MATCH "amg_version_string_[ \t]+=[ \t]+([^ \t]+)" _amg_version_string_match "${_amg_const_mod}")
set(AMGSTRING  "${CMAKE_MATCH_1}")

endmacro(_amg_read_const)

