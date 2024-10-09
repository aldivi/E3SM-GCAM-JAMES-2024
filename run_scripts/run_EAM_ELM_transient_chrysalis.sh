# ----- Directory paths -----
export USERID=ac.eva.sinha
export BASE_DIR=/home/${USERID}/
export E3SM_DIR=${BASE_DIR}/iesm
export E3SM_CASE_DIR=${E3SM_DIR}/cime/scripts
export E3SM_OUTPUT_DIR=/lcrc/group/e3sm/${USER}
export SURFACE_DATA_DIR=/lcrc/group/e3sm/data/inputdata/lnd/clm2/surfdata_map

# ------ Create new case -----
export RES=ne30pg2_f09_oEC60to30v3 # non-default grids are: atm:ne30np4.pg2  lnd:0.9x1.25  ocnice:oEC60to30v3  rof:null  glc:null  wav:null   mask is: oEC60to30v3
export COMPSET=20TR_EAM%CMIP6_ELM%CNPRDCTCBC_CICE%PRES_DOCN%DOM_SROF_SGLC_SWAV
export COMPSET_alias=I20TREAMELMCNPRDCTCBC
export CASEID=$(date '+%Y%m%d')
export CASE_NAME=${CASEID}_${COMPSET_alias}_${RES}
export CASE_ARCHIVE_DIR=${E3SM_OUTPUT_DIR}/${CASE_NAME}/archive

export CASE_SCRIPTS_DIR=${E3SM_OUTPUT_DIR}/${CASE_NAME}/case_scripts

# select based on the model resolution
# 2 deg : surfdata_iESM_dyn_hist_simyr2015_c230516.nc
# 1 deg : landuse.timeseries_0.9x1.25_HIST_simyr2015_c201021.nc
#export iesm_dyn_source=surfdata_iESM_dyn_hist_simyr2015_c230516.nc
export iesm_dyn_source=landuse.timeseries_0.9x1.25_HIST_simyr2015_c201021.nc

# Delete old case and run directory
rm -rf ${CASE_SCRIPTS_DIR}
rm -rf ${E3SM_OUTPUT_DIR}/${CASE_NAME}

cd ${E3SM_CASE_DIR}
./create_newcase \
 --case ${CASE_NAME} \
 --compset ${COMPSET} \
 --res ${RES} \
 --output-root ${E3SM_OUTPUT_DIR} \
 --script-root ${CASE_SCRIPTS_DIR} \
 --project e3sm \
 --machine chrysalis

# ----- Modify user_nl_elm -----
cd ${CASE_SCRIPTS_DIR}

export finidat_COMPSET_alias=I1850EAMELMCNPRDCTCBC
export finidat_CASEID=20240102
export finidat_case=${finidat_CASEID}_${finidat_COMPSET_alias}_${RES}
export RUN_REFDATE=0401-01-01
export finidat=${E3SM_OUTPUT_DIR}/${finidat_case}/run/${finidat_case}.elm.r.${RUN_REFDATE}-00000.nc
export fsurdat=${SURFACE_DATA_DIR}/surfdata_0.9x1.25_HIST_simyr1850_c201021.nc
export flanduse=${SURFACE_DATA_DIR}/landuse.timeseries_0.9x1.25_HIST_simyr1850-2015_c201021.nc

export domainpath=/lcrc/group/e3sm/data/inputdata/share/domains
export lnd_domainfile=domain.lnd.0.9x1.25_oEC60to30v3.231108.nc
export atm_domainfile=domain.lnd.ne30pg2_oEC60to30v3.200220.nc

# suplphos = 'ALL' sets supplemental phosphorus as active for all vegetation types

cat >> user_nl_elm << EOF
&elm_inparm
 hist_mfilt = 1, 1
 hist_nhtfrq = 0, 0
 hist_dov2xy = .true., .false.
 hist_fincl2 = 'GPP', 'ER', 'HR', 'NPP'
 fsurdat = '$fsurdat'
 flanduse_timeseries = '$flanduse'
 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 do_budgets = .true.
 model_year_align_ndep = 1850
 stream_year_first_ndep = 1850
 stream_year_last_ndep = 2005
 model_year_align_pdep = 2000
 stream_year_first_pdep = 2000
 stream_year_last_pdep = 2000
 model_year_align_popdens = 1850
 stream_year_first_popdens = 1850
 stream_year_last_popdens = 2010
EOF

# ----- Case setup -----
./xmlchange SAVE_TIMING=FALSE
./xmlchange RUN_TYPE=hybrid
./xmlchange GET_REFCASE=TRUE
./xmlchange RUN_REFDIR=/lcrc/group/e3sm/ac.eva.sinha/E3SM_GCAM_lnd_init/${finidat_case}
./xmlchange RUN_REFCASE=${finidat_case}
./xmlchange RUN_REFDATE=${RUN_REFDATE}
./xmlchange RUN_STARTDATE=1850-01-01
./xmlchange ATM_DOMAIN_PATH=${domainpath}
./xmlchange LND_DOMAIN_PATH=${domainpath}
./xmlchange ATM_DOMAIN_FILE=${atm_domainfile}
./xmlchange LND_DOMAIN_FILE=${lnd_domainfile}
./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}
./xmlchange ATM_NCPL=48
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=20
./xmlchange REST_N=5
./xmlchange RESUBMIT=2
./xmlchange JOB_QUEUE=slurm
./xmlchange JOB_WALLCLOCK_TIME=12:00:00
./xmlchange MAX_MPITASKS_PER_NODE=64
./xmlchange MAX_TASKS_PER_NODE=64
./xmlchange NTASKS_ATM=5120
./xmlchange NTASKS_CPL=5120
./xmlchange NTASKS_OCN=5120
./xmlchange NTASKS_ICE=5120
./xmlchange NTASKS_LND=5120
./xmlchange SSTICE_YEAR_ALIGN=1869
./xmlchange SSTICE_YEAR_START=1869
./xmlchange SSTICE_YEAR_END=1879

./xmlchange ROOTPE=0
./xmlchange NTHRDS=1

./case.setup

# ----- Case build -----
./case.build

# ----- Run model -----
./case.submit
