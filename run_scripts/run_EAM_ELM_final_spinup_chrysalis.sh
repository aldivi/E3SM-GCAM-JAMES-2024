# ----- Directory paths -----
export USERID=ac.eva.sinha
export BASE_DIR=/home/${USERID}/
export E3SM_DIR=${BASE_DIR}/iesm
export E3SM_CASE_DIR=${E3SM_DIR}/cime/scripts
export E3SM_OUTPUT_DIR=/lcrc/group/e3sm/${USER}
export SURFACE_DATA_DIR=/lcrc/group/e3sm/data/inputdata/lnd/clm2/surfdata_map

# ------ Create new case -----
export RES=ne30pg2_f09_oEC60to30v3 # non-default grids are: atm:ne30np4.pg2  lnd:0.9x1.25  ocnice:oEC60to30v3  rof:null  glc:null  wav:null   mask is: oEC60to30v3
export COMPSET=1850_EAM%CMIP6_ELM%CNPRDCTCBC_CICE%PRES_DOCN%DOM_SROF_SGLC_SWAV
export COMPSET_alias=I1850EAMELMCNPRDCTCBC
export CASEID=20240102 # $(date '+%Y%m%d')
export CASE_NAME=${CASEID}_${COMPSET_alias}_${RES}

export CASE_SCRIPTS_DIR=${E3SM_OUTPUT_DIR}/${CASE_NAME}/case_scripts

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

export finidat_case=${CASE_NAME}_ad_spinup
export finidat=${E3SM_OUTPUT_DIR}/${finidat_case}/run/${finidat_case}.elm.r.0201-01-01-00000.nc
export fsurdat=${SURFACE_DATA_DIR}/surfdata_0.9x1.25_HIST_simyr1850_c201021.nc

export domainpath=/lcrc/group/e3sm/data/inputdata/share/domains
export lnd_domainfile=domain.lnd.0.9x1.25_oEC60to30v3.231108.nc
export atm_domainfile=domain.lnd.ne30pg2_oEC60to30v3.200220.nc

# suplphos = 'ALL' sets supplemental phosphorus as active for all vegetation types

cat >> user_nl_elm << EOF
&elm_inparm
 hist_mfilt = 1
 hist_nhtfrq = 0
 finidat = '$finidat'
 fsurdat = '$fsurdat'
 nyears_ad_carbon_only = 25
 spinup_mortality_factor = 10
 do_budgets = .false.
EOF

# ----- Case setup -----
./xmlchange SAVE_TIMING=FALSE
./xmlchange RUN_REFDATE=0201-01-01
./xmlchange ATM_DOMAIN_PATH=${domainpath}
./xmlchange LND_DOMAIN_PATH=${domainpath}
./xmlchange ATM_DOMAIN_FILE=${atm_domainfile}
./xmlchange LND_DOMAIN_FILE=${lnd_domainfile}
./xmlchange ATM_NCPL=48
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=80
./xmlchange REST_N=10
./xmlchange RESUBMIT=5
./xmlchange JOB_QUEUE=slurm
./xmlchange JOB_WALLCLOCK_TIME=48:00:00
./xmlchange MAX_MPITASKS_PER_NODE=64
./xmlchange MAX_TASKS_PER_NODE=64
./xmlchange NTASKS_ATM=5120
./xmlchange NTASKS_CPL=5120
./xmlchange NTASKS_OCN=5120
./xmlchange NTASKS_ICE=5120
./xmlchange NTASKS_LND=5120

./xmlchange ROOTPE=0
./xmlchange NTHRDS=1

./case.setup

# ----- Case build -----
./case.build

# ----- Run model -----
./case.submit
