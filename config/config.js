{
    "general": {
        "output_global_attributes": {
            "run_name": "Demo dvmdostem run",
            "description": "A basic test run. Data is not scientifically sound."
        }
    },
    "IO": {
        "parameter_dir": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/parameters/",
        "hist_climate_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/historic-climate.nc",
        "proj_climate_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/projected-climate.nc",
        "veg_class_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/vegetation.nc",
        "drainage_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/drainage.nc",
        "soil_texture_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/soil-texture.nc",
        "co2_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/co2.nc",
        "proj_co2_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/projected-co2.nc",
        "topo_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/topo.nc",
        "fri_fire_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/fri-fire.nc",
        "hist_exp_fire_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/historic-explicit-fire.nc",
        "proj_exp_fire_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/projected-explicit-fire.nc",
        "fire_fuel_model_file": "/home/ext_ejafarov_woodwellclimate_org/fire/dvm-dos-tem/external/Fireweed/Input_Files/Standard_Fuel_Models_Table.csv",
        "fire_fosberg_a_file": "/home/ext_ejafarov_woodwellclimate_org/fire/dvm-dos-tem/external/Fireweed/Input_Files/NWGC_Fosberg_Table_A.txt",
        "fire_fosberg_b_file": "/home/ext_ejafarov_woodwellclimate_org/fire/dvm-dos-tem/external/Fireweed/Input_Files/NWGC_Fosberg_Table_B.txt",
        "fire_fosberg_c_file": "/home/ext_ejafarov_woodwellclimate_org/fire/dvm-dos-tem/external/Fireweed/Input_Files/NWGC_Fosberg_Table_C.txt",
        "fire_fosberg_d_file": "/home/ext_ejafarov_woodwellclimate_org/fire/dvm-dos-tem/external/Fireweed/Input_Files/NWGC_Fosberg_Table_D.txt",
        "runmask_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/input/run-mask.nc",
        "output_dir": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/output/",
        "restart_from": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/output/restart-sp.nc",
        "output_spec_file": "/mnt/exacloud/ext_ejafarov_woodwellclimate_org/Exp_bgc_noWetland_split/batch_0/config/output_spec.csv",
        "output_monthly": 1,
        "output_nc_eq": 0,
        "output_nc_sp": 0,
        "output_nc_tr": 1,
        "output_nc_sc": 1,
        "output_interval": 1
    },
    "calibration-IO": {
        "unique_pid_tag": "",
        "caldata_tree_loc": "/tmp/batch_0"
    },
    "stage_settings": {
        "inter_stage_pause": false,
        "pr": {
            "env": true,
            "bgc": false,
            "nfeed": false,
            "avlnflg": false,
            "baseline": false,
            "dsb": false,
            "dsl": false,
            "dyn_lai": false
        },
        "eq": {
            "env": true,
            "bgc": true,
            "nfeed": true,
            "avlnflg": true,
            "baseline": true,
            "dsb": false,
            "dsl": true,
            "fire": false,
            "dyn_lai": true
        },
        "sp": {
            "env": true,
            "bgc": true,
            "nfeed": true,
            "avlnflg": true,
            "baseline": true,
            "dsb": false,
            "dsl": true,
            "fire": false,
            "dyn_lai": true
        },
        "tr": {
            "env": true,
            "bgc": true,
            "nfeed": true,
            "avlnflg": true,
            "baseline": true,
            "dsb": false,
            "dsl": true,
            "fire": true,
            "dyn_lai": true
        },
        "sc": {
            "env": true,
            "bgc": true,
            "nfeed": true,
            "avlnflg": true,
            "baseline": true,
            "dsb": false,
            "dsl": true,
            "fire": true,
            "dyn_lai": true
        }
    },
    "module_settings": {
        "dsb": {
            "fire": {
                "ignition_tr": 0,
                "ignition_sc": 0,
                "fire_process_wildfire": true,
                "fire_max_snow": 1.0,
                "fire_dynamic_fuel": true,
                "fire_moss_as_dead_fuel": true,
                "fire_calculate_delta": false,
                "fire_heat_frac_to_soil": 0.25,
                "fire_gf_heat_loss_factor": 0.5,
                "fire_gf_surface_pd": 1.0,
                "fire_gf_smolder_pd": 3.0,
                "fire_gf_layer_thickness": 0.1,
                "fire_temp_fm": -1,
                "fire_windspeed": 650.0,
                "fire_tempair": -300.0,
                "fire_rh_pct": -1.0,
                "fire_dwd_combust_frac": 1.0
            }
        }
    },
    "model_settings": {
        "cell_timelimit": 0,
        "dynamic_lai": 1,
        "baseline_start": 1901,
        "baseline_end": 1931
    }
}
