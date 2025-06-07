# ScenarioReduction
This repository contains the code for the proposed scenario reduction methodology. The methodology is demonstrated using case studies from the associated paper. These case studies are implemented in the main script scenario_reduction_casestudies_milanvercammen.jl. This script serves as the main entry point and should be used to reproduce the results presented in the paper. All other function files and data must be placed in the same folder for the script to run correctly.

Repository Structure:

- scenario_reduction_casestudies_milanvercammen.jl: Main script that reproduces the case studies from the paper.

- scenario_reduction_functions.jl: Contains the core functions implementing the scenario reduction methodology.

- Other .jl files: Supporting scripts such as scenario generation (load_scenario_generation.jl), PMD distribution (CreatePMDDictionary.jl), etc. These are based on existing tools and methods as referenced in the corresponding files.

Data:

- Distribution network (DN):
The DN used in this study is the Spanish POLA network, represented by four JSON files (65019_74478_Mod_branches.json, 65019_74478_Mod_buses.json, 65019_74478_Mod_devices.json, 65019_74478_mod_configuration.json).
Additionally, the file matching_POLAnetwork_fluviusconsumers.csv provides the mapping between POLA network nodes and Fluvius consumer profiles.

- Load profiles:
Type 1 consumer profiles from Fluvius must be downloaded by the user from the following location:
https://github.com/umar-hashmi/Public-Load-profile-Datasets/blob/67649f157b0b1bc96a64f57504a0847c38bd3f85/Fluvius_Consumer_profiles_2022/Type1_consumer_Fluvius.zip
Instructions on where to place these files are provided at the top of the scenario_reduction_casestudies_milanvercammen.jl script.
Additionally, the file consumertype_IDs.csv contains the IDs of the Type 1 consumer profiles used in the case studies.
