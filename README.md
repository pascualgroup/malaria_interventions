# Organization
Experiments are organized in 4 levels:

1. **Parameter space**: This is the set of parameters such as biting rates, interventions, etc.
2. **Scenario**: Each parameter space can be run for the 3 scenarios.
3. **Experiment**: Experiments are variations on the parameter space. By definition '00' denotes the basic parameter space that is used to reach a steady-state and create a checkpoint. '01' is control. It has the *exact same parameters* as '00'. It loads the checkpoint created by 00 and continues the simulation. Any other number (e.g., 02, 03...) denotes a variation from 00 (e.g., to test interventions, different biting rates or different seasonal patterns).
4. **Runs**: Each combination of paramter space, scenario and experiment can be run multiple times. Importantly, the seed for the random number ()

# File names
Paramter files names are denoted as: PSxx_y_Ezz.py, where 
# Pipeline
Run im
