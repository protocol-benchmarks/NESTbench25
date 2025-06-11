#ifndef ENGINE_ISING_H
#define ENGINE_ISING_H

#ifdef __cplusplus
extern "C" {
#endif

//functions exposed from engine_ising.c.

void final_answer(void); //100 instances of the order parameter calculated using 10^4 trajectories
void load_protocol(void);//reads input_control_parameters.dat
void output_protocol(void);
void output_optimal_protocol(void);
void visualize_protocol(void); //visualize the effect of the current protocol
void load_default_protocol(void); //loads default protocol
double calculate_order_parameter(int n_traj); //computes order parameter over n_traj trajectories


#ifdef __cplusplus
}
#endif

#endif // ENGINE_ISING_H
