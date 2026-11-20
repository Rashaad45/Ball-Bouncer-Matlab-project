%% RUN_ALL  One-click reproduction of every result in the report.
%  Typical run time: ~1-2 min in MATLAB R2021a+, ~10 min in GNU Octave 8+.
%  Works without any toolbox (no Control/Robotics/Optimization toolbox).

test_3rrs;                 % 1. self-tests (stops here if the maths is broken)
run_kinematics_analysis;   % 2. workspace, parasitic motion, IK-error maps
run_closed_loop;           % 3. 4 scenarios x 4 configurations
run_latency_study;         % 4. ablation vs camera latency
run_monte_carlo;           % 5. robustness to tolerances / servo / camera

% 6. animation of the proposed controller on the figure-8
S  = load('results/closed_loop_runs.mat');
P  = params_3rrs();
animate_3rrs(S.runs{3,4}, P, struct('file', 'results/figure8_proposed.mp4'));
