# 3-RRS Ball Balancer — Inverse Kinematics & Latency-Aware Control (MATLAB)

Exact kinematics of a 3-servo (3-RRS) ball-balancing parallel robot, a nonlinear
ball-on-moving-plate simulation, and a vision-in-the-loop controller
(LQR + integral, flatness feed-forward, delay-compensated Kalman predictor,
kinematics-aware governor). No toolboxes required.

## Quick start
```matlab
cd BallBalancer3RRS
test_3rrs        % 7 self-tests, ~5 s
run_all          % everything: kinematics, closed loop, latency study, Monte Carlo, animation
```
MATLAB R2019b+ or GNU Octave 8+. In MATLAB `run_all` also writes `results/figure8_proposed.mp4`.

## Headline results (simulation, 60 FPS camera, 33 ms latency)
| Scenario | Standard (exact IK, naive KF) | Proposed |
|---|---|---|
| Circle r=50 mm | 9.07 mm RMS | **1.23 mm** |
| Figure-8 | 8.49 mm | **1.30 mm** |
| Fast 120 mm move | ball falls off | **8.25 mm** |
| Latency 100 ms (circle) | ball falls off | **1.97 mm** |
| Monte Carlo, 40 perturbed plants (median) | 9.26 mm | **3.59 mm** |

Exact vs shortcut IK changes closed-loop error by < 1 %; shortcut IKs carry 0.5–1.1° static
tilt error, and the fixed-centre IK used in several open-source builds breaks plate rigidity by up to 0.9 %.

## Files
| File | Purpose |
|---|---|
| `params_3rrs.m` | all parameters — edit this for your CAD |
| `ik_3rrs.m`, `parasitic_3rrs.m`, `fk_3rrs.m`, `jacobian_3rrs.m`, `rot_3rrs.m` | exact kinematics |
| `ik_linear_3rrs.m`, `ik_fixedcentre_3rrs.m` | baseline (shortcut) IKs for comparison |
| `workspace_3rrs.m` | reachable-tilt lookup table |
| `design_controller_3rrs.m`, `dlqr_riccati.m` | LQR + integral + Kalman design |
| `governor_3rrs.m` | workspace + servo-rate governor |
| `sim_3rrs.m` | closed-loop simulation (600 Hz plant / 60 Hz control) |
| `scenarios_3rrs.m`, `configs_3rrs.m` | benchmarks and configurations |
| `run_*.m`, `animate_3rrs.m`, `test_3rrs.m` | studies, animation, tests |
| `results/` | figures, CSVs, per-run logs (generated in Octave 8) |
| `refs/` | cloned open-source ball-balancer repos (code only, `.git` and media stripped) |

## Reference repositories (in `refs/`)
- https://github.com/I-M-Robotics-Lab/Ball-Balancing-Robot
- https://github.com/thomasonzhou/ball-balancing-software
- https://github.com/Martin-Ansteensen/ball-platform

Their licences apply to their code.
