function C = configs_3rrs()
% CONFIGS_3RRS  Controller/kinematics configurations compared in the study.
%   B1/B2 reproduce what typical open-source builds do; P is the proposed
%   pipeline; A1/A2 are the ablation steps between them.
mk = @(n,ik,dc,gv,ff) struct('name',n,'ik',ik,'delay_comp',dc,'governor',gv,'ff',ff);
C.B1 = mk('B1 Linear IK',        'linear',      false, false, false);
C.B2 = mk('B2 Fixed-centre IK',  'fixedcentre', false, false, false);
C.A1 = mk('A1 Exact IK',         'exact',       false, false, false);
C.A2 = mk('A2 Exact+FF',         'exact',       false, false, true);
C.A3 = mk('A3 Exact+FF+DelayKF', 'exact',       true,  false, true);
C.P  = mk('P  Proposed (full)',  'exact',       true,  true,  true);
end
