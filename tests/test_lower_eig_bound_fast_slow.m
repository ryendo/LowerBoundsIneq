% p  = profile('info');
% FT = p.FunctionTable;

% fn   = string({FT.FunctionName})';
% tot  = [FT.TotalTime]';
% calls = [FT.NumCalls]';

% % 子関数に使った合計時間を引いて SelfTime を作る
% childTot = zeros(numel(FT),1);
% if isfield(FT,'Children')
%     for k = 1:numel(FT)
%         ch = FT(k).Children;
%         if ~isempty(ch) && isfield(ch,'TotalTime')
%             childTot(k) = sum([ch.TotalTime]);
%         else
%             childTot(k) = 0;
%         end
%     end
% end
% self = tot - childTot;
% self(self < 0) = 0;  % 数値誤差ガード

% T = table(fn, tot, self, calls, 'VariableNames', ...
%     {'FunctionName','TotalTime','SelfTime','NumCalls'});

% % TotalTime が重い順
% T1 = sortrows(T,'TotalTime','descend');
% disp(T1(1:min(20,height(T1)), :));

% % SelfTime（その関数の中身）が重い順
% T2 = sortrows(T,'SelfTime','descend');
% disp(T2(1:min(20,height(T2)), :));




format long infsup

nEig = 1;
lagrange_order = 3;

meshCR = make_mesh_by_gmsh(intval('0.5'), sqrt(intval('3'))/2, 0.1);
meshCG = make_mesh_by_gmsh(intval('0.5'), sqrt(intval('3'))/2, 0.2);


tic
eig_bounds = lower_eig_bound(lagrange_order, meshCR, meshCG, nEig)
toc


tic
eig_bounds = lower_eig_bound_fast(lagrange_order, meshCR, meshCG, nEig)
toc


vert = I_intval(meshCG.nodes);
edge = meshCG.edges;
tri  = meshCG.elements;
bd   = meshCG.boundary_edges;

tic
[eig_value, eig_func_no_bdry, eig_func_with_bdry, A_grad, A_L2, A_xx, A_xy, A_yy, bd_dof_idx] = laplace_eig_lagrange_detailed(lagrange_order, vert, edge, tri, bd, nEig);
toc

I_sup(eig_value)'

