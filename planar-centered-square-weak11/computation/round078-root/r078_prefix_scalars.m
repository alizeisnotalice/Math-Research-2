% Exact symbolic identities and one whole-cell witness, not halo integration.
tic;
syms n x y real
h=(n+1)/n^2; e=1/n^2; k=x+2; l=y+2;
exprs=[expand(4*(k-1)*(l-1)-k*l); simplify(h^2-e); expand(4*(k-1)^2-k)];
expected=[3*x*y+2*x+2*y; (2*n+1)/n^4; 4*x^2+7*x+2];
assert(all(isAlways(simplify(exprs-expected)==0)));
N=20; E=sym(1)/400; H=sym(21)/400; b=sym(37)/100; B=2*b;
lo=-sym(61)/1000; hi=-sym(59)/1000; r=sym(89)/200;
[I,J]=ndgrid(0:N-1,0:N-1); px=sym(I(:))*H; py=sym(J(:))*H;
% All old coordinates are nonnegative: farthest E corner is lo, nearest hi.
allOld=isAlways(px-lo<r & py-lo<r);
possibleOld=isAlways(px-hi<r & py-hi<r);
assert(isequal(allOld,possibleOld)); assert(sum(allOld)==64);
assert(all(allOld==logical(I(:)<8 & J(:)<8)));
outsideAllDist=max([sym(1)/2+hi,-lo]);
assert(isAlways(outsideAllDist<r));
fullmass=64*E+B; fullslack=simplify(fullmass-4*r^2);
assert(isAlways(fullslack==sym(1079)/10000));
minOutsideRadius=sym(439)/1000;
thinMax=B+10*E; singleMax=b+100*E;
gaps=simplify(4*minOutsideRadius^2-[B,thinMax,singleMax]);
assert(all(isAlways(gaps>0)));
assert(isAlways(-hi>1/(2*sym(N))));
cellArea=simplify((hi-lo)^2);assert(isAlways(cellArea==sym(1)/250000));
report=struct;
report.scope='Exact scalar identities and whole-cell fixed-input witness only; no 440 caps, menu union, numerical search or general C2 bound';
report.matlab=version; v=ver('symbolic');report.symbolic=v.Version;
report.pid=feature('getpid');report.identity_lhs=arrayfun(@(z)char(z),exprs,'UniformOutput',false);
report.identity_expected=arrayfun(@(z)char(z),expected,'UniformOutput',false);
report.cell={char(lo),char(hi),char(lo),char(hi)};report.radius=char(r);
report.always_captured_old_count=sum(allOld);report.possibly_captured_old_count=sum(possibleOld);
report.common_radius_outside_max_distance=char(outsideAllDist);report.full_mass=char(fullmass);report.full_slack=char(fullslack);
report.minimum_any_outside_radius=char(minOutsideRadius);report.residual_thin_single_squared_gaps=arrayfun(@(z)char(z),gaps,'UniformOutput',false);
report.cell_area=char(cellArea);report.pass=true;report.seconds=toc;
file=fullfile(fileparts(mfilename('fullpath')),'matlab-prefix-scalars.json');fid=fopen(file,'w');assert(fid>0);fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(jsonencode(report));
