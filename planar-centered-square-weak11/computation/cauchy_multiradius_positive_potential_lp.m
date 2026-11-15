function cauchy_multiradius_positive_potential_lp(output_root)
% Numerical convex LP proposal. Cell objective is a Jensen lower approximation.
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
started=tic;radii=[1.125,1.25,1.5,2,3,4,6,8,12,16];angles=[0,1/8,1/4,3/8,1/2];orbits=struct([]);J=0;
for R=radii
 for t=angles
  a=R*(1-t);b=R*t;pts=[];
  for sx=[-1,1],for sy=[-1,1],pts=[pts;sx*a,sy*b;sx*b,sy*a];end,end %#ok<AGROW>
  J=J+1;orbits(J).R=R;orbits(J).t=t;orbits(J).points=unique(pts,'rows');
 end
end
L=max(radii);all_coords=[];for j=1:J,all_coords=[all_coords;abs(orbits(j).points(:))];end %#ok<AGROW>
y_nodes=unique([0:1/16:4,4.25:.25:8,8.5:.5:L,all_coords']);
x_arm=unique([L:L/16:2*L,2*L+2*L/16:2*L/16:4*L]);T=x_arm(end);ny=numel(y_nodes)-1;
[X,Y]=ndgrid(y_nodes,y_nodes);[XA,YA]=ndgrid(x_arm,y_nodes);
ncentral=ny^2;narm=(numel(x_arm)-1)*ny+ny;ncell=ncentral+narm;
C=zeros(ncell,J);
for j=1:J
 pts=orbits(j).points;P=primitive(X,Y);PA=primitive(XA,YA);M=-signed_xlogx(y_nodes);
 for k=1:size(pts,1)
  P=P-primitive(X-pts(k,1),Y-pts(k,2))/size(pts,1);
  PA=PA-primitive(XA-pts(k,1),YA-pts(k,2))/size(pts,1);
  M=M+signed_xlogx(y_nodes-pts(k,2))/size(pts,1);
 end
 central=diff(diff(P,1,1),1,2);arm=diff(diff(PA,1,1),1,2);
 tail=diff(M-PA(end,:));
 C(:,j)=[central(:);arm(:);tail(:)];
end
half_multiplicity=[2*ones(ncentral,1);4*ones(narm,1)];
[dx,dy]=diamond_points(32,128);WD=zeros(numel(dx),J);
for j=1:J,WD(:,j)=orbit_potential(dx,dy,orbits(j).points);end
Amat=[sparse(-WD),sparse(numel(dx),ncell);sparse(C),-speye(ncell)];
bvec=[-ones(numel(dx),1);zeros(ncell,1)];f=[zeros(J,1);half_multiplicity];
opts=optimoptions('linprog','Display','iter','ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8);
[v,fval,exitflag,lpout]=linprog(f,Amat,bvec,[],[],zeros(J+ncell,1),[],opts);
assert(exitflag>0,'LP did not converge');a=max(v(1:J),0);
[fx,fy]=diamond_points(128,512);WF=zeros(numel(fx),J);
for j=1:J,WF(:,j)=orbit_potential(fx,fy,orbits(j).points);end
fine_min=min(WF*a);assert(fine_min>0,'No positive fine sample minimum');a=a/fine_min;
cell_cost=half_multiplicity'*max(C*a,0);
integrand=@(t,y)4*L*max(total_potential(y+L*t./(1-t),y,a,orbits),0)./(1-t).^2;
actual_cost=integral2(integrand,0,1,0,L,'Method','tiled','AbsTol',1e-5,'RelTol',1e-4);
for j=1:J,orbits(j).coefficient=a(j);end
out=struct('status','completed_numerical_multiradius_potential_LP','matlab_version',version,'workers',1,'radii',radii,'angles',angles,'orbit_count',J,'cell_count',ncell,'coarse_D_point_count',numel(dx),'fine_D_point_count',numel(fx),'y_nodes',y_nodes,'arm_x_nodes',x_arm,'unbounded_tail_cells',ny,'LP_exitflag',exitflag,'LP_output',lpout,'LP_objective_lower_approximation',fval,'fine_min_before_rescale',fine_min,'rescaled_cell_lower_approximation',cell_cost,'rescaled_numerical_full_cost',actual_cost,'total_source_mass',sum(a),'orbits',orbits,'integral_AbsTol',1e-5,'integral_RelTol',1e-4,'elapsed_seconds',toc(started),'limits','Floating convex LP and quadrature only. Coefficients nonnegative and each orbit compensated exactly in the model, so continuum C4 follows abstract interface, but full D domination/integral/LP optimality are not certified. No new C2 upper bound.');
fid=fopen(fullfile(output_root,'multiradius-positive-potential-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('LP=%.15g cell=%.15g integral=%.15g fineMin=%.15g\n',fval,cell_cost,actual_cost,fine_min);
end
function [x,y]=diamond_points(N,Nb)
x=[];y=[];
for i=0:N,for j=0:i,if i+j<=N&&i+j>0,x(end+1,1)=i/N;y(end+1,1)=j/N;end,end,end %#ok<AGROW>
t=(0:Nb)'/(2*Nb);x=[x;1-t];y=[y;t];pts=unique([x,y],'rows');x=pts(:,1);y=pts(:,2);
end
function out=total_potential(x,y,a,orbits)
out=zeros(size(x));for j=1:numel(a),if a(j)>0,out=out+a(j)*orbit_potential(x,y,orbits(j).points);end,end
end
function w=orbit_potential(x,y,pts)
r=x+y;w=zeros(size(r));for k=1:size(pts,1)
 d=delta(x,pts(k,1))+delta(y,pts(k,2));w=w+d./(r.*(r+d))/size(pts,1);
end
end
function d=delta(u,c)
if c>=0,d=c-2*min(u,c);else,d=(-c)+zeros(size(u));end
end
function h=primitive(x,y)
a=abs(x);b=abs(y);h=sign(x).*sign(y).*(xlogx(a+b)-xlogx(a)-xlogx(b));
end
function f=xlogx(x)
f=zeros(size(x));idx=x>0;f(idx)=x(idx).*log(x(idx));
end
function f=signed_xlogx(x)
f=zeros(size(x));idx=x~=0;f(idx)=x(idx).*log(abs(x(idx)));
end
